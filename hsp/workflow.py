import os
import sys
import json
import xml.etree.ElementTree as ET
import time
from datetime import datetime, timedelta
import multiprocessing
import numpy as np
import cv2

from hsp.modules.config import orderjson_to_config
from hsp.utils.filewatcher import TimeoutFileWatcher
from hsp.utils import parallel, common
from hsp.modules import object_tracking, geometric_locating

from trajectory_evaluate import FJ_target, DD_target, map_values, read_json

def xml_to_dict(element : ET):
    result = {}
    
    if element.attrib:
        result = element.attrib
        
    if element.text and element.text.strip():
        result = element.text.strip()
        
    for child in element:
        child_data = xml_to_dict(child)

        if child.tag in result:
            if isinstance(result[child.tag], list):
                result[child.tag].append(child_data)
            else:
                result[child.tag] = [result[child.tag], child_data]
        else:
            result[child.tag] = child_data

    return result

def xml_string_to_json(xml_string):
    root = ET.fromstring(xml_string)
    data_dict = xml_to_dict(root)
    return json.dumps({root.tag:data_dict}, indent=2, ensure_ascii=False)

# def generate_object_detection_xml(json : dict):

# def json_to_tasks(json : str, taskdir : os.path):

def read_order_file(orderpath, findcfg = False):
    cfg = orderjson_to_config(orderpath)
    if cfg is None:
        return None
    os.makedirs(cfg['output_dir'], exist_ok=True)
    os.makedirs(cfg['temporary_dir'], exist_ok=True)

    configpath = os.path.join(cfg['temporary_dir'],'hsp.json')
    if findcfg and os.path.exists(configpath):
        with open(configpath, 'r') as f:
            cfg = json.load(f)
    else:
        with open(configpath, 'w') as f:
            f.write(json.dumps(cfg, indent=2, ensure_ascii=False))

    return cfg

def while_loop_with_events( process, total, batchsize, event, logger):
    progress = 0
    while progress < total:
        num = min(batchsize, total-progress)

        process(progress, num)

        progress += num
        logger.progress_update(progress)

        if event.is_terminated():
            return True, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return True, "用户终止"

def evaluation(parameters):
    with open(parameters['OrderPath'],'r') as f:
        data = json.load(f)
        eval = data['InterfaceFile']["Evaluation"]
        ref = read_json(eval['RealTrajectory'])
        target = read_json(eval['UnvalTrajectory'])
        type = eval['Type']
        choice = eval['Indicates']
        if type == "2":
            point_error, mean_speed_error, mean_deg_error = FJ_target(target, ref)
            mapping = {
                '0': point_error,
                '1': mean_speed_error,
                '2': mean_deg_error
            }
            results = map_values(choice, mapping)
            return results
        else:
            position_mean_error, mean_speed_error = DD_target(target, ref)
            mapping = {
                '0': position_mean_error,
                '1': mean_speed_error,
            }
            results = map_values(choice, mapping)
            return results

    return None

def background_extraction( output, directory, imagelist, config):
    if not config['overwritten'] and os.path.exists(output):
        return
    images = common.concatenate_images(imagelist, directory)
    bg = object_tracking.extract_background(images, config['background_method'])
    common.rasterio_write(output, bg)

    if config['debug']:
        with open(os.path.join(os.path.splitext(output)[0]+'.txt'), 'w') as f:
            f.write('\n'.join(imagelist))

def anomaly_detection( output_prefix, image, background, config):
    if not config['overwritten'] and os.path.exists(output_prefix+config['detection_output_postfix'][0]):
        return

    # define environment variables
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(config['omp_num_threads'])

    object_tracking.detect_anomaly( image, os.path.dirname(output_prefix), background, output_prefix, config['debug'])

def anomaly_tracking( frameid : int, trackers : list,  seedfiles : list, params : dict, load_file):
    trackers_finished = []
    for i, file in enumerate(seedfiles):
        seeds = load_file(file)
        trackers = object_tracking.pointwise_tracking(seeds[0], trackers, frameid+i, params['kalman'])
        for tracker in trackers:
            if not tracker.good :
                tracker.remove_all_missings()
                if tracker.is_valid(params['min_frame_number'], params['min_speed'], params['max_acceleration'], params.get('curvature','spline')):
                    trackers_finished.append(tracker)
        trackers = [tracker for tracker in trackers if tracker.good]
    return trackers_finished, trackers

def trajectory_file_refinement( trajectorydir, trajectoryname, framedir, config):
    trajectoryfile = os.path.join(trajectorydir, trajectoryname)
    start, cnt, points = object_tracking.load_tracking(trajectoryfile)
    trajectory = object_tracking.refine_trajectory(points, framedir, config)
    with open(trajectoryfile, 'w') as f:
        f.write('{} {}\n'.format(start, cnt))
        for pt in trajectory:
            f.write('{}\n'.format('\t'.join(f"{x:.2f}" if isinstance(x,float) else str(x) for x in pt)))


def trajectory_locating( output, trajectory, type, framedir, frametime, config):
    if not config['overwritten'] and os.path.exists(output):
        return
    with open(trajectory, 'r') as f:
        start, cnt, points = object_tracking.load_tracking(trajectory)
        
        if cnt <= 0:
            return 

        header_info = {
            "Category": {
                "Name": "DD" if type else "FJ",
                "Confidence": min(len(points)/100, 1),
                "Classification":{
                    "Name" : "",
                    "BoostStage": "",
                    "Confidence": 0
                }
            }
        }

        output_info = {
            "LaunchAzimuth": 0,
            "ID": os.path.splitext(os.path.basename(output))[0],
            "DataTime": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "SatelliteList": {"Satellite": [{"ID": "1"}]}, "PointList": {"Point": []}}

        sofa_transformer = config['sofa']
        points_info = output_info["PointList"]["Point"]
        locating = geometric_locating.Locating(os.path.join(framedir, points[0][0]))
        if type:
            height_range = 2*locating.height_scale()
            acceleration = min(2.5*9.8, 2*height_range/(cnt**2)) 
            z0 = locating.height_off()-locating.height_scale()
        else:
            acceleration = 0.0
            z0 = locating.height_off()
        for i,point in enumerate(points):
            frame = os.path.join(framedir, point[0])
            info = {
                "Time": frametime[start+i] if frametime is not None else datetime.fromtimestamp(os.path.getctime(frame)).strftime("%Y-%m-%d %H:%M:%S"),
                "ImageCoordinates": [point[1],point[2]],
                "DigitalNumber": point[3],
                "Energy": point[3]*0.01
                }

            if locating.load_rpc(frame):

                z = z0+acceleration*(i)**2/2

                graphic, proj, centric = locating.transform(point[1], point[2], z)

                info['Location'] = [graphic[0],graphic[1],z]
                info['Projection'] = [proj[0], proj[1], z]
                info['CGCS2000'] = [centric[0],centric[1],centric[2]]
            
                if sofa_transformer is not None:
                    cgcs2j = np.linalg.inv(sofa_transformer.j2000_to_cgcs2000_matrix(info['Time'])) 

                    j2000 = cgcs2j @ np.array(centric).reshape(3,1)
                    info['J2000'] = [j2000[0,0],j2000[1,0],j2000[2,0]]

            points_info.append(info)

        points_info[0]['WarningStatus'] = '01H'
        points_info[-1]['WarningStatus'] = '03H'
        for i in range(1,len(points_info)-1):
            points_info[i]['WarningStatus'] = '02H'

        proj = np.array([x['Projection'] for x in points_info])
        vx = np.gradient(proj, axis=0)

        for i in range(len(points_info)):
            points_info[i]['Velocity'] = list(vx[i])
            points_info[i]['Type'] = 'Observed'

        with open(output, 'w') as fout: 
            json.dump({"Header":header_info, "Trajectory":output_info}, fout, indent=2, ensure_ascii=False)

def trajectory_predicting(output, trajectory, config):
    if not config['overwritten'] and os.path.exists(output):
        return
    import hsp.modules.trajectory_prediction as tp
    import pyproj
    from hsp.rpcm import compute_epsg

    def generate_points_dict(points, time, geographic_to_proj, sofa):
        infos = []
        for point in points:
            info = {
                "Time": time.strftime("%Y-%m-%d %H:%M:%S"),
                "ImageCoordinates": [-1,-1],
                "DigitalNumber": 0,
                "Energy": 0
                }
            graphic = geometric_locating.transform_geocentric_to_geographic(*point)
            proj = geographic_to_proj.transform(*graphic)

            info['Location'] = list(graphic)
            info['Projection'] = list(proj)
            info['CGCS2000'] = list(point)

            if sofa is not None:
                cgcs2j = np.linalg.inv(sofa_transformer.j2000_to_cgcs2000_matrix(time)) 

                j2000 = cgcs2j @ np.array(point).reshape(3,1)
                info['J2000'] = [j2000[0,0],j2000[1,0],j2000[2,0]]

            info['WarningStatus'] = '02H'
            infos.append(info)

            time += timedelta(seconds=1)

        for i in range(len(infos)):
            infos[i]['Velocity'] = [0,0,0]
            infos[i]['Type'] = 'Predicted'
        return infos

    with open(trajectory, 'r') as f:

        try:
            data = json.load(f)
            if data['Header']['Category']['Name'].find('DD') < 0:
                return
        except:
            return 

        points = [x['CGCS2000'] for x in data['Trajectory']['PointList']['Point']]
        time = [datetime.strptime(x['Time'], "%Y-%m-%d %H:%M:%S")
                 for x in data['Trajectory']['PointList']['Point']]
    
        if points is None or len(points) < 5:
            return 

        missile_model_params = [35.0, 70.0, 50.0]
        
        velocities = np.gradient(points, axis=0)
        accelerations = np.linalg.norm(np.gradient(velocities, axis=0), axis=1)

        is_two_stage, jump_index = tp.detect_acceleration_jump(accelerations)

        if is_two_stage:
            # 根据观测数据动态调整参数
            max_acc = np.max(accelerations)
            avg_acc = np.mean(accelerations)

            two_stage_params = {
                'a0_1': avg_acc * 0.6,  # 一级起始加速度
                'a1_1': avg_acc * 0.9,  # 一级结束加速度
                'T1': jump_index + 10,  # 基于检测到的跳跃点
                'a0_2': max_acc * 0.7,  # 二级起始加速度
                'a1_2': max_acc * 1.1,  # 二级结束加速度
                'T2': len(accelerations) - jump_index + 20,  # 二级助推时间
                'T_total': len(accelerations) + 30  # 总助推时间
            }

            # 匹配二级助推模型
            match_start_time, correlation, mse, stage_info = tp.find_boost_phase_position(
                accelerations, missile_model_params, 'two_stage', two_stage_params)

            # 使用二级助推参数进行后续计算
            effective_model_params = two_stage_params
            stage_type = 'two_stage'

        else:
            # 单级助推
            match_start_time, correlation, mse, stage_info = tp.find_boost_phase_position(
                accelerations, missile_model_params, 'single')
            effective_model_params = missile_model_params
            stage_type = 'single'

        stop_altitude = config['ground_altitude']
        launch_points = tp.estimate_launch_point_kinematic(
            points, velocities, match_start_time, stop_altitude,
            effective_model_params, stage_type)

        points_to_shutdown, shutdown_velocity_ecef = tp.estimate_shutdown_point_kinematic(
            points, velocities, accelerations, match_start_time,
            effective_model_params, stage_type)

        time_interval = 1
        landing_points = tp.missile_impact_prediction(
            points[-1] if len(points_to_shutdown) == 0 else points_to_shutdown[-1], shutdown_velocity_ecef, stop_altitude, time_interval)

        points_info = []
        sofa_transformer = config['sofa']
        lonlat = data['Trajectory']['PointList']['Point'][0]['Location']
        geographic_to_proj = pyproj.Transformer.from_crs(4326, compute_epsg(lonlat[0], lonlat[1]), always_xy=True)

        if len(launch_points) > 0:
            launch_points.reverse()
            infos = generate_points_dict(launch_points, time[0]-timedelta(seconds=len(launch_points)), geographic_to_proj, sofa_transformer)
            points_info += infos

        points_info += data['Trajectory']['PointList']['Point']

        if len(points_to_shutdown) > 0:
            infos = generate_points_dict(points_to_shutdown, time[-1], geographic_to_proj, sofa_transformer)
            points_info += infos

        points_info[-1]['Type'] += '|BurnOut'

        if len(landing_points)>0:
            infos = generate_points_dict(landing_points, time[-1]+timedelta(seconds=len(points_to_shutdown)), geographic_to_proj, sofa_transformer)
            points_info += infos

        proj = np.array([x['Projection'] for x in points_info])
        vxs = np.gradient(proj, axis=0)

        for vx,pt in zip(vxs,points_info):
            pt['Velocity'] = list(vx)

        data['Trajectory']['PointList']['Point'] = points_info
        data['Header']['Category']['Classification']['BoostStage'] = "2" if is_two_stage else "1"
        with open(output, 'w') as fout: 
            json.dump(data, fout, indent=2, ensure_ascii=False)

def automatic_processing(cfg, event, logger, share, start_from = 0):
    def newfile_callback(directory, filename):
        if any(filename.lower().endswith(ext) for ext in newfile_callback.extensions):
            newfile_callback.files.append(filename)
        return True

    common.print_elapsed_time.t0 = datetime.now()
    logger.set_step(0) 

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes'] is not None:
        nb_workers = cfg['max_processes']
    
    batchsize = cfg['background_frame_number']*nb_workers
    if cfg['batchsize'] is not None:
        batchsize = cfg['batchsize']

    dir_input = cfg['input_dir']
    dir_output = cfg['output_dir']
    extensions = cfg['input_file_extensions']
    wsdir = cfg['temporary_dir']

    if start_from <= 1:
        logger.set_step(1) 
        print('1) tracking objects in frames ...')
        dir_detection = os.path.join(wsdir, object_tracking.DETECTION_DIRNAME)
        os.makedirs(dir_detection, exist_ok=True)
        dir_tracking = os.path.join(wsdir, object_tracking.TRACKING_DIRNAME)
        os.makedirs(dir_tracking, exist_ok=True)
        dir_geoloc = os.path.join(wsdir, geometric_locating.GEOLOCATION_DIRNAME)
        os.makedirs(dir_geoloc, exist_ok=True)

        if not object_tracking.generate_detection_configs(dir_detection, cfg['detection_minimum_size_pointwise'],
                                         cfg['detection_minimum_size_linewise'], cfg['detection_maximum_size_linewise']):
            return False, "生成配置文件失败"
        tracking_params = object_tracking.generate_tracking_configs(dir_tracking, cfg['tracking_missing_frames'],
                                         cfg['tracking_minimum_frames'], cfg['target_minimum_speed'], cfg['target_maximum_speed'])
        if tracking_params is None :
            return False, "生成配置文件失败"

        if not geometric_locating.init(cfg): #or share.get('input_frames') is None or len(share['input_frames'])==0:
            return False, "轨迹生成配置初始化失败"

        bg_batchsize = cfg['background_frame_number']

        files = list()
        newfile_callback.files = files
        newfile_callback.extensions = extensions

        common.scan_existing_files(dir_input, newfile_callback)
        files.sort()
        file_watcher = TimeoutFileWatcher(newfile_callback, None, cfg['filewatcher_timeout'])

        file_watcher.start(dir_input)

        progress = 0
        bg_batches = []
        bg_idx_st = 0
        bg_idx_ed = 0
        share['ongoing_trajectories'] = [list(),list()]
        share['finished_trajectories'] = [list(),list()]
        share['3d_trajectories'] = []
        trackers = share['ongoing_trajectories']
        finished_trajectories = share['finished_trajectories']
        trajectory3d = share['3d_trajectories']
        num_trajectories = [0,0]
        while file_watcher.is_alive() or progress < len(files):
            filecount = min(len(files), progress + batchsize)
            newfile_count = filecount-progress
            if bg_idx_ed+bg_batchsize < filecount:
                bg_batch = (filecount-bg_idx_ed) // bg_batchsize
                bg_batches_new = []
                for i in range(bg_batch):
                    st = bg_idx_ed+i*bg_batchsize
                    bg_batches_new.append(
                        (os.path.join(dir_detection, f'bg_{st}_{st+bg_batchsize-1}.tif'), dir_input, files[st:st+bg_batchsize] ))
                print('1a) extracting background ...')
                parallel.launch_calls(background_extraction, bg_batches_new, nb_workers, cfg, timeout=cfg['timeout'])
                bg_batches += bg_batches_new
                bg_idx_ed += bg_batchsize*bg_batch
                i = min((progress-bg_idx_st) // bg_batchsize, len(bg_batches)-1)
                bg_idx_st += i*bg_batchsize
                bg_batches = bg_batches[i:]

            if len(bg_batches) > 0 and newfile_count > 0:
                anomaly_batches = []
                for i in range(newfile_count):
                    file = files[progress+i]
                    idx = min((progress+i-bg_idx_st)//bg_batchsize, len(bg_batches)-1)
                    anomaly_batches.append((os.path.join(dir_detection, os.path.basename(file)), os.path.join(
                        dir_input, file), bg_batches[idx][0]))
                print('1b) detecting anomaly objects ...')
                parallel.launch_calls(anomaly_detection, anomaly_batches, nb_workers, cfg, timeout=cfg['timeout'])

                print('1c) tracking anomaly objects ...')
                for i, detection_postfix in enumerate(cfg['detection_output_postfix']):
                    seedfiles = [x[0]+detection_postfix for x in anomaly_batches]
                    finished, trackers[i] = anomaly_tracking(progress, trackers[i], seedfiles, tracking_params, object_tracking.load_detection)
                    for tracker in finished:
                        name = 'M{}{}'.format(len(finished_trajectories[i])+1, cfg['tracking_output_postfix'][i])
                        finished_trajectories[i].append(object_tracking.save_tracking(dir_tracking, name, tracker, files))

                if not file_watcher.is_alive() and progress+newfile_count >= len(files):
                    for i, trajectories in enumerate(trackers):
                        if len(trajectories) == 0:
                            continue
                        for tracker in trajectories:
                            tracker.remove_all_missings()
                            if tracker.is_valid(tracking_params['min_frame_number'], tracking_params['min_speed'], tracking_params['max_acceleration'], tracking_params.get('curvature','spline')):
                                name = 'M{}{}'.format(len(finished_trajectories[i])+1, cfg['tracking_output_postfix'][i])
                                finished_trajectories[i].append(object_tracking.save_tracking(dir_tracking, name, tracker, files))
                        trackers[i].clear()

                new_trajectories = list()
                for st, x in zip(num_trajectories, finished_trajectories):
                    new_trajectories += x[st:]

                print('\t#NEW Targets: {}/{}'.format(len(new_trajectories), num_trajectories[0]+num_trajectories[1]+len(new_trajectories)))
                common.print_elapsed_time()

                if len(new_trajectories) > 0:

                    print('1d) refining trajectories ...')
                    batches = [(dir_tracking, info[-1], dir_input) for info in new_trajectories]
                    parallel.launch_calls(trajectory_file_refinement, batches, nb_workers, cfg, timeout=cfg['timeout'])

                    print('1f) geolocating trajectories ...')
                    st = num_trajectories[0]+num_trajectories[1]+1
                    batches = [( os.path.join(dir_output, f'M{st+i}.json' ), os.path.join( dir_tracking, trajectory[-1]), trajectory[-2] < 3 ) for i,trajectory in enumerate(new_trajectories) ]
                    parallel.launch_calls(trajectory_locating, batches, nb_workers, dir_input, None, cfg, timeout=cfg['timeout'])

                    new_trajectory3d = [(batch[0],batch[-1]) for batch in batches if os.path.exists(batch[0])]

                    print('1e) predicting trajectories ...')
                    batches = [( os.path.join(dir_output, f'M{st+i}_predicted.json' ), trajectory[0] ) for i,trajectory in enumerate(new_trajectory3d) ]
                    parallel.launch_calls(trajectory_predicting, batches, nb_workers, cfg, timeout=cfg['timeout'])
                    
                    new_trajectory3d = [ (batch[0],) if os.path.exists(batch[0]) else newt for batch, newt in zip(batches,new_trajectory3d) ]
                    trajectory3d += new_trajectory3d

                    num_trajectories = [len(x) for x in finished_trajectories]

                    for batch in new_trajectory3d:
                        logger.report(batch[0])
                
                    common.print_elapsed_time()

                progress += newfile_count
                logger.set_total(len(files)) 
                logger.progress_update(progress)

            if event.is_terminated():
                return True, "用户终止"

            while event.is_paused():
                time.sleep(0.5)
                if event.is_terminated():
                    return True, "用户终止"

        file_watcher.join()

    common.print_elapsed_time(True)
    return True, {'targets':trajectory3d}

def detection_processing(cfg, event, logger, share):

    common.print_elapsed_time.t0 = datetime.now()
    logger.set_step(0) 

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes'] is not None:
        nb_workers = cfg['max_processes']
    
    batchsize = cfg['background_frame_number']

    directory = cfg['input_dir']
    extensions = cfg['input_file_extensions']
    wsdir = cfg['temporary_dir']

    dir_detection = os.path.join(wsdir, object_tracking.DETECTION_DIRNAME)
    os.makedirs(dir_detection, exist_ok=True)

    if not object_tracking.generate_detection_configs(dir_detection, cfg['detection_minimum_size_pointwise'],
                                     cfg['detection_minimum_size_linewise'], cfg['detection_maximum_size_linewise']):
        return False, "生成配置文件失败"

    print(f'1) searching frames from {directory} ...')

    def newfile_callback(directory, filename):
        if any(filename.lower().endswith(ext) for ext in newfile_callback.extensions):
            newfile_callback.files.append(filename)
        return True

    share['input_frames'] = list()
    files = share['input_frames']
    newfile_callback.files = files
    newfile_callback.extensions = extensions

    common.scan_existing_files(directory, newfile_callback)
    files.sort()

    if len(files) < batchsize:
        print(f"[ERROR]: Not enough frames found in {directory} for background extraction")
        return False, f"影像数量不足，无法进行目标检测"

    frame_st = 0
    frame_ed = len(files)
    if len(files[0]) == len(cfg['start_frame']):
        frame_st = files.index(cfg['start_frame'])
    if len(files[0]) == len(cfg['end_frame']):
        frame_ed = files.index(cfg['end_frame'])+1

    if frame_ed-frame_st < batchsize:
        print(f"[ERROR]: Not enough frames found in {directory} for background extraction")
        return False, f"影像数量不足，无法进行目标检测"

    files = files[frame_st:frame_ed]

    print(f'2) detecting anomaly objects from #{len(files)} frames ...')
    logger.set_step(1) 
    logger.set_total(len(files)) 
    bg_batches = []
    output_prefix = []

    # def batch_process(i, sz):
    #     if sz == batchsize:
    #         ed = i+sz
    #         bg_batches = []
    #         bg_batches.append(
    #             (os.path.join(dir_detection, f'bg_{i}_{ed-1}.tif'), directory, files[i:ed] ))
    #         print('2a) extracting background ...')
    #         parallel.launch_calls(background_extraction, bg_batches, nb_workers, cfg, timeout=cfg['timeout'])

    #     anomaly_batches = [(os.path.join(dir_detection, os.path.basename(file)), os.path.join(
    #             directory, file), bg_batches[0][0]) for file in files[frame_st:frame_ed]]
    #     print('2b) detecting anomaly objects ...')
    #     parallel.launch_calls(anomaly_detection, anomaly_batches, nb_workers, cfg, timeout=cfg['timeout'])

    #     output_prefix += [x[0] for x in anomaly_batches]

    # success, msg = while_loop_with_events(batch_process, len(files), batchsize, event, logger)

    # if not success:
    #     return success, msg

    progress = 0
    while progress < len(files):
        frame_st = progress
        if frame_st+batchsize <= len(files):
            frame_ed = frame_st+batchsize
            bg_batches = []
            bg_batches.append(
                (os.path.join(dir_detection, f'bg_{frame_st}_{frame_ed-1}.tif'), directory, files[frame_st:frame_ed] ))
            print('2a) extracting background ...')
            parallel.launch_calls(background_extraction, bg_batches, nb_workers, cfg, timeout=cfg['timeout'])
        else:
            frame_ed = len(files)

        anomaly_batches = [(os.path.join(dir_detection, os.path.basename(file)), os.path.join(
                directory, file), bg_batches[0][0]) for file in files[frame_st:frame_ed]]
        print('2b) detecting anomaly objects ...')
        parallel.launch_calls(anomaly_detection, anomaly_batches, nb_workers, cfg, timeout=cfg['timeout'])

        output_prefix += [x[0] for x in anomaly_batches]

        progress = frame_ed
        logger.progress_update(progress)

        if event.is_terminated():
            return True, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return True, "用户终止"

    print('3) collecting detection results ...')
    logger.set_step(2) 
    logger.set_total(2) 

    summary_files = [os.path.join(cfg['output_dir'], 'detection_point_list.txt'), os.path.join(cfg['output_dir'], 'detection_line_list.txt')]
    for summary,postfix in zip( summary_files, cfg['detection_output_postfix']):
        with open(summary, 'w') as f:
            for prefix in output_prefix:
                f.write('{}\n'.format(prefix+postfix))

    common.print_elapsed_time(True)
    return True, {'points':summary_files[0], 'lines':summary_files[1]}

def tracking_processing(cfg, event, logger, share):

    common.print_elapsed_time.t0 = datetime.now()
    logger.set_step(0) 

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes'] is not None:
        nb_workers = cfg['max_processes']
    batchsize = cfg['background_frame_number']
    if cfg['batchsize'] is not None:
        batchsize = cfg['batchsize']

    wsdir = cfg['temporary_dir']
    dir_detection = os.path.join(wsdir, object_tracking.DETECTION_DIRNAME)
    dir_tracking = os.path.join(wsdir, object_tracking.TRACKING_DIRNAME)
    dir_input = cfg['input_dir']

    os.makedirs(dir_tracking, exist_ok=True)

    tracking_params = object_tracking.generate_tracking_configs(dir_tracking, cfg['tracking_missing_frames'],
                                     cfg['tracking_minimum_frames'], cfg['target_minimum_speed'], cfg['target_maximum_speed'])
    if tracking_params is None :
        return False, "生成配置文件失败"

    summary_files = [os.path.join(cfg['output_dir'], 'detection_point_list.txt'), os.path.join(cfg['output_dir'], 'detection_line_list.txt')]

    for i,file in enumerate(summary_files):
        try:
            with open(file, 'r') as f:
                summary_files[i] = [ l.strip() for l in f]
                if len(summary_files[i]) < cfg['tracking_minimum_frames']:
                    raise ValueError('frame number is less than tracking tolerance')
        except Exception as e :
            print(f"[ERROR]: detection file list {file}: {e}!")
            return False, "目标检测文件列表载入失败"

    assert(len(summary_files[0]) == len(summary_files[1]))
    filecount = len(summary_files[0])
    print(f'1) checking detection results: #{filecount} point files and line files ...')
    for i, files in enumerate(summary_files):
        if os.path.exists(files[0]):
            continue
        file = os.path.join(dir_detection, files[0])
        if not os.path.exists(file):
            print(f"[ERROR]: detection files NOT exist!")
            return False, "目标检测文件载入失败"
        summary_files[i] = [os.path.join(dir_detection, file) for file in files]

    print('2) tracking anomaly objects ... ')
    logger.set_step(1)
    logger.set_total(filecount)
    share['ongoing_trajectories'] = [list(),list()]
    share['finished_trajectories'] = [list(),list()]
    trackers = share['ongoing_trajectories']
    finished_trajectories = share['finished_trajectories']

    frames = [os.path.basename(x)[:-len(cfg['detection_output_postfix'][0])] for x in summary_files[0]]

    progress = 0
    while progress < filecount:
        num = min(batchsize, filecount-progress)

        for i, postfix in enumerate(cfg['tracking_output_postfix']):
            seedfiles = summary_files[i][progress:progress+num]
            finished, trackers[i] = anomaly_tracking(progress, trackers[i], seedfiles, tracking_params, object_tracking.load_detection)
            for tracker in finished:
                name = 'M{}{}'.format(len(finished_trajectories[i])+1, postfix)
                finished_trajectories[i].append(object_tracking.save_tracking(dir_tracking, name, tracker, frames))
        
        progress += num
        logger.progress_update(progress)

        if event.is_terminated():
            return True, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return True, "用户终止"

    for i, trajectories in enumerate(trackers):
        if len(trajectories) == 0:
            continue
        for tracker in trajectories:
            tracker.remove_all_missings()
            if tracker.is_valid(tracking_params['min_frame_number'], tracking_params['min_speed'], tracking_params['max_acceleration'], tracking_params.get('curvature','spline')):
                name = 'M{}{}'.format(len(finished_trajectories[i])+1, cfg['tracking_output_postfix'][i])
                finished_trajectories[i].append(object_tracking.save_tracking(dir_tracking, name, tracker, frames))
        trackers[i].clear()
    
    print('3) collecting tracking results ...')
    trajectories = [item for sublist in finished_trajectories for item in sublist]
    progress = 0
    filecount = len(trajectories)
    logger.set_step(2)
    logger.set_total(filecount)

    while progress < filecount:
        num = min(batchsize, filecount-progress)

        batches = [(dir_tracking, info[-1], dir_input) for info in trajectories[progress:progress+num]]
        parallel.launch_calls(trajectory_file_refinement, batches, nb_workers, cfg, timeout=cfg['timeout'])

        progress += num

        if event.is_terminated():
            return True, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated(): 
                return True, "用户终止"

    summary_files = [os.path.join(cfg['output_dir'], 'tracking_point_list.csv'), os.path.join(cfg['output_dir'], 'tracking_line_list.csv')]
    for summary, trajectories in zip( summary_files, finished_trajectories):
        with open(summary, 'w') as f:
            for trajectory in trajectories:
                f.write('{}\n'.format(','.join(f"{x:.2f}" if isinstance(x,float) else str(x) for x in trajectory)))

    common.print_elapsed_time(True)
    return True, {'points':summary_files[0], 'lines':summary_files[1]}

def geolocating_processing(cfg, event, logger, share):

    common.print_elapsed_time.t0 = datetime.now()
    logger.set_step(0) 

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes'] is not None:
        nb_workers = cfg['max_processes']
    batchsize = cfg['background_frame_number']
    if cfg['batchsize'] is not None:
        batchsize = cfg['batchsize']

    wsdir = cfg['temporary_dir']
    dir_tracking = os.path.join(wsdir, object_tracking.TRACKING_DIRNAME)
    dir_geoloc = os.path.join(wsdir, geometric_locating.GEOLOCATION_DIRNAME)
    dir_input = cfg['input_dir']
    dir_output = cfg['output_dir']

    os.makedirs(dir_geoloc, exist_ok=True)

    trajectories = [os.path.join(cfg['output_dir'], 'tracking_point_list.csv'), os.path.join(cfg['output_dir'], 'tracking_line_list.csv')]
    for i,file in enumerate(trajectories):
        try:
            with open(file, 'r') as f:
                trajectories[i] = [l.strip().split(',')[-2:] for l in f]
        except Exception as e:
            print(f"[ERROR]: trajectory file list {file}: {e}!")
            return False, "目标跟踪没有结果"

    filecount = len(trajectories[0])+len(trajectories[1])
    print('1) checking #{} tracking results: #{} point files and #{} line files ...'.format(filecount, len(trajectories[0]), len(trajectories[1])))

    if filecount == 0:
        print(f"[ERROR]: No tracking results found!")
        return False, "目标跟踪没有结果"

    if not geometric_locating.init(cfg): #or share.get('input_frames') is None or len(share['input_frames'])==0:
        return False, "轨迹生成配置初始化失败"

    for i, files in enumerate(trajectories):
        if len(files) == 0 or os.path.exists(files[0][-1]):
            continue
        file = os.path.join(dir_tracking, files[0][-1])
        if not os.path.exists(file):
            print(f"[ERROR]: tracking files NOT exist!")
            return False, "目标跟踪文件载入失败"
        trajectories[i] = [ [ float(info[0]), os.path.join(dir_tracking, info[-1])] for info in files]

    frametime = None

    print('2) geolocating trajectories ... ')
    logger.set_step(1)
    logger.set_total(filecount)
    share['3d_trajectories'] = []
    trajectory3d = share['3d_trajectories']

    trajectories = trajectories[0]+trajectories[1]
    progress = 0
    filecount = len(trajectories)
    while progress < filecount:
        num = min(batchsize, filecount-progress)

        batches = [( os.path.join(dir_output, f'M{progress+i+1}.json' ), trajectory[1], trajectory[0] < 3 ) for i,trajectory in enumerate(trajectories[progress:progress+num]) ]
        parallel.launch_calls(trajectory_locating, batches, nb_workers, dir_input, frametime, cfg, timeout=cfg['timeout'])

        trajectory3d += [(batch[0],) for batch in batches if os.path.exists(batch[0])]

        progress += num
        logger.progress_update(progress)

        if event.is_terminated():
            return True, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return True, "用户终止"

    print('3) collecting geolocating results ...')
    with open(os.path.join(dir_output, 'targets.csv'), 'w') as f:
        for trajectory in trajectory3d:
            f.write(','.join([str(x) for x in trajectory])+'\n')

    logger.report([x[0] for x in trajectory3d], cfg['orderjson'])

    common.print_elapsed_time(True)
    return True, {'3d_trajectories':trajectory3d }

def trajectory_predicting_processing(cfg, event, logger, share):
    def newfile_callback(directory, filename):
        if any(filename.lower().endswith(ext) for ext in ['.json']):
            newfile_callback.files.append(os.path.join(directory, filename))
        return True

    common.print_elapsed_time.t0 = datetime.now()
    logger.set_step(0) 

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes'] is not None:
        nb_workers = cfg['max_processes']
    batchsize = cfg['background_frame_number']
    if cfg['batchsize'] is not None:
        batchsize = cfg['batchsize']

    dir_input = cfg['input_dir']
    dir_output = cfg['output_dir']

    share['predicted_3d_trajectories'] = list()
    if not share.get('3d_trajectories'):
        share['3d_trajectories'] = list()
    trajectory3d = share['3d_trajectories']
    predicted_trajectory3d = share['predicted_3d_trajectories']

    if len(trajectory3d) == 0:
        files = os.path.join(dir_output, 'targets.csv')
        if not os.path.exists(files):
            newfile_callback.files = list()
            common.scan_existing_files(dir_input, newfile_callback)
            newfile_callback.files.sort()
            if len(newfile_callback.files) == 0:
                print(f"[ERROR]: No 3D trajectories found!")
                return False, "三维轨迹没有结果"
            trajectory3d += [(file, ) for file in newfile_callback.files]
        else:
            with open(files, 'r') as f:
                trajectory3d += [(file.strip(),) for file in f]


    if not geometric_locating.init(cfg): #or share.get('input_frames') is None or len(share['input_frames'])==0:
        return False, "轨迹生成配置初始化失败"

    print(f'1) predicting #{len(trajectory3d)} trajectories ... ')

    progress = 0
    filecount = len(trajectory3d)
    while progress < filecount:
        num = min(batchsize, filecount-progress)

        batches = [( os.path.join(dir_output, f'M{progress+i+1}_predicted.json' ), trajectory[0] ) for i,trajectory in enumerate(trajectory3d[progress:progress+num]) ]
        parallel.launch_calls(trajectory_predicting, batches, nb_workers, cfg, timeout=cfg['timeout'])

        predicted_trajectory3d += [(batch[0],) for batch in batches if os.path.exists(batch[0])]
        
        progress += num
        logger.progress_update(progress)

        if event.is_terminated():
            return True, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return True, "用户终止"

    print('2) collecting predicting results ...')

    with open(os.path.join(dir_output, 'predicted_targets.csv'), 'w') as f:
        for trajectory in predicted_trajectory3d:
            f.write(','.join([str(x) for x in trajectory])+'\n')

    logger.report([x[0] for x in predicted_trajectory3d], cfg['orderjson'])

    common.print_elapsed_time(True)
    return True, {'predicted_3d_trajectories':predicted_trajectory3d}

def draw_seeds( framepath, seedfile):
    image = common.rasterio_read_as_rgb24(framepath)
    seeds = object_tracking.load_detection(seedfile)[0]

    if len(seeds) > 0:
        common.draw(image, seeds, (0,255,0), kwargs={'showCoordinates':True})
        cv2.putText(image, f'{os.path.basename(framepath)}', (0,image.shape[0]-10), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,255,0))

    return image
    
def draw_trajectory(frames, trajectory, label = None, color = (0, 255, 0)):
    bg = frames
    if isinstance(trajectory, str):
        if label is None:
            label = os.path.splitext(os.path.basename(trajectory))[0]
        st, _, trajectory = object_tracking.load_tracking(trajectory)
        trajectory = [(x[1], x[2]) for x in trajectory]
        bg = frames[st:]

    n = min(len(bg), len(trajectory))
    for i in range(2, n):
        common.draw( bg[i], trajectory[:i-1], color, kwargs={'linked':True, 'radius':1})
        common.draw( bg[i], [trajectory[i]], color)

    return frames

def trajectory_video( output:str, trajectory: str, framedir : str, clip = True, margin = 256):
    label = os.path.splitext(os.path.basename(trajectory))[0]
    st, cnt, trajectory = object_tracking.load_tracking(trajectory)

    if cnt == 0:
        return

    win = None
    xy = np.array([[x[1], x[2]] for x in trajectory])

    if clip:
        width, height,_ = common.get_image_shape(os.path.join(framedir, trajectory[0][0]))
        win = list(common.bounding_box2D(xy))
        win[0] = max(0, win[0]-margin)
        win[1] = max(0, win[1]-margin)
        win[2] = min(width-win[0], win[2]+2*margin)
        win[3] = min(height-win[1], win[3]+2*margin)
        xy -= np.array([win[0], win[1]])

    frames = []
    for info in trajectory:
        frames.append(common.rasterio_read_as_rgb24(os.path.join(framedir, info[0]), win))
        cv2.putText(frames[-1], f'{info[0]}', (0,frames[-1].shape[0]-10), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,255,0))

    frames = draw_trajectory(frames, xy, label, common.generate_random_color(cnt))

    height, width = frames[0].shape[:2]
    video = cv2.VideoWriter(output, cv2.VideoWriter_fourcc(*'mp4v'), 1, (width, height))
    for frame in frames:
        video.write(frame)

    video.release()


if __name__ == "__main__":
    import sys
    with open(sys.argv[1], 'r') as f:
        xml_string = f.read()
        with open(sys.argv[2], 'w') as f:
            f.write(xml_string_to_json(xml_string))
