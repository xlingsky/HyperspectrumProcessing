import os
import sys
import json
import xml.etree.ElementTree as ET
import time
from datetime import datetime
import multiprocessing
import numpy as np

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
            return False, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return False, "用户终止"

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
        output_info = {"trajectory_id": os.path.splitext(os.path.basename(output))[ 0],
                       "datatime": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), 
                       "SatelliteList": { "Satellite":[{"ID": "1"} ]}, "PointList": { "Point": [] } }
        start, cnt, points = object_tracking.load_tracking(trajectory)
        
        if cnt <= 0:
            return 

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
                "ImageCoordinates": f"{point[1]:.1f},{point[2]:.1f}",
                "DigitalNumber": f"{point[3]:.1f}",
                "Energy": f"{point[3]*0.01:.1f}"
                }

            if locating.load_rpc(frame):

                z = z0+acceleration*(i)**2/2

                graphic, proj, centric = locating.transform(point[1], point[2], z)

                info['Location'] = f"{graphic[0]:.6f},{graphic[1]:.6f},{z:.1f}"
                info['Projection'] = [proj[0], proj[1], z]
                info['CGCS2000'] = f"{centric[0]:.1f},{centric[1]:.1f},{centric[2]:.1f}"
            
                if sofa_transformer is not None:
                    cgcs2j = np.linalg.inv(sofa_transformer.j2000_to_cgcs2000_matrix(info['Time'])) 

                    j2000 = cgcs2j @ np.array(centric).reshape(3,1)
                    info['J2000'] = f"{j2000[0,0]:.1f},{j2000[1,0]:.1f},{j2000[2,0]:.1f}"

            points_info.append(info)

        points_info[0]['WarningStatus'] = '01H'
        points_info[-1]['WarningStatus'] = '03H'
        for i in range(1,len(points_info)-1):
            points_info[i]['WarningStatus'] = '02H'

        dist = 0
        for i in range(0,len(points_info)-1):
            diff = np.array(points_info[i+1]['Projection'])-np.array(points_info[i]['Projection'])
            dist += np.linalg.norm(diff)
            points_info[i]['Velocity'] = f"{diff[0]:.1f},{diff[1]:.1f},{diff[2]:.1f}"

        output_info['total_distance'] = f"{dist:.1f}"

        with open(output, 'w') as fout: 
            json.dump({"Trajectory":output_info}, fout, indent=2, ensure_ascii=False)

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
                        name = 'M{}{}'.format(len(finished_trajectories[i]+1), cfg['tracking_output_postfix'][i])
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

                print('\t#NEW Targets: {}'.format(len(new_trajectories)))
                common.print_elapsed_time()

                if len(new_trajectories) > 0:

                    print('1d) refining trajectories ...')
                    batches = [(dir_tracking, info[-1], dir_input) for info in new_trajectories]
                    parallel.launch_calls(trajectory_file_refinement, batches, nb_workers, cfg, timeout=cfg['timeout'])

                    print('1f) geolocating trajectories ...')
                    st = num_trajectories[0]+num_trajectories[1]+1
                    batches = [( os.path.join(dir_output, f'M{st+i}.json' ), trajectory[1], trajectory[0] < 3 ) for i,trajectory in enumerate(new_trajectories) ]
                    parallel.launch_calls(trajectory_locating, batches, nb_workers, dir_input, None, cfg, timeout=cfg['timeout'])

                    trajectory3d += [(batch[0],) for batch in batches]

                    num_trajectories = [len(x) for x in finished_trajectories]
                
                    common.print_elapsed_time()

                progress += newfile_count
                logger.set_total(len(files)) 
                logger.progress_update(progress)

            if event.is_terminated():
                return False, "用户终止"

            while event.is_paused():
                time.sleep(0.5)
                if event.is_terminated():
                    return False, "用户终止"

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
            return False, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return False, "用户终止"

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
            return False, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return False, "用户终止"

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
            return False, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return False, "用户终止"

    summary_files = [os.path.join(cfg['output_dir'], 'tracking_point_list.csv'), os.path.join(cfg['output_dir'], 'tracking_line_list.csv')]
    for summary, trajectories in zip( summary_files, finished_trajectories):
        with open(summary, 'w') as f:
            for trajectory in trajectories:
                f.write('{}\n'.format(','.join(f"{x:.2f}" if isinstance(x,float) else str(x) for x in trajectory)))

    common.print_elapsed_time(True)
    return True, {'points':summary_files[0], 'lines':summary_files[1]}

def geolocating_process(cfg, event, logger, share):

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
            return False, "二维轨迹文件列表载入失败"

    filecount = len(trajectories[0])+len(trajectories[1])
    print('1) checking #{} tracking results: #{} point files and #{} line files ...'.format(filecount, len(trajectories[0]), len(trajectories[1])))

    if filecount == 0:
        print(f"[ERROR]: No tracking results found!")
        return False, "二维轨迹文件列表载入失败"

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

    progress = 0
    filecount = len(trajectories[0])
    while progress < filecount:
        num = min(batchsize, filecount-progress)

        batches = [( os.path.join(dir_output, f'M{progress+i+1}.json' ), trajectory[1], trajectory[0] < 3 ) for i,trajectory in enumerate(trajectories[0][progress:progress+num]) ]
        parallel.launch_calls(trajectory_locating, batches, nb_workers, dir_input, frametime, cfg, timeout=cfg['timeout'])

        trajectory3d += [(batch[0],) for batch in batches]

        progress += num
        logger.progress_update(progress)

        if event.is_terminated():
            return False, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return False, "用户终止"

    filecount = len(trajectories[1])
    while progress < filecount:
        num = min(batchsize, filecount-progress)

        batches = [( os.path.join(dir_output, f'M{progress+i+1}.json' ), trajectory[1], True ) for i,trajectory in enumerate(trajectories[1][progress:progress+num]) ]
        parallel.launch_calls(trajectory_locating, batches, nb_workers, dir_input, frametime, cfg, timeout=cfg['timeout'])

        trajectory3d += [(batch[0],) for batch in batches]

        progress += num
        logger.progress_update(progress)

        if event.is_terminated():
            return False, "用户终止"

        while event.is_paused():
            time.sleep(0.5)
            if event.is_terminated():
                return False, "用户终止"

    print('3) collecting geolocating results ...')
    with open(os.path.join(dir_output, 'targets.csv'), 'w') as f:
        for trajectory in trajectory3d:
            f.write(','.join([str(x) for x in trajectory])+'\n')

    common.print_elapsed_time(True)
    return True, {'3d_trajectories':trajectory3d }

if __name__ == "__main__":
    import sys
    with open(sys.argv[1], 'r') as f:
        xml_string = f.read()
        with open(sys.argv[2], 'w') as f:
            f.write(xml_string_to_json(xml_string))