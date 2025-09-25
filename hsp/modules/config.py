import json
import os

DEFAULT_CFG = {
    'output_dir':'hsp_output',
    'temporary_dir': 'hsp_tmp',
    'clean_tmp': True,
    'max_processes': None,
    'timeout': 1800,
    'overwritten': True,
    'input_file_extensions': ['.tif', '.png'],
    'filewatcher_timeout': 30,
    'omp_num_threads': 1,
    'batchsize': 100,

    'background_frame_number': 10,
    'background_method': 'mean',
    'detection_output_postfix': ['_points.txt', '_lines.txt'],
    'detection_minimum_size_pointwise': 1,
    'detection_minimum_size_linewise': 4,
    'detection_maximum_size_linewise': 10,

    'tracking_missing_frames': 30,
    'tracking_minimum_frames': 30,
    'target_minimum_speed': 0.1,
    'target_maximum_speed': 2,
    'tracking_output_postfix': ['_points.txt', '_lines.txt'],

    'recognition_airplane_angular_velocity': 3,
    'ground_altitude': 0,

    'sofa': 'EOP.txt',
    'orderjson': "order.json",

    'debug': False
}

def load_config(config_file):
    """
    Load a configuration file in JSON format, and merge it with the default
    configuration.

    Args:
        config_file: path to the configuration file

    Returns:
        dictionary containing the configuration
    """
    cfg = DEFAULT_CFG.copy()
    if config_file:
        with open(config_file) as f:
            user_cfg = json.load(f)
            cfg.update(user_cfg)
    return cfg

def save_config(config_file, config):
    """
    Save a configuration file in JSON format.

    Args:
        config_file: path to the configuration file
        config: dictionary containing the configuration
    """
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=4)
        return True
    return False

def orderjson_to_config(orderxml):
    try:
        with open(orderxml, encoding='utf-8') as f:
            data = json.load(f)
            cfg = DEFAULT_CFG.copy()
            params = data['InterfaceFile']['FileBody']
            cfg['input_dir'] = params['InputFile']
            cfg['output_dir'] = params['WorkDir']
            cfg['temporary_dir'] = params['LocalWorkDir']
            cfg['start_frame'] = params.get('StartFrame')
            cfg['end_frame'] = params.get('EndFrame')
            if cfg['start_frame'] is None:
                cfg['start_frame'] = ''
            if cfg['end_frame'] is None:
                cfg['end_frame'] = ''

            params = data['InterfaceFile']['Parameters']
            try:
                cfg['background_frame_number'] = int(params['BackgroundFrames'])
                cfg['tracking_missing_frames'] = int(params['MaxMissingFrames'])
                cfg['tracking_minimum_frames'] = int(params['MinDetectionFrames'])
                cfg['target_maximum_size'] = float(params['TargetMaximumSize'])/400
                cfg['target_minimum_speed'] = float(params['MinSpeed'])/400
                cfg['target_maximum_speed'] = float(params['MaxSpeed'])/400
            except Exception as e:
                print(f"[WARNING]: {e}")

            if not os.path.exists(cfg['input_dir']):
                return None
            return cfg
    except Exception as e:
        print(f"[ERROR]: {e}")
    return None