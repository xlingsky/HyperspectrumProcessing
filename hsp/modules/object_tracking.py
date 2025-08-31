import numpy as np
import os
from scipy.ndimage import median_filter
from scipy.spatial import KDTree
from sklearn.decomposition import PCA
from typing import List
import cv2
import numpy as np
import json

from hsp.utils import common

DETECTION_DIRNAME = 'detection'
TRACKING_DIRNAME = 'tracking'

def extract_background(frames, method='median', pca_components=3):
    """
    Extract background from multiple frames using different methods
    
    Parameters:
    frames -- Input video/hyperspectral sequence ( T x H x W )
              H = height, W = width, T = number of frames
    method -- Background extraction method:
              'median' - temporal median (default)
              'mean' - temporal mean
              'pca' - principal component reconstruction
              'morph' - morphological filtering
    pca_components -- Number of components to keep for PCA method
    
    Returns:
    background -- Estimated background (H x W )
    """
    T, H, W = frames.shape
    dimT = 0
    
    if method == 'median':
        # Temporal median - robust to transient anomalies
        background = np.median(frames, axis=dimT)
        
    elif method == 'mean':
        # Simple temporal mean
        background = np.mean(frames, axis=dimT)
        
    elif method == 'pca':
        # PCA-based background modeling
        pca = PCA(n_components=pca_components)
        flattened = frames.reshape(T, -1)  # Flatten spatial dimensions
        
        # Fit PCA and reconstruct background
        pca.fit(flattened)
        reconstructed = pca.inverse_transform(pca.transform(flattened))
        background = np.median(reconstructed, axis=0).reshape(H, W)
        
    elif method == 'morph':
        # Morphological approach (median filtering in space and time)
        # First apply temporal median
        temp_median = np.median(frames, axis=dimT)

        # Then spatial median filtering
        background = median_filter(temp_median, size=3)
    
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'median', 'mean', 'pca', or 'morph'")
    
    return background

def generate_detection_configs(configdir, point_size, line_min, line_max):
    templatedir = os.path.join(common.parent_dir, 'template')
    try:
        rx = json.load(open(os.path.join(templatedir, 'tsk_rx.json')))
        detection = json.load(open(os.path.join(templatedir, 'tsk_detection.json')))
        all = json.load(open(os.path.join(templatedir, 'tsk_rx_detection.json')))

        detection['HSP']['task'][0]['length_min'] = line_min
        detection['HSP']['task'][0]['length_max'] = line_max

        all['HSP']['task'][2]['length_min'] = line_min
        all['HSP']['task'][2]['length_max'] = line_max

        with open(os.path.join(configdir, 'tsk_rx.json'), 'w') as f:
            f.write(json.dumps( rx, indent=2, ensure_ascii=False))
        with open(os.path.join(configdir, 'tsk_detection.json'), 'w') as f:
            f.write(json.dumps( detection, indent=2, ensure_ascii=False))
        with open(os.path.join(configdir, 'tsk_rx_detection.json'), 'w') as f:
            f.write(json.dumps( all, indent=2, ensure_ascii=False))

        return True
    except Exception as e:
        print(f"[ERROR]: {e}")
        return False


def generate_tracking_configs(configdir, tracking_missing_frames, tracking_minimum_frames, target_minimum_speed, target_maximum_speed):
    templatedir = os.path.join(common.parent_dir, 'template')
    try:
        tracking = json.load(open(os.path.join(templatedir, 'tsk_tracking.json')))
        tracking['min_frame_number'] = tracking_minimum_frames
        tracking['min_speed'] = target_minimum_speed
        tracking['kalman']['max_missing_frames'] = tracking_missing_frames
        tracking['kalman']['search_radius'] = target_maximum_speed
        with open(os.path.join(configdir, 'tsk_tracking.json'), 'w') as f:
            f.write(json.dumps( tracking, indent=2, ensure_ascii=False))
        return tracking
    except Exception as e:
        print(f"[ERROR]: {e}")
        return None

def detect_anomaly( image, configdir, bg, output_prefix, output_feature = True):
    if output_feature:
        rx = output_prefix+'_rx.tif'
        common.run('{0} -task {1} -o {2} {3} {4} -v 0'.format('hsp',
               os.path.join(configdir, 'tsk_rx.json'), rx, image, bg))
        common.run('{0} -task {1} -o {2} {3} -v 0'.format('hsp',
               os.path.join(configdir, 'tsk_detection.json'), output_prefix, rx ))
    else:
        common.run('{0} -task {1} -o {2} {3} {4} -v {5}'.format('hsp',
               os.path.join(configdir, 'tsk_rx_detection.json'), output_prefix, image, bg, 0))

class KalmanTracker:
    def __init__(self, frame_start: int, ptid: int, pt: np.ndarray, params: dict):
        # State vector: [x, y, vx, vy] - position (x,y) and velocity (vx,vy)
        self.state_size = 4
        # Measurement vector: [x, y] - we can only observe position
        self.measurement_size = 2
        
        # Initialize Kalman filter
        self.kf = cv2.KalmanFilter(self.state_size, self.measurement_size, 0)
        
        # Transition Matrix (A)
        self.kf.transitionMatrix = np.array(params['transition_matrix'], dtype=np.float32)
        
        # Measurement Matrix (H)
        self.kf.measurementMatrix = np.array([[1, 0, 0, 0],
                                             [0, 1, 0, 0]], dtype=np.float32)
        
        # Process Noise Covariance Matrix (Q)
        self.kf.processNoiseCov = np.eye(self.state_size, dtype=np.float32) * params['process_noise']
        
        # Measurement Noise Covariance Matrix (R)
        self.kf.measurementNoiseCov = np.eye(self.measurement_size, dtype=np.float32) * params['measurement_noise']
        
        # Error Covariance Matrix (P)
        self.kf.errorCovPost = np.eye(self.state_size, dtype=np.float32) * params['error_cov_post']
        
        # Initial State Posterior
        self.kf.statePost = np.array([[pt[0]], [pt[1]], [0], [0]], dtype=np.float32)
        
        # Tracker-specific properties
        self._point_ids = [ptid]
        self._points = [pt]
        self._frame_start = frame_start
        self._missings = 0
        self._tolerance = params['max_missing_frames']
        self._search_radius = params['search_radius']
    
    @property
    def good(self) -> bool:
        return self._missings <= self._tolerance
    
    def search_radius(self) -> float:
        return self._search_radius if self._missings == 0 else self._search_radius * 2
    
    def update(self, id: int, pt: np.ndarray):
        self._point_ids.append(id)
        self._points.append(pt) 
        self.kf.correct(np.array([[pt[0]],[pt[1]]], dtype=np.float32))
        self._missings = 0
    
    def missing(self, pt: np.ndarray):
        self._point_ids.append(-1)
        self._points.append(pt)
        self.kf.correct(np.array([[pt[0]],[pt[1]]], dtype=np.float32))
        self._missings += 1
    
    def remove_all_missings(self):
        self._point_ids = self._point_ids[:len(self._point_ids) - self._missings]
        self._points = self._points[:len(self._points) - self._missings]
    
    def estimate(self) -> np.ndarray:
        prediction = self.kf.predict()
        return prediction.reshape(-1)[:2]
    
    def is_valid(self, min_frame_number: int, min_speed: float, max_acceleration: float) -> bool:
        if len(self._points) < min_frame_number:
            return False
        
        distance = 0.0
        speeds = []
        
        for i in range(1, len(self._points)):
            pt1 = self._points[i - 1]
            pt2 = self._points[i]
            dx = pt2[0] - pt1[0]
            dy = pt2[1] - pt1[1]
            v = dx if abs(dx) > abs(dy) else dy
            speeds.append(v)
            distance += abs(v)
            
            if (distance - i * min_speed) < -1:
                return False
        
        if len(self._points) > 5:
            data = np.array(speeds, dtype=np.float32).reshape(1, -1)
            min_val, max_val = np.min(data), np.max(data)
            
            if min_val < 0 and max_val > 0:
                hist_size = max(5, int(np.ceil(max_val - min_val)))
                hist_range = (min_val, max_val)
                hist = cv2.calcHist(data, [0], None, [hist_size], hist_range)
                
                id0 = int(-hist_size * min_val / (max_val - min_val))
                if hist[id0] > max_acceleration * len(self._points):
                    return False
                
                if hist[0] + hist[hist_size - 1] > max_acceleration * len(self._points):
                    return False
        
        return True

# def tracking(seeds: List[np.ndarray], trackers: List[KalmanTracker], frameid: int, params: dict, ):
def pointwise_tracking(seeds : list, trackers : list, frameid: int, params : dict, num_nn: int = 1):
    """Track objects using Kalman filter and KD-tree"""
    status = [0] * len(seeds)
    
    if len(trackers) > 0:
        predicted_points = [tracker.estimate() for tracker in trackers]
    
        if len(seeds) == 0:
            for tracker, pt in zip(trackers,predicted_points):
                tracker.missing(pt)
            return trackers
    
        # Build KD-tree
        data = np.array(seeds)
        kdtree = KDTree(data[:,:2])
    
        for tracker, pt in zip(trackers,predicted_points):
            # Find nearest neighbors
            dist, idx = kdtree.query(pt, k=num_nn)

            search_radius = tracker.search_radius()

            if status[idx] == 0 and dist < search_radius:
                seed_pt = data[idx]
                tracker.update(idx, seed_pt)
                status[idx] = 1
            else:
                tracker.missing(pt)
    
    # Add new trackers for unused seeds
    for i, used in enumerate(status):
        if used == 0:
            trackers.append(KalmanTracker(frameid, i, seeds[i], params))

    return trackers

def load_detection(file: str):
    all = []
    try:
        with open(file, 'r') as f:
            i = 0
            lines = f.readlines()
            layers = int(lines[i].strip())
            i += 1
            for j in range(layers):
                cnt = int(lines[i].strip())
                i += 1
                points = [[float(x) for x in line.strip().split()] for line in lines[i:i+cnt]]
                i += cnt
                all.append(points)
                    
    except Exception as e:
        print(f"[ERROR]: {e}")
    return all

def save_tracking(directory: str, file: str, tracker: KalmanTracker, frames:list):
    try:
        with open(os.path.join(directory, file), 'w') as f:
            f.write('{} {}\n'.format(tracker._frame_start, len(tracker._points)))
            for i, pt in enumerate(tracker._points):
                f.write('{}\t{}\t{}\n'.format( frames[tracker._frame_start+i], '\t'.join(str(x) for x in pt), tracker._point_ids[i]))
    except Exception as e:
        print(f"[ERROR]: {e}")

    return [tracker._frame_start, len(tracker._points), file]