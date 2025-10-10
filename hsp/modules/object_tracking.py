import numpy as np
import os
from scipy.ndimage import median_filter
from scipy.spatial import KDTree
from sklearn.decomposition import PCA
from typing import List
import cv2
import numpy as np
import json
import rasterio

from hsp.utils import common, curvature

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
    try:
        if output_feature:
            rx = output_prefix+'_rx.tif'
            common.run('{0} -task {1} -o {2} {3} {4} -v 0'.format('hsp',
                   os.path.join(configdir, 'tsk_rx.json'), rx, image, bg))
            common.run('{0} -task {1} -o {2} {3} -v 0'.format('hsp',
                   os.path.join(configdir, 'tsk_detection.json'), output_prefix, rx ))
        else:
            common.run('{0} -task {1} -o {2} {3} {4} -v {5}'.format('hsp',
                   os.path.join(configdir, 'tsk_rx_detection.json'), output_prefix, image, bg, 0))
    except Exception as e:
        print(f"[ERROR]: {e}")

class GlobalOffsetList:
    def __init__(self, frame_id: int = 0, offset: np.ndarray = np.zeros(2), point_set: List = []):
        self.start_frame_id = frame_id
        self.start_offset = offset
        self.last_point_set = point_set
        self.last_offset_id = 0
        self.offsets = [np.zeros(2)]

    def offset(self, frame_id: int) -> np.ndarray:
        assert frame_id >= self.start_frame_id and frame_id < self.start_frame_id + len(self.offsets)
        return self.offsets[frame_id-self.start_frame_id]

    def global_offset(self, frame_id: int) -> np.ndarray:
        return self.start_offset + self.offset(frame_id)
    
    def append(self, point_set: List):
        if len(point_set) < 3 :
            self.offsets.append(self.offsets[self.last_offset_id])
            return
        elif len(self.last_point_set) < 3 :
            offset = np.zeros(2)
        else:
            points = np.concatenate( (np.array(self.last_point_set)[:,:2], np.array(point_set)[:,:2]), axis=0 )
            xmin, ymin = points.min(axis=0)
            xmax, ymax = points.max(axis=0)
            width = np.ceil(xmax - xmin+1).astype(np.int32)
            height = np.ceil(ymax - ymin+1).astype(np.int32)
            ref = np.zeros((height, width), dtype=np.uint8)
            for p in self.last_point_set:
                x = int(np.round(p[0]-xmin))
                y = int(np.round(p[1]-ymin))
                ref[y,x] = 1
            src = np.zeros((height, width), dtype=np.uint8)
            for p in point_set:
                x = int(np.round(p[0]-xmin))
                y = int(np.round(p[1]-ymin))
                src[y,x] = 1
            pad_width = max(len(self.offsets)-self.last_offset_id, 5)
            padded = np.pad(src, pad_width=pad_width, mode='constant', constant_values=0)
            res = cv2.matchTemplate(padded, ref, cv2.TM_SQDIFF)
            min_val, _, min_loc, _ = cv2.minMaxLoc(res)
            if min_val < 0.2*max(len(self.last_point_set), len(point_set)):
                offset = np.array([pad_width, pad_width]) - np.array(min_loc)
            else:
                offset = np.zeros(2)
        self.offsets.append(offset+self.offsets[self.last_offset_id])
        self.last_point_set = point_set
        self.last_offset_id = len(self.offsets)-1

    def save(self, filename: str):
        with open(filename, 'w') as f:
            for offset in self.offsets:
                f.write('{0} {1}\n'.format(offset[0], offset[1]))

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
        self._response_dieout_ratio = params['response_dieout_ratio']
        self._response = 0
    
    @property
    def good(self) -> bool:
        return self._missings <= self._tolerance
    
    def search_radius(self) -> float:
        return self._search_radius if self._missings == 0 else self._search_radius * 2
    
    def update(self, id: int, pt: np.ndarray):
        if pt[2] < self._response*self._response_dieout_ratio:
            return False
        self._point_ids.append(id)
        self._points.append(pt) 
        self.kf.correct(np.array([[pt[0]],[pt[1]]], dtype=np.float32))
        self._missings = 0
        self._response = pt[2]
        return True
    
    def missing(self, pt: np.ndarray):
        self._point_ids.append(-1)
        self._points.append(pt)
        self.kf.correct(np.array([[pt[0]],[pt[1]]], dtype=np.float32))
        self._missings += 1
    
    def remove_all_missings(self, valid_end : int = 2):
        self._point_ids = self._point_ids[:len(self._point_ids) - self._missings]
        self._points = self._points[:len(self._points) - self._missings]
        ridx = reversed(range(valid_end, len(self._points)+1))
        for i in ridx:
            if all(self._point_ids[j] >= 0 for j in range(i - valid_end, i)):
                self._point_ids = self._point_ids[:i]
                self._points = self._points[:i]
                return True

        self._point_ids = []
        self._points = []
        return False
    
    def estimate(self) -> np.ndarray:
        prediction = self.kf.predict()
        return prediction.reshape(-1)[:2]

    def curvature(self, method = 'spline') -> np.ndarray:
        if len(self._points) < 5:
            return np.zeros(len(self._points))
        data = np.array([[x[0], x[1]] for x in self._points], dtype=np.float32)
        return curvature.calculate_angular_acceleration(data)
    
    def is_valid(self, offsets : GlobalOffsetList, min_frame_number: int, min_speed: float, max_acceleration: float, curvature: str) -> bool:
        if len(self._points) < min_frame_number:
            return False
        
        rx = [x[2] for x,id in zip(self._points, self._point_ids) if id >= 0]
        if len(rx) < 0.5 * len(self._points):
            return False

        if np.count_nonzero(np.array(rx)>100) < 0.6 * len(rx):
            return False

        spacing = []
        xy = []
        for i in range(len(self._points)):
            if self._point_ids[i] >= 0:
                spacing.append(i)
                xy.append(np.array(self._points[i][:2])+offsets.global_offset(self._frame_start+i))
        xy = np.array(xy, dtype=np.float32)
        xmin, ymin = xy.min(axis=0)
        xmax, ymax = xy.max(axis=0) 
        xlen = np.ceil(xmax - xmin)
        ylen = np.ceil(ymax - ymin)

        if len(xy) > 2*max(xlen, ylen):
            return False

        v = np.gradient(xy, spacing, axis=0)
        vn = np.linalg.norm(v, axis=1)
        if np.count_nonzero(vn < min_speed) > 0.3 * len(vn):
            return False

        self._curvature = self.curvature(curvature)
        
        return True

def find_seed(  trackerid, seedidx, neighbors, occupied):
    for id in seedidx:
        if id < len(neighbors) and not occupied[id]:
            tid = np.argmax(neighbors[id] >= trackerid)
            if tid < len(neighbors[id]) and neighbors[id][tid] == trackerid:
                return id

    return -1

# def tracking(seeds: List[np.ndarray], trackers: List[KalmanTracker], frameid: int, params: dict, ):
def pointwise_tracking(seeds : list, trackers : list, frameid: int, params : dict, num_nn: int = 3):
    """Track objects using Kalman filter and KD-tree"""
    occupied = [False] * len(seeds)
    
    if len(trackers) > 0:
        predicted_points = [tracker.estimate() for tracker in trackers]
    
        if len(seeds) == 0:
            for tracker, pt in zip(trackers,predicted_points):
                tracker.missing(pt)
            return trackers
    
        # Build KD-tree
        kdtree_trackers = KDTree(np.array(predicted_points))
        kdtree_seeds = KDTree(np.array(seeds)[:,:2])
    
        seed_neighbors = [kdtree_trackers.query(seed[:2], k=num_nn, distance_upper_bound=params['search_radius']*3)[1] for seed in seeds]
        tracker_neighbors = [kdtree_seeds.query(predicted_points[i], k=num_nn, distance_upper_bound=params['search_radius']*3)[1] for i in range(len(trackers))]
    
        for i, pts in enumerate(tracker_neighbors):
            flag = True 
            for id in pts:
                if id < len(seeds) and not occupied[id]:
                    tid = np.argmax(seed_neighbors[id] >= i)
                    if tid < num_nn and seed_neighbors[id][tid] == i:
                        if trackers[i].update(id, seeds[id]):
                            occupied[id] = True
                            flag = False
                            break
            if flag:
                trackers[i].missing(predicted_points[i])
    
    # Add new trackers for unused seeds
    for i, used in enumerate(occupied):
        if not used:
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
                f.write('{}\t{}\t{}\n'.format( frames[tracker._frame_start+i], '\t'.join(f'{x:.2f}' for x in pt), tracker._point_ids[i]))
    except Exception as e:
        print(f"[ERROR]: {e}")

    return [tracker._frame_start, len(tracker._points), np.max(tracker._curvature[0]), file]

def load_tracking_header(f):
    try:
        if isinstance(f, str) :
            if os.path.isfile(f):
                file = open(f, 'r')
                f = file.readline()
        else:
            f = f.readline()

        t = f.strip().split()
        return [int(t[0]), int(t[1])] + [float(x) for x in t[2:]]
    except Exception as e:
        print(f"[ERROR]: {e}")
    return None

def load_tracking( file:str ):
    hdrs = []
    points = []
    try:
        with open(file, 'r') as f:
            lines = f.readlines()
            hdrs = load_tracking_header(lines[0])
            for i in range(hdrs[1]):
                t = lines[i+1].strip().split()
                t[1:-1] = [float(x) for x in t[1:-1]]
                t[-1] = int(t[-1])
                points.append(t)
    except Exception as e:
        print(f"[ERROR]: {e}")
    return hdrs, points

def compute_geometric_center(image_window):
    """
    Compute the geometric center (centroid) of a grayscale image window.

    Args:
        image_window (numpy.ndarray): 2D array representing the image window.

    Returns:
        tuple: (x_center, y_center) coordinates of the centroid.
    """
    # Get the dimensions of the image window
    height, width = image_window.shape

    # Create coordinate grids
    y_coords, x_coords = np.mgrid[0:height, 0:width]

    # Compute total intensity (sum of all pixel values)
    total_intensity = np.sum(image_window)

    if total_intensity == 0:
        # Avoid division by zero; return center of window if uniform
        return width / 2, height / 2

    # Compute weighted centroids
    x_center = np.sum(x_coords * image_window) / total_intensity
    y_center = np.sum(y_coords * image_window) / total_intensity

    return x_center, y_center

def interpolate(image, x, y):
    height, width = image.shape
    x0 = np.floor(x).astype(int)
    x1 = np.minimum(x0+1, width-1)
    y0 = np.floor(y).astype(int)
    y1 = np.minimum(y0+1, height-1)

    x_frac = x-x0
    y_frac = y-y0

    top = image[y0, x0]*(1-x_frac)+image[y0, x1]*x_frac
    bottom = image[y1, x0]*(1-x_frac)+image[y1, x1]*x_frac
    return top*(1-y_frac)+bottom*y_frac

def refine_trajectory(trajectory: list, directory: str, config : dict):
    sz = int(7)
    half_sz = sz // 2
    pdf_size = 3
    margin = (sz-pdf_size) // 2
    flattened_images = []
    for i, pt in enumerate(trajectory):
        frame = os.path.join(directory, pt[0])
        try:
            img = rasterio.open(frame)
            x = int(pt[1])-half_sz
            y = int(pt[2])-half_sz
            if x >= 0 and y >= 0:
                xe = x+sz
                ye = y+sz
                if xe <= img.width and ye <= img.height:
                    win = rasterio.windows.Window( x, y, xe-x, ye-y)
                    data = img.read(1, window=win)
                    xo, yo = compute_geometric_center(data[margin:-margin, margin:-margin])
                    xo += margin
                    yo += margin
                    vmean = np.mean(data)
                    trajectory[i] = [pt[0], x+xo, y+yo, interpolate(data, xo, yo)-vmean, pt[-1]]
                    flattened_images.append(data.flatten())
                    continue
            data = img.read(1, window = rasterio.windows.Window(int(pt[1]),int(pt[2]), 1, 1) )
            trajectory[i] = [pt[0], pt[1], pt[2], data[0,0], pt[-1]]
        except:
            continue
    if len(flattened_images) < 3:
        return trajectory, [0,0]

    image_matrix = np.vstack(flattened_images)

    # Perform PCA
    pca = PCA()
    pca.fit(image_matrix) # Fit PCA on our image data

    # Calculate the cumulative explained variance ratio
    # This tells us how much variance is captured by the first N components.
    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)

    # Find the number of components needed for 95% variance
    num_components_for_95 = np.argmax(cumulative_variance >= 0.95) + 1

    return trajectory, [pca.n_components_, num_components_for_95]
