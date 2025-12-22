import cv2
import numpy as np
from scipy import ndimage
from scipy.spatial import KDTree
from sklearn.cluster import DBSCAN
from collections import deque
from hsp.utils.common import rasterio_write, rasterio_read
from hsp.utils.detector import Detector
from dataclasses import dataclass
from enum import Enum
import warnings
from typing import Dict, List, Tuple, Optional
warnings.filterwarnings('ignore')

@dataclass
class TrackingConfig:
    
    # Template matching
    template_extending_size: int = 5
    search_margin: int = 25
    correlation_threshold: float = 0.8
    seed_correlation_threshold: float = 0.6
    distance_to_merge: int = 2
    
    # Tracking management
    seed_confidence_threshold: float = 0.1
    max_missing_frames: int = 15
    min_track_length: int = 10
    new_track_threshold: float = 0.1
    max_tracks: int = 10000

    # State constraints (optional)
    min_velocity: float = 0.0    # pixels/frame
    max_velocity: float = 100.0  # pixels/frame
    max_acceleration: float = 2.0  # pixels/frame²
    std_innovation: float = 2  # pixels²
    
@dataclass
class KalmanConfig:
    """Configuration for OpenCV Kalman Filter"""
    # Motion model type
    motion_model: str = 'constant_velocity'  # 'constant_velocity', 'constant_acceleration'
    
    # Noise parameters (higher values = more uncertainty)
    process_noise_pos: float = 1e-3      # Process noise for position
    process_noise_vel: float = 1e-2      # Process noise for velocity
    process_noise_acc: float = 1e-1      # Process noise for acceleration
    measurement_noise: float = 1e-2      # Measurement noise
    
    # Time step (seconds per frame)
    dt: float = 1.0
    
    # Adaptation for small objects
    adapt_for_small_objects: bool = False
    small_object_max_size: int = 32  # pixels
    adaptive_noise_scaling: float = 2.0  # Scale noise for small objects
    
class Tracker:
    def __init__(self, config: Optional[KalmanConfig] = None):
        self.config = config or KalmanConfig()
        
       # Initialize OpenCV KalmanFilter
        self._init_kalman()
        
        # Track state
        self.is_initialized = False
        self.start_frame = 0
        self.last_prediction: Optional[Tuple[float, float]] = None
        self.last_update: Optional[Tuple[float, float]] = None
        
        # History for analysis
        self.history: List[Tuple[float, float]] = []
        self.seed_id : List[int] = []
        self.seed_history: List[Dict] = []
        self.innovation_history: List[Tuple[float, float]] = []
        
        # For adaptive noise
        self.object_size: Optional[np.ndarray] = None
        self.missing_frames = 0
        self.confidence = 1.0 

    def _init_kalman(self):
        """Initialize OpenCV KalmanFilter based on motion model"""
        if self.config.motion_model == 'constant_velocity':
            # State: [x, y, vx, vy] (4)
            # Measurement: [x, y] (2)
            self.kf = cv2.KalmanFilter(4, 2)
            
            # State transition matrix (F)
            dt = self.config.dt
            self.kf.transitionMatrix = np.array([
                [1, 0, dt, 0],
                [0, 1, 0, dt],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ], dtype=np.float32)
            
            # Measurement matrix (H)
            self.kf.measurementMatrix = np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0]
            ], dtype=np.float32)
            
            # Process noise covariance (Q)
            self.kf.processNoiseCov = np.diag([
                self.config.process_noise_pos,
                self.config.process_noise_pos,
                self.config.process_noise_vel,
                self.config.process_noise_vel
            ]).astype(np.float32)
            
            # Measurement noise covariance (R)
            self.kf.measurementNoiseCov = np.eye(2, dtype=np.float32) * self.config.measurement_noise
            
            # Initial state
            self.kf.statePre = np.zeros((4, 1), dtype=np.float32)
            self.kf.statePost = np.zeros((4, 1), dtype=np.float32)
            
            # Error covariance (P)
            self.kf.errorCovPre = np.eye(4, dtype=np.float32)
            self.kf.errorCovPost = np.eye(4, dtype=np.float32) * 10.0
            
        elif self.config.motion_model == 'constant_acceleration':
            # State: [x, y, vx, vy, ax, ay] (6)
            # Measurement: [x, y] (2)
            self.kf = cv2.KalmanFilter(6, 2)
            
            dt = self.config.dt
            dt2 = 0.5 * dt * dt
            
            # State transition matrix (F)
            self.kf.transitionMatrix = np.array([
                [1, 0, dt, 0, dt2, 0],
                [0, 1, 0, dt, 0, dt2],
                [0, 0, 1, 0, dt, 0],
                [0, 0, 0, 1, 0, dt],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1]
            ], dtype=np.float32)
            
            # Measurement matrix (H)
            self.kf.measurementMatrix = np.array([
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0]
            ], dtype=np.float32)
            
            # Process noise covariance (Q)
            self.kf.processNoiseCov = np.diag([
                self.config.process_noise_pos,
                self.config.process_noise_pos,
                self.config.process_noise_vel,
                self.config.process_noise_vel,
                self.config.process_noise_acc,
                self.config.process_noise_acc
            ]).astype(np.float32)
            
            # Measurement noise covariance (R)
            self.kf.measurementNoiseCov = np.eye(2, dtype=np.float32) * self.config.measurement_noise
            
            # Initial state
            self.kf.statePre = np.zeros((6, 1), dtype=np.float32)
            self.kf.statePost = np.zeros((6, 1), dtype=np.float32)
            
            # Error covariance (P)
            self.kf.errorCovPre = np.eye(6, dtype=np.float32)
            self.kf.errorCovPost = np.eye(6, dtype=np.float32) * 10.0
        
        else:
            raise ValueError(f"Unknown motion model: {self.config.motion_model}")
        
    def initialize(self, 
                  start_frame: int, seed_id: int, template: np.ndarray, seed: dict,
                  initial_position: Tuple[float, float],
                  initial_velocity: Optional[Tuple[float, float]] = None,
                  initial_size: Optional[np.ndarray] = None):
        """
        Initialize the Kalman filter with initial state
        
        Args:
            initial_position: (x, y) initial position in pixels
            initial_velocity: (vx, vy) initial velocity in pixels/frame (optional)
            initial_size: (width, height) object size for adaptive noise (optional)
        """
        # Store object size for adaptive noise
        self.object_size = np.array(seed['bbox'][2:], dtype=np.float32) if initial_size is None else initial_size
        
        # Set initial state based on motion model
        if self.config.motion_model == 'constant_velocity':
            if initial_velocity is not None:
                self.kf.statePre = np.array([
                    [initial_position[0]],
                    [initial_position[1]],
                    [initial_velocity[0]],
                    [initial_velocity[1]]
                ], dtype=np.float32)
            else:
                self.kf.statePre = np.array([
                    [initial_position[0]],
                    [initial_position[1]],
                    [0.0],
                    [0.0]
                ], dtype=np.float32)
            
            self.kf.statePost = self.kf.statePre.copy()
            
            # Set initial covariance higher for unknown velocity
            if initial_velocity is None:
                self.kf.errorCovPost[2, 2] = 100.0
                self.kf.errorCovPost[3, 3] = 100.0
        
        elif self.config.motion_model == 'constant_acceleration':
            if initial_velocity is not None:
                self.kf.statePre = np.array([
                    [initial_position[0]],
                    [initial_position[1]],
                    [initial_velocity[0]],
                    [initial_velocity[1]],
                    [0.0],  # acceleration
                    [0.0]
                ], dtype=np.float32)
            else:
                self.kf.statePre = np.array([
                    [initial_position[0]],
                    [initial_position[1]],
                    [0.0],
                    [0.0],
                    [0.0],
                    [0.0]
                ], dtype=np.float32)
            
            self.kf.statePost = self.kf.statePre.copy()
            
            # Higher covariance for unknown velocity/acceleration
            if initial_velocity is None:
                self.kf.errorCovPost[2, 2] = 100.0
                self.kf.errorCovPost[3, 3] = 100.0
                self.kf.errorCovPost[4, 4] = 100.0
                self.kf.errorCovPost[5, 5] = 100.0
        
        # Apply adaptive noise if enabled and object size is known
        if self.config.adapt_for_small_objects and self.object_size is not None:
            self._apply_adaptive_noise()
        
        self.start_frame = start_frame
        self.template = template

        # Initialize history
        self.history = [initial_position]
        self.seed_id = [seed_id]
        self.seed_history = [{'id': seed_id, 'seed': seed, 'distance': 0.0, 'confidence': 1.0}]
        self.innovation_history = [(0.0, 0.0)]
        
        self.is_initialized = True
        self.missing_frames = 0
        self.confidence = 1.0
        self.last_prediction = initial_position
        self.last_update = initial_position
    
    def predict(self) -> Tuple[float, float]:
        """
        Predict the next state
        
        Returns:
            Predicted position (x, y)
        """
        # Predict next state
        predicted_state = self.kf.predict()
        
        # Get predicted position
        predicted_position = (float(predicted_state[0]), float(predicted_state[1]))
        
        # Store prediction
        self.last_prediction = predicted_position
        
        return predicted_position
    
    def update(self, measurement: Tuple[float, float], id : int, info: dict, template : Optional[np.ndarray] = None) -> Tuple[float, float]:
        # Prepare measurement vector
        measurement_np = np.array([[measurement[0]], [measurement[1]]], dtype=np.float32)
        
        # Update the filter
        corrected_state = self.kf.correct(measurement_np)
        
        # Get corrected position
        corrected_position = (float(corrected_state[0]), float(corrected_state[1]))
        
        # Calculate innovation (measurement - prediction)
        if self.last_prediction is not None:
            innovation = (
                measurement[0] - self.last_prediction[0],
                measurement[1] - self.last_prediction[1]
            )
            self.innovation_history.append(innovation)
        
        # Store update
        self.last_update = corrected_position
        self.history.append(measurement)
        self.seed_id.append(id)
        
        if id >= -1:
            self.missing_frames = 0
            self.confidence = 1.0
        else:
            # Update missing frames counter
            self.missing_frames += 1

            # Decay confidence when predicting without updates
            if self.missing_frames > 0:
                self.confidence *= 0.95
        
        # Adapt noise based on innovation if enabled
        if self.config.adapt_for_small_objects:
            self._adapt_noise_based_on_innovation()

        if template is not None:
            self.template = template

        self.seed_history.append(info)

        self.object_size = (self.object_size+np.array(info['seed']['bbox'][2:]))/2
        
        return corrected_position
    
    def get_velocity(self) -> Tuple[float, float]:
        """Get current estimated velocity"""
        return (float(self.kf.statePost[2]), float(self.kf.statePost[3]))

    def _apply_adaptive_noise(self):
        """Adapt noise parameters for small objects"""
        if self.object_size is None:
            return
        
        width, height = self.object_size
        max_dimension = max(width, height)
        
        if max_dimension < self.config.small_object_max_size:
            # Increase process noise for small objects (they move more erratically)
            scale_factor = self.config.adaptive_noise_scaling * \
                         (self.config.small_object_max_size / max(1, max_dimension))
            
            # Scale process noise
            if self.config.motion_model == 'constant_velocity':
                self.kf.processNoiseCov[0, 0] *= scale_factor
                self.kf.processNoiseCov[1, 1] *= scale_factor
                self.kf.processNoiseCov[2, 2] *= scale_factor
                self.kf.processNoiseCov[3, 3] *= scale_factor
            else:  # constant_acceleration
                self.kf.processNoiseCov[0, 0] *= scale_factor
                self.kf.processNoiseCov[1, 1] *= scale_factor
                self.kf.processNoiseCov[2, 2] *= scale_factor
                self.kf.processNoiseCov[3, 3] *= scale_factor
                self.kf.processNoiseCov[4, 4] *= scale_factor
                self.kf.processNoiseCov[5, 5] *= scale_factor

    def get_uncertainty(self) -> float:
        """Get position uncertainty (trace of position covariance)"""
        return float(self.kf.errorCovPost[0, 0] + self.kf.errorCovPost[1, 1])
    
    def get_prediction_uncertainty(self) -> float:
        """Get prediction uncertainty"""
        return float(self.kf.errorCovPre[0, 0] + self.kf.errorCovPre[1, 1])

    def get_acceleration(self) -> Tuple[float, float]:
        """Get current estimated acceleration (if applicable)"""
        if self.config.motion_model == 'constant_acceleration':
            return (float(self.kf.statePost[4]), float(self.kf.statePost[5]))
        else:
            return (0.0, 0.0)

    def _adapt_noise_based_on_innovation(self):
        """Adapt noise based on recent innovation history"""
        if len(self.innovation_history) < 3:
            return
        
        # Calculate innovation statistics
        innovations = np.array(self.innovation_history[-10:])  # Last 10 innovations
        innovation_norms = np.sqrt(np.sum(innovations**2, axis=1))
        
        if len(innovation_norms) < 3:
            return
        
        mean_innovation = np.mean(innovation_norms)
        std_innovation = np.std(innovation_norms)
        
        # Adjust measurement noise based on innovation consistency
        if std_innovation > mean_innovation * 2.0:
            # Inconsistent measurements, increase measurement noise
            self.kf.measurementNoiseCov *= 1.5
        elif std_innovation < mean_innovation * 0.5:
            # Very consistent, decrease measurement noise
            self.kf.measurementNoiseCov *= 0.8
        
        # Ensure noise doesn't go below minimum
        min_noise = 1e-4
        self.kf.measurementNoiseCov = np.maximum(self.kf.measurementNoiseCov, min_noise)

    def get_performance_metrics(self) -> dict:
        """Get performance metrics for monitoring"""
        if not self.is_initialized or len(self.innovation_history) < 2:
            return {}
        
        innovations = np.array(self.innovation_history[1:])  # Skip first zero
        innovation_norms = np.sqrt(np.sum(innovations**2, axis=1))
        
        return {
            'avg_innovation': float(np.mean(innovation_norms)),
            'std_innovation': float(np.std(innovation_norms)),
            'max_innovation': float(np.max(innovation_norms)),
            'missing_frames': self.missing_frames,
            'confidence': self.confidence,
            'position_uncertainty': self.get_uncertainty(),
            'prediction_uncertainty': self.get_prediction_uncertainty(),
            'total_frames': len(self.history)
        }

    def get_template(self) -> np.ndarray:
        return self.template

    def finish(self):
        if self.missing_frames == 0:
            return
        num_frames = len(self.history)
        self.history = self.history[:num_frames-self.missing_frames]
        self.seed_id = self.seed_id[:num_frames-self.missing_frames]
        self.seed_history = self.seed_history[:num_frames-self.missing_frames]
        self.innovation_history = self.innovation_history[:num_frames-self.missing_frames]
        self.missing_frames = 0

    def valid_frame_number(self) -> int:
        return np.count_nonzero(np.array(self.seed_id) >= 0)

    def save(self, filepath:str, frames:list[str]):
        try:
            with open(filepath, 'w') as f:
                metrics = self.get_performance_metrics()
                f.write(f'{self.start_frame} {len(self.history)} {metrics["avg_innovation"]:.2f} {metrics["std_innovation"]:.2f}\n')#{self.confidence:.2f}
                for i, info in enumerate(self.seed_history):
                    seed = info['seed']
                    f.write(f'{frames[self.start_frame+i]}\t{seed['centroid'][0]:.2f}\t{seed['centroid'][1]:.2f}\t{seed['intensity']:.1f}\t{seed['area']:.0f}\t{info['id']}\n')
        except Exception as e:
            print(f"Error saving tracker: {e}")

class Manager:
    def __init__(self, start_frame: int, config: Optional[TrackingConfig] = None):
        self.config = config or TrackingConfig()
        
        # Tracking state
        self.tracks: List[Tracker] = []
        self.history_tracks: List[Tracker] = []
        self.frame_count = 0
        self.frame_id = start_frame
        
        # Camera motion
        self.camera_translation_history: List[np.ndarray] = []
        
        # Statistics
        self.stats = {
            'seeds_detected': 0,
            'tracks_created': 0,
            'tracks_lost': 0,
            'matches_found': 0,
            'processing_time': 0.0
        }

    def _extract_template(self, 
                        frame: np.ndarray,
                        seed: Dict) -> np.ndarray:
        """Extract template from frame"""
        center = seed['centroid']
        size = (seed['bbox'][2]+self.config.template_extending_size, seed['bbox'][3]+self.config.template_extending_size)
        x, y = center
        w, h = size
        
        x1 = max(0, int(x - w // 2))
        y1 = max(0, int(y - h // 2))
        x2 = min(frame.shape[1], int(x + w // 2))
        y2 = min(frame.shape[0], int(y + h // 2))
        
        template = frame[y1:y2, x1:x2].copy()
        
        # Pad if needed
        if template.shape[0] < h or template.shape[1] < w:
            pad_y = (h - template.shape[0]) // 2
            pad_x = (w - template.shape[1]) // 2
            padded = cv2.copyMakeBorder(
                template, pad_y, pad_y, pad_x, pad_x,
                cv2.BORDER_REFLECT_101
            )
            return padded
        
        return template

    def _match_template(self,
                      frame: np.ndarray,
                      template: np.ndarray,
                      search_center: Tuple[float, float],
                      search_radius: int,
                      scales: List[float] = [1.0]
                      ) -> Dict:
        """Template matching with multi-scale support"""
        best_match = {
            'position': search_center,
            'confidence': 0.0,
            'scale': 1.0
        }
        
        for scale in scales:
            # Scale template
            if scale != 1.0:
                scaled_template = cv2.resize(template.astype(np.float32), None,
                                           fx=scale, fy=scale,
                                           interpolation=cv2.INTER_AREA)
                scaled_template = scaled_template.astype(np.uint16)
            else:
                scaled_template = template
            
            t_h, t_w = scaled_template.shape
            
            # Define search region
            search_x, search_y = search_center
            x1 = max(0, int(search_x - search_radius))
            y1 = max(0, int(search_y - search_radius))
            x2 = min(frame.shape[1], int(search_x + search_radius + t_w))
            y2 = min(frame.shape[0], int(search_y + search_radius + t_h))
            
            if x2 <= x1 or y2 <= y1:
                continue
            
            search_region = frame[y1:y2, x1:x2]
            
            if search_region.shape[0] < t_h or search_region.shape[1] < t_w:
                continue
            
            # Template matching
            try:
                result = cv2.matchTemplate(search_region.astype(np.float32), scaled_template.astype(np.float32), 
                                          cv2.TM_CCOEFF_NORMED)
                min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(result)
                
                if max_val > best_match['confidence']:
                    match_x = max_loc[0] + x1 + t_w // 2
                    match_y = max_loc[1] + y1 + t_h // 2
                    
                    best_match = {
                        'position': (float(match_x), float(match_y)),
                        'confidence': float(max_val),
                        'scale': scale
                    }
            except Exception as e:
                continue
        
        return best_match

    def create_track_from_seed(self, 
                              seed_id: int,
                              seed: Dict,
                              frame: np.ndarray) -> Tracker:
        """Create a new track from a seed detection"""
        
        # Extract template
        template = self._extract_template(frame, seed)
        
        # Create track
        track = Tracker()
        track.initialize( self.frame_id, seed_id, template, seed=seed, initial_position=seed['centroid'] )
        
        # Store
        self.tracks.append(track)
        
        self.stats['tracks_created'] += 1
        
        return track

    def _query_seed(self, seeds: List[Dict], seedtree: KDTree, position: Tuple[float, float], seed : Dict, distance_upper_bound : float,  k=3, min_confidence=0.5) -> Tuple[float, int]:
        if seedtree is None:
            return min_confidence, -1
        _, indices = seedtree.query(position, k=k, distance_upper_bound=distance_upper_bound)
        best_id = -1
        for id in indices:
            if id >= len(seeds):
                break
            intensity0 = seed['intensity']
            intensity1 = seeds[id]['intensity']
            confidence = 1 - abs(intensity1-intensity0)/max(intensity0, intensity1)
            if confidence > min_confidence:
                min_confidence = confidence
                best_id = id
        return min_confidence, best_id

    def update_tracks(self,
                      seeds: List[Dict],
                      tracks: List[Tracker],
                      frame: np.ndarray,
                      camera_motion: np.ndarray = np.array([0.0,0.0])) -> Tuple[Dict, List]:
        """Associate detected seeds to existing tracks"""
        associations = {}

        seed_occupied = [False] * len(seeds)
        kdtree_seeds = KDTree(np.array([np.array(seed['centroid']) for seed in seeds])) if len(seeds) > 0 else None

        for track_id, track in enumerate(tracks):
            
            # Predict position
            predicted_pos = track.predict()
            predicted_frame_pos = (predicted_pos[0] - camera_motion[0], predicted_pos[1] - camera_motion[1])
            confidence, best_seed_idx = self._query_seed(seeds, kdtree_seeds, predicted_frame_pos, track.seed_history[-1]['seed'], self.config.distance_to_merge+np.max(track.object_size))
            # kdtree_seeds.query(predicted_frame_pos) if kdtree_seeds is not None else (float('inf'), -1)
            if best_seed_idx >=0 :
                associations[track_id] = {
                    'id': best_seed_idx,
                    'seed': seeds[best_seed_idx],
                    'distance': 0,
                    'confidence': confidence #self.config.correlation_threshold
                }
                seed_occupied[best_seed_idx] = True
                continue
            
            # Find best matching seed
            best_seed_idx = -1
            best_match_confidence = 0.0

            match_result = self._match_template(
                frame,
                track.get_template(),
                predicted_frame_pos,
                (self.config.search_margin + max(track.get_template().shape[0], track.get_template().shape[1]) )// 2,
                scales=[1.0]
            )

            if match_result['confidence'] > self.config.seed_correlation_threshold:
                confidence, best_seed_idx = self._query_seed(seeds, kdtree_seeds, match_result['position'], track.seed_history[-1]['seed'], self.config.distance_to_merge+np.max(track.object_size))
                # dist, best_seed_idx = kdtree_seeds.query(match_result['position']) if kdtree_seeds is not None else (float('inf'), -1)
                if best_seed_idx >= 0:
                    best_match_confidence = match_result['confidence']
                    seed_occupied[best_seed_idx] = True
                    associations[track_id] = {
                        'id': best_seed_idx,
                        'seed': seeds[best_seed_idx],
                        'distance': 0,
                        'confidence': best_match_confidence
                    }
                elif match_result['confidence'] > self.config.correlation_threshold:
                    last_seed = track.seed_history[-1]['seed'].copy()
                    last_seed['centroid'] = match_result['position']
                    associations[track_id] = {
                        'id': -1,
                        'seed': last_seed,
                        'distance': 0,
                        'confidence': match_result['confidence']
                    }
        return associations, seed_occupied

    def process(self, seeds: List[Dict], frame: np.ndarray, camera_motion: np.ndarray = np.array([0.0,0.0]) ) -> bool:
        self.frame_count += 1
        
        self.stats['seeds_detected'] = len(seeds)
        
        self.camera_translation_history.append(camera_motion)
        
        filtered_seeds = [
            seed for seed in seeds 
            if seed['confidence'] >= self.config.seed_confidence_threshold
        ]
        
        associations, seed_occupied = self.update_tracks(
            filtered_seeds, self.tracks, frame, camera_motion
        )

        seed_to_tracks = {}
        for track_id, association in associations.items():
            if association['id'] == -1:
                continue
            if association['id'] in seed_to_tracks:
                seed_to_tracks[association['id']].append(track_id)
            else:
                seed_to_tracks[association['id']] = [track_id]

        for _, track_ids in seed_to_tracks.items():
            if len(track_ids) < 2:
                continue
            best_id, best_confidence = -1, 0.0
            for track_id in track_ids:
                confidence = associations[track_id]['confidence'] + len(self.tracks[track_id].history)/self.config.min_track_length
                if  confidence > best_confidence:
                    best_id = track_id
                    best_confidence = confidence
            for track_id in track_ids:
                if track_id != best_id:
                    associations.pop(track_id)

        missing_tracks = [True] * len(self.tracks)
        for track_id, association in associations.items():
            track = self.tracks[track_id]
            track.update(
                (association['seed']['centroid'][0] + camera_motion[0], association['seed']['centroid'][1] + camera_motion[1]), 
                association['id'], 
                info=association, template=self._extract_template(frame, association['seed'])
                )
            missing_tracks[track_id] = False
        

        for track_id, track in enumerate(self.tracks):
            if missing_tracks[track_id]:
                predicted_pos = track.last_prediction
                info = {
                    'id': -2,
                    'seed': track.seed_history[-1]['seed'].copy(),
                    'distance': 0,
                    'confidence': 0
                }
                info['seed']['centroid'] = (predicted_pos[0]-camera_motion[0], predicted_pos[1]-camera_motion[1])
                track.update( (predicted_pos[0], predicted_pos[1]), -2, info=info, template=None)

        for i, seed in enumerate(filtered_seeds):
            if seed_occupied[i]:
                continue
            # Check if we should create new track
            if seed['confidence'] >= self.config.new_track_threshold:
                # Create new track
                self.create_track_from_seed(i, seed, frame)
        
        self.frame_id += 1
        return True

    def _evaluate_track_quality(self, track: Tracker) -> float:
        valid_frame_number = track.valid_frame_number() 
        if track.missing_frames > self.config.max_missing_frames:
            return 0.0
        if len(track.history) >= min(self.config.min_track_length//2, 5):
            if valid_frame_number/ len(track.history) < 0.5:
                return 0.0
            speed = np.linalg.norm(np.array(track.history[-1]) - np.array(track.history[0])) / (len(track.history)-1)
            if speed < self.config.min_velocity or speed > self.config.max_velocity:
                return 0.0
            speed = np.linalg.norm(np.array(track.get_velocity()))
            if speed < self.config.min_velocity or speed > self.config.max_velocity:
                return 0.0
        # if 2*self.config.min_track_length >= len(track.history) >= self.config.min_track_length :
        #     metrics = track.get_performance_metrics()
        #     if len(metrics)>0 and metrics['max_innovation'] > 10:
        #         return 0.0
        return valid_frame_number / self.config.min_track_length

    def check(self) -> Tuple[List[Tracker], List[Tracker]]:
        """Check and remove lost tracks"""
        active_tracks = []
        terminated_tracks = []
        new_tracks = []
        for track in self.tracks:
            quality = self._evaluate_track_quality(track)
            if quality > 0:
                active_tracks.append(track)
                if quality == 1.0:
                    new_tracks.append(track)
            else:
                terminated_tracks.append(track)

        last_history_count = len(self.history_tracks)
        for track in terminated_tracks:
            track.finish()
            if self._evaluate_track_quality(track) >= 1.0:
                self.history_tracks.append(track)
            self.stats['tracks_lost'] += 1
        self.tracks = active_tracks
        
        return new_tracks, self.history_tracks[last_history_count:]

    def terminate(self):
        last_history_count = len(self.history_tracks)
        for track in self.tracks:
            track.finish()
            if self._evaluate_track_quality(track) >= 1.0:
                self.history_tracks.append(track)
        return self.history_tracks[last_history_count:]

class DenoisingModule:
    """Advanced denoising module for noisy frames"""
    
    def __init__(self, config):
        self.config = config
        
    def process(self, frame):
        """Apply multiple denoising techniques"""

        config = self.config['gaussian_blur']
        if config['apply']:
            frame = ndimage.gaussian_filter(
                frame, 
                sigma=config['sigma'] if config['sigma'] else np.array(config['radius'])/config['truncate'],
                mode='mirror',
                truncate=config['truncate']
            )
        
        # # Method 2: Non-local means denoising
        # frame = cv2.fastNlMeansDenoising(
        #     frame,
        #     h=self.config['non_local_means']['h'],
        #     templateWindowSize=self.config['non_local_means']['templateWindowSize'],
        #     searchWindowSize=self.config['non_local_means']['searchWindowSize']
        # )
        
        # # Method 4: Wavelet denoising (simplified)
        # frame = self.wavelet_denoise(frame)
        
        return frame 
    
    def wavelet_denoise(self, image):
        """Simplified wavelet-based denoising"""
        # Discrete wavelet transform approximation
        coeffs = cv2.dct(np.float32(image))
        
        # Soft thresholding
        threshold = self.config['wavelet_threshold'] * np.max(np.abs(coeffs))
        coeffs_thresh = np.sign(coeffs) * np.maximum(np.abs(coeffs) - threshold, 0)
        
        # Inverse transform
        denoised = cv2.idct(coeffs_thresh)
        return np.uint8(np.clip(denoised, 0, 255))

if __name__ == '__main__':
    import os
    import glob
    from sklearn.cluster import DBSCAN
    from skimage import segmentation, filters, feature, morphology, exposure, measure
    from scipy import ndimage as ndi
    import matplotlib.pyplot as plt
    srcdir = '/Users/xlingsky/Desktop/demo/test'
    dstdir = '/Users/xlingsky/Desktop/demo/test'

    image_filess = sorted(glob.glob(os.path.join(srcdir, 'fg') + '/*.tif'))
    image_path = image_filess[1]

    image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    mask = cv2.imread(image_path.replace('fg', 'mask'), cv2.IMREAD_UNCHANGED)
    detector = Detector()
    objects = detector.process(image, mask, max_region_area=25, min_region_area=2)

    mask = np.zeros_like(image, dtype=np.uint8)
    for obj in objects:
        mask[obj.coords[:,0], obj.coords[:,1]] = 255
        
    cv2.imwrite(os.path.join(dstdir, 'segmentation.tif'), mask)
    

    # segments = segmentation.find_boundaries(labels, mode='thick')
    # names = ['image', 'segmentation']

    # image = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    # image = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)

    # color_mask = cv2.applyColorMap(mask, cv2.COLORMAP_JET)
    # cv2.imshow('label', color_mask)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()


    # tracker = MovingObjectTracker()
    # tracker.process(srcdir=srcdir, dstdir=dstdir, start_from=0)