import numpy as np
from scipy.ndimage import median_filter
from sklearn.decomposition import PCA
import cv2

def extract_background(frames, method='median', pca_components=3):
    """
    Extract background from multiple frames using different methods
    
    Parameters:
    frames -- Input video/hyperspectral sequence ( H x W X T)
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
    H, W, T = frames.shape
    dimT = 2
    
    if method == 'median':
        # Temporal median - robust to transient anomalies
        background = np.median(frames, axis=dimT)
        
    elif method == 'mean':
        # Simple temporal mean
        background = np.mean(frames, axis=dimT)
        
    elif method == 'pca':
        # PCA-based background modeling
        pca = PCA(n_components=pca_components)
        flattened = frames.reshape(-1, T)  # Flatten spatial dimensions
        
        # Fit PCA and reconstruct background
        pca.fit(flattened)
        reconstructed = pca.inverse_transform(pca.transform(flattened))
        background = np.median(reconstructed, axis=1).reshape(H, W)
        
    elif method == 'morph':
        # Morphological approach (median filtering in space and time)
        # First apply temporal median
        temp_median = np.median(frames, axis=dimT)

        # Then spatial median filtering
        background = median_filter(temp_median, size=3)
    
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'median', 'mean', 'pca', or 'morph'")
    
    return background

def remove_background(frames, interval = 30, method = 'median'):
    """
    Remove background from an image using the estimated background

    Parameters:
    frames -- Input video sequence (H x W X T)

    Returns:
    foreground -- Image with background removed (H x W X T)
    """
    T = frames.shape[2]

    if interval is None:
        interval = T

    if interval > T:
        raise ValueError("Interval must be less than or equal to the number of frames")

    foreground = np.zeros_like(frames)

    for i in range(0, T, interval):

        ed = min(i + interval, T)
        background = extract_background(frames[ ..., i:ed], method=method)

        foreground[  ..., i:ed ] = frames[  ..., i:ed] - background[  ..., np.newaxis ]

    return foreground

class BackgroundSubtractor:
    def __init__(self, history=10, var_threshold=2.5, n_mixtures = 5, detect_shadows=False, warmup_frames=None):
        
        # Initialize MOG2
        self.mog = cv2.createBackgroundSubtractorMOG2(
            history=history,
            varThreshold=var_threshold,
            detectShadows=detect_shadows
        )

        self.mog.setNMixtures(n_mixtures)
        
        # State management
        self.is_initialized = False
        self.warmup_frames = warmup_frames

    def warmup(self, frames, warmup_learning_rate=0.8):
        """Warm-up MOG2 with initial frames"""
        dtype = frames[0].dtype
        if self.warmup_frames is not None:
            frames = frames[:self.warmup_frames]
        for frame in frames:
            self.mog.apply(frame if dtype == np.uint8 or dtype == np.float32 else frame.astype(np.float32), learningRate=warmup_learning_rate)
        self.is_initialized = True
    
    def process(self, frame, learning_rate=0.1):
        """
        Process a single frame through the pipeline
        
        Args:
            frame: Input frame (BGR or grayscale)
            
        Returns:
            foreground_mask: Binary foreground mask
            background_model: Current background model
        """
        
        dtype = frame.dtype
        
        # Apply MOG2
        foreground_mask = self.mog.apply(frame if dtype == np.uint8 or dtype == np.float32 else frame.astype(np.float32), learningRate=learning_rate)
        
        # Get background model
        background_model = self.mog.getBackgroundImage()
        
        return foreground_mask, background_model.astype(dtype)

if __name__ == '__main__':
    import glob, os
    image_files = sorted(glob.glob("/Users/xlingsky/Desktop/demo/FJ12/*.tif"))
    frames = [cv2.imread(f, cv2.IMREAD_UNCHANGED) for f in image_files]

    bg_subtractor = BackgroundSubtractor(warmup_frames=10)

    output_dir = '/Users/xlingsky/Desktop/demo/test'
    os.makedirs(os.path.join(output_dir, 'fg'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'mask'), exist_ok=True)

    for path, frame in zip( image_files, frames):
        filename = os.path.splitext( os.path.basename(path) )[0]
        fg, bg = bg_subtractor.process(frame)
        cv2.imwrite( os.path.join(output_dir, 'mask', f'{filename}_mask.tif'), fg)
        cv2.imwrite( os.path.join(output_dir, 'fg', f'{filename}_fg.tif'), frame.astype(np.int16) - bg.astype(np.int16) )

