import numpy as np
from scipy.ndimage import median_filter
from sklearn.decomposition import PCA

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