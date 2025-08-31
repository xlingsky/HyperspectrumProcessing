import numpy as np
from scipy.linalg import inv
from background import remove_background
import sys
import os
from osgeo import gdal
import cv2 as cv

def global_rx_detector(image, background_mask=None):
    """
    Global RX (Reed-Xiaoli) anomaly detection algorithm
    
    Parameters:
    image -- Input hyperspectral image (H x W x B), where:
             H = height, W = width, B = number of spectral bands
    background_mask -- Optional binary mask (H x W) specifying background pixels
    
    Returns:
    anomaly_map -- Detection result (H x W) with RX statistic values
    """
    # Get image dimensions
    h, w, b = image.shape
    
    # Use entire image as background if no mask provided
    if background_mask is None:
        background_mask = np.ones((h, w), dtype=bool)
    
    # Reshape image to 2D matrix (N x B), where N = H*W
    pixels = image.reshape(-1, b)
    
    # Extract background pixels
    background_pixels = pixels[background_mask.reshape(-1)]
    
    # Calculate background statistics
    mean = np.mean(background_pixels, axis=0)
    cov = np.cov(background_pixels, rowvar=False)
    
    # Handle potential singular covariance matrix
    if cov.shape == ():
        cov = np.array([[cov]])
    try:
        inv_cov = inv(cov)
    except np.linalg.LinAlgError:
        inv_cov = inv(cov + np.eye(b) * 1e-6)  # Regularization
    
    # Center all pixels
    centered_pixels = pixels - mean
    
    # Compute RX statistic: (x - μ)^T * Σ^-1 * (x - μ)
    rx_values = np.sum(centered_pixels @ inv_cov * centered_pixels, axis=1)
    
    # Reshape results to original image dimensions
    return rx_values.reshape(h, w)

def local_rx_detector(image, window_size=15, inner_window=3):
    """
    Local RX anomaly detection with dual-window approach
    
    Parameters:
    image -- Input hyperspectral image (H x W x B)
    window_size -- Size of outer background window (odd integer)
    inner_window -- Size of inner guard window (odd integer, < window_size)
    
    Returns:
    anomaly_map -- Detection result (H x W)
    """
    h, w, b = image.shape
    anomaly_map = np.zeros((h, w))
    half_win = window_size // 2
    half_inner = inner_window // 2
    
    # Pad image to handle borders
    padded_img = np.pad(image, ((half_win, half_win), 
                                (half_win, half_win), 
                                (0, 0)), mode='reflect')
    
    # Process each pixel with local windows
    for i in range(h):
        for j in range(w):
            # Get outer window coordinates
            i_pad, j_pad = i + half_win, j + half_win
            outer = padded_img[i_pad-half_win:i_pad+half_win+1, 
                              j_pad-half_win:j_pad+half_win+1, :]
            
            # Create mask excluding inner window
            mask = np.ones((window_size, window_size), dtype=bool)
            mask[half_win-half_inner:half_win+half_inner+1, 
                 half_win-half_inner:half_win+half_inner+1] = False
            
            # Extract background pixels
            background = outer[mask].reshape(-1, b)
            
            # Skip if insufficient background pixels
            if len(background) < b:
                anomaly_map[i, j] = 0
                continue
            
            # Calculate local statistics
            mean = np.mean(background, axis=0)
            cov = np.cov(background, rowvar=False)
            
            # Handle singular covariance
            try:
                inv_cov = inv(cov)
            except np.linalg.LinAlgError:
                inv_cov = inv(cov + np.eye(b) * 1e-6)
            
            # Compute RX statistic for center pixel
            center_pixel = np.max(padded_img[i_pad-half_inner:i_pad+half_inner+1, j_pad-half_inner:j_pad+half_inner+1, :], axis=(0, 1))
            centered = center_pixel - mean
            anomaly_map[i, j] = centered @ inv_cov @ centered
    
    return anomaly_map

def visualize_results(original, detection, title="RX Anomaly Detection", output=None, dpi=200):
    """
    Visualize original image and detection results with output control
    
    Parameters:
    original -- Original input image (H x W x B)
    detection -- Detection result map (H x W)
    title -- Title for the detection plot (default: "RX Anomaly Detection")
    output -- Output control: 
              None (default) - shows interactive plot
              'filename.ext' - saves to file (e.g., 'results.png')
    dpi -- Resolution for saved images (default: 100)
    """
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(15, 6))
    
    # Original image (first band)
    plt.subplot(1, 2, 1)
    plt.imshow(original[:, :, 0], cmap='gray')
    plt.title('Original Image (First Band)')
    plt.colorbar()
    
    # Detection results
    plt.subplot(1, 2, 2)
    plt.imshow(detection, cmap='hot')
    plt.title(title)
    plt.colorbar()
    
    plt.tight_layout()
    
    # Handle output
    if output is None:
        plt.show()
    else:
        # Extract file extension
        if '.' not in output:
            output += '.png'  # Default to PNG if no extension
            
        plt.savefig(output, dpi=dpi, bbox_inches='tight')
        plt.close()
        print(f"Results saved to {output}")

def test():
    # Simulate hyperspectral data (100x100x10)
    np.random.seed(42)
    h, w, b = 100, 100, 10
    
    # Create background data (multivariate normal)
    mean = np.zeros(b)
    cov = np.eye(b) * 0.5 + np.ones((b, b)) * 0.5
    background = np.random.multivariate_normal(mean, cov, size=h*w).reshape(h, w, b)
    
    # Add anomalies (10 anomalous pixels)
    anomalies = np.random.multivariate_normal(mean + 5, np.eye(b)*0.1, size=10)
    for i in range(10):
        background[i*10, i*10] = anomalies[i]
    
    # Run Global RX
    global_rx_result = global_rx_detector(background)
    visualize_results(background, global_rx_result, "Global RX Result")
    
    # Run Local RX (slower but better for local anomalies)
    local_rx_result = local_rx_detector(background, window_size=15, inner_window=3)
    visualize_results(background, local_rx_result, "Local RX Result (15x15 window)")

# Example Usage
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="RX Anomaly Detection")

    parser.add_argument("input_image", help="Path to the input image")
    parser.add_argument("-method", choices=["global", "local"], default="global", help="Detector to use")
    parser.add_argument("-bb", type=int, default=30, help="Batch size for background sampling")
    parser.add_argument("-brx", type=int, default=0, help="Batch size for RX")
    parser.add_argument("-bm", choices=["none", "median", "mean", "pca", "morph"], default="median", help="Background method to use")
    parser.add_argument("-o", help="Prefix for output files")
    parser.add_argument("-enhance", type=str, default="none", help="Enhance contrast before detection")
    parser.add_argument("-local_window", type=int, default=15, help="Window size for local RX")
    parser.add_argument("-local_inner", type=int, default=3, help="Inner window size for local RX")
    args = parser.parse_args()
    
    if args.input_image is None:
        print("Please provide an input image.")
        sys.exit(1)

    if args.o is not None:
        output_prefix = args.o
    else:
        output_prefix= os.path.splitext(args.input_image)[0]

    # Load the hyperspectral image
    print("Loading image...")
    dataset = gdal.Open(args.input_image)
    if dataset is None:
        sys.exit(1)
    image = dataset.ReadAsArray().transpose(1, 2, 0).astype(np.float32)  # Convert to H x W x B
    dataset = None  # Close the dataset

    if args.enhance in ["median", "mean", "tophat"]:
        print(f"Enhancing by {args.enhance}...")
        if args.enhance == "median":
            for i in range(image.shape[2]):
                image[:, :, i] = image[:,:,i] - cv.medianBlur(image[:, :, i], args.local_inner)
        elif args.enhance == "mean":
            for i in range(image.shape[2]):
                image[:, :, i] = image[:,:,i] - cv.GaussianBlur(image[:, :, i], args.local_inner, 0)
        elif args.enhance == "tophat":
            kernel = cv.getStructuringElement(cv.MORPH_RECT, (args.local_window, args.local_window))
            for i in range(image.shape[2]):
                image[:, :, i] = cv.morphologyEx(image[:, :, i], cv.MORPH_TOPHAT, kernel)

        driver = gdal.GetDriverByName('GTiff')
        out_ds = driver.Create( output_prefix + '_enhance.tif', image.shape[1], image.shape[0], image.shape[2], gdal.GDT_Float32)
        for i in range(image.shape[2]):
            out_ds.GetRasterBand(i + 1).WriteArray(image[:, :, i])
        out_ds.FlushCache()
        out_ds = None

    # Run detection
    if args.bm in ['median', 'mean', 'pca', 'morph']:
        print(f"Removing background by {args.bm}...")
        foreground = remove_background(image, args.bb if args.bb < image.shape[2] else image.shape[2], args.bm)
        # Save the foreground image
        driver = gdal.GetDriverByName('GTiff')
        out_ds = driver.Create( output_prefix + '_foreground.tif', image.shape[1], image.shape[0], image.shape[2], gdal.GDT_Float32)
        for i in range(image.shape[2]):
            out_ds.GetRasterBand(i + 1).WriteArray(foreground[:, :, i])
        out_ds.FlushCache()
        out_ds = None
        image = foreground

    batch_rx = args.brx
    if batch_rx < 1:
        batch_rx = args.bb if args.bb < image.shape[2] else image.shape[2]

    detection_result = np.zeros((image.shape[0], image.shape[1], np.ceil(image.shape[2]/batch_rx).astype(int)))

    print(f"Batch RX with size {batch_rx}...")
    # Batch RX detection
    for i in range(detection_result.shape[2]):
        print(f"\t[{i+1}/{detection_result.shape[2]}] Detecting by {args.method}...")
        start = i*batch_rx
        end = (start + batch_rx) if i < detection_result.shape[2] - 1 else image.shape[2]
        batch_image = image[:, :, start:end]
        if args.method == 'global':
            detection_result[:,:, i] = global_rx_detector(batch_image, background_mask=None)
        elif args.method == 'local':
            detection_result[:,:, i] = local_rx_detector(batch_image, args.local_window, args.local_inner)
        detection_result[:,:, i] /= (end - start)

    print(f"Saving to {output_prefix+ '_rx_result.tif'}...")
    # Save the detection result
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_prefix + '_rx_result.tif', 
                           detection_result.shape[1], detection_result.shape[0], detection_result.shape[2], gdal.GDT_Float32)
    for i in range(detection_result.shape[2]):
        out_ds.GetRasterBand(i + 1).WriteArray(detection_result[:, :, i])
    out_ds.FlushCache()
    out_ds = None
    
    # Visualize results
    visualize_results(image, np.sum(detection_result, axis=2),"RX Anomaly Detection Result", output = output_prefix + '_rx_result.png')