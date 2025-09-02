import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def calculate_curvature_finite_difference(points):
    """
    Calculate curvature using finite differences
    Formula: κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
    """
    x = points[:, 0]
    y = points[:, 1]
    
    # First derivatives
    dx = np.gradient(x)
    dy = np.gradient(y)
    
    # Second derivatives
    d2x = np.gradient(dx)
    d2y = np.gradient(dy)
    
    # Calculate curvature
    numerator = np.abs(dx * d2y - dy * d2x)
    denominator = (dx**2 + dy**2)**1.5
    
    # Avoid division by zero
    curvature = np.zeros_like(numerator)
    mask = denominator > 1e-10
    curvature[mask] = numerator[mask] / denominator[mask]
    
    return curvature

def calculate_curvature_spline(points, smoothing=0):
    """
    Calculate curvature using spline interpolation
    """
    x = points[:, 0]
    y = points[:, 1]
    
    # Create parameter (arc length)
    t = np.zeros(len(x))
    for i in range(1, len(x)):
        t[i] = t[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
    
    # Fit splines
    tck_x = interpolate.splrep(t, x, s=smoothing)
    tck_y = interpolate.splrep(t, y, s=smoothing)
    
    # Evaluate derivatives
    t_eval = np.linspace(t[0], t[-1], len(t))
    
    # First derivatives
    dx_dt = interpolate.splev(t_eval, tck_x, der=1)
    dy_dt = interpolate.splev(t_eval, tck_y, der=1)
    
    # Second derivatives
    d2x_dt2 = interpolate.splev(t_eval, tck_x, der=2)
    d2y_dt2 = interpolate.splev(t_eval, tck_y, der=2)
    
    # Calculate curvature
    numerator = np.abs(dx_dt * d2y_dt2 - dy_dt * d2x_dt2)
    denominator = (dx_dt**2 + dy_dt**2)**1.5
    
    curvature = np.zeros_like(numerator)
    mask = denominator > 1e-10
    curvature[mask] = numerator[mask] / denominator[mask]
    
    return curvature

def calculate_curvature_moving_window(points, window_size=3):
    """
    Calculate curvature using a moving window approach
    """
    curvature = np.zeros(len(points))
    
    for i in range(window_size//2, len(points) - window_size//2):
        # Get points in window
        window_points = points[i-window_size//2:i+window_size//2+1]
        
        # Fit circle to points
        try:
            center, radius = fit_circle(window_points)
            if radius > 0:
                curvature[i] = 1.0 / radius
        except:
            curvature[i] = 0
    
    return curvature

def fit_circle(points):
    """
    Fit a circle to points using least squares
    """
    x = points[:, 0]
    y = points[:, 1]
    
    # Setup linear system
    A = np.column_stack([2*x, 2*y, np.ones_like(x)])
    b = x**2 + y**2
    
    # Solve for circle parameters
    c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    
    center_x, center_y = c[0], c[1]
    radius = np.sqrt(c[2] + center_x**2 + center_y**2)
    
    return (center_x, center_y), radius

def calculate_angular_acceleration(points, time=None):
    """
    Calculate angular acceleration from 2D points
    Returns: angular_velocity, angular_acceleration
    """
    x = points[:, 0]
    y = points[:, 1]
    
    # If time not provided, assume uniform sampling
    if time is None:
        time = np.arange(len(points))
    
    # Calculate angles (in radians)
    angles = np.arctan2(y, x)
    
    # Handle angle wrapping (unwrap phase)
    angles_unwrapped = np.unwrap(angles)
    
    # Calculate angular velocity (first derivative)
    angular_velocity = np.gradient(angles_unwrapped, time)
    
    # Calculate angular acceleration (second derivative)
    angular_acceleration = np.gradient(angular_velocity, time)
    
    return angular_velocity, angular_acceleration

def calculate_curvature(points, method='spline', **kwargs):
    if method == 'spline': 
        return calculate_curvature_spline(points, kwargs.get('smoothing'))
    elif method == 'moving_window':
        return calculate_curvature_moving_window(points, kwargs.get('window_size', 3))
    else:
        return calculate_curvature_finite_difference(points)