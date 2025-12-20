import numpy as np
from scipy.optimize import curve_fit

def fit_quadratic_2d(grid):
    """
    Fit a 2D quadratic surface to the given grid of points.

    Args:
        grid (np.array): 2D array of shape (m, n) representing the surface values.

    Returns:
        Coefficients of the fitted quadratic surface in the form:
        f(x, y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
    """
    m, n = grid.shape
    X, Y = np.meshgrid(np.arange(n), np.arange(m))
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    Z_flat = grid.flatten()
    
    A = np.column_stack((X_flat**2, Y_flat**2, X_flat*Y_flat, X_flat, Y_flat, np.ones_like(X_flat)))
    
    coeffs, residual, _, _ = np.linalg.lstsq(A, Z_flat, rcond=None)
    
    return coeffs, residual  # a, b, c, d, e, f

def find_peak_from_quadratic(coeffs):
    """
    Find the peak (maximum) of the fitted quadratic surface.

    Args:
        coeffs (np.array): Coefficients of the quadratic surface.

    Returns:
        (x_peak, y_peak): Coordinates of the peak.
    """
    a, b, c, d, e, f = coeffs
    
    denom = 4*a*b - c**2
    if denom == 0:
        raise ValueError("The quadratic surface does not have a unique peak.")
    
    x_peak = (c*e - 2*b*d) / denom
    y_peak = (c*d - 2*a*e) / denom
    
    return (x_peak, y_peak)

def fit_gaussian_2d(grid):
    """
    Fit a 2D Gaussian function to the given grid of points.

    Args:
        grid (np.array): 2D array of shape (m, n) representing the surface values.

    Returns:
        Coefficients of the fitted Gaussian function including amplitude, center, rotation, standard deviation and offset.
    """
    def gaussian_2d(xy, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
        """2D Gaussian function"""
        x, y = xy
        x0 = float(x0)
        y0 = float(y0)
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    
        g = offset + amplitude * np.exp(-(a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2))
        return g.ravel()

    m, n = grid.shape
    X, Y = np.meshgrid(np.arange(n), np.arange(m))
    initial_guess = (np.max(grid)-np.min(grid), n//2, m//2, 1, 1, 0, np.min(grid))

    popt, pcov = curve_fit(gaussian_2d, (X.ravel(), Y.ravel()), grid.ravel(), p0=initial_guess, maxfev=5000)
    return popt, pcov  # amplitude, x0, y0, sigma_x, sigma_y, theta, offset
