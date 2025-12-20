import cv2
import numpy as np

from dataclasses import dataclass
from hsp.utils.fitting import fit_gaussian_2d, fit_quadratic_2d, find_peak_from_quadratic

def subpixel_parabolic(image, shift):
    """Estimate subpixel shift using parabolic fitting"""
    x, y = int(shift[0]), int(shift[1])
    assert x > 0 and x < image.shape[1]-1 and y > 0 and y < image.shape[0]-1

    dx = (image[y, x+1] - image[y, x-1]) / 2.0
    dy = (image[y+1, x] - image[y-1, x]) / 2.0
    dxx = image[y, x+1] - 2*image[y, x] + image[y, x-1]
    dyy = image[y+1, x] - 2*image[y, x] + image[y-1, x]

    if dxx == 0 or dyy == 0:
        return np.array([0,0])

    subpixel_x = -dx / dxx
    subpixel_y = -dy / dyy

    return np.array([subpixel_x, subpixel_y])

def subpixel_gaussian(image, shift):
    """Estimate subpixel shift using Gaussian fitting"""
    x, y = int(shift[0]), int(shift[1])
    assert x > 0 and x < image.shape[1]-1 and y > 0 and y < image.shape[0]-1

    log_x1 = np.log(image[y, x-1])
    log_x2 = np.log(image[y, x])
    log_x3 = np.log(image[y, x+1])

    log_y1 = np.log(image[y-1, x])
    log_y2 = np.log(image[y, x])
    log_y3 = np.log(image[y+1, x])

    denom_x = 2*log_x2 - log_x1 - log_x3
    denom_y = 2*log_y2 - log_y1 - log_y3

    if denom_x == 0 or denom_y == 0:
        return np.array([0,0])

    subpixel_x = 0.5 * (log_x1 - log_x3) / denom_x
    subpixel_y = 0.5 * (log_y1 - log_y3) / denom_y

    return np.array([subpixel_x, subpixel_y])

class ShiftEstimator:
    @dataclass
    class Image:
        id: int
        data: np.ndarray

    @dataclass
    class Shift:
        shift : np.ndarray
        confidence : float

    def __init__(self, search_radius = (10,10), method = cv2.TM_CCOEFF_NORMED, retry = 2, fitting = None):
        self.method = method
        self.fitting = fitting
        self.radius = search_radius
        self.retry = retry

    def begin(self, image):
        self.src_x = self.Image(0, image)
        self.src_y = self.src_x
        self.shifts = [self.Shift(np.array([0,0]), 1.0)]
        self.progress = 0
        self.last_shift = np.array([0,0])
        return self


    def _findTranslation(self, img1, img2, ccmethod, search_radius, retry, fitting = None ):

        if fitting == 'phase_correlation':
            shift, response = self.phase_correlation(img1, img2)
            return shift

        radius = search_radius
        margin = 1
        for _ in range(retry):
            res = cv2.matchTemplate(img1, img2[radius[1]:-radius[1], radius[0]:-radius[0]], ccmethod)
            min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
            if ccmethod in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
                loc = min_loc
            else:
                loc = max_loc

            if loc[0] > margin and loc[0] < (img1.shape[1] - margin) and loc[1] > margin and loc[1] < (img1.shape[0] - margin):
                if fitting:
                    if fitting == 'ecc':
                        # subpixel fitting
                        warpmat = np.array([[1, 0, loc[0]], [0, 1, loc[1]]], dtype=np.float32)
                        try:
                            (cc,warpmat) = cv2.findTransformECC(img2[radius[1]:-radius[1], radius[0]:-radius[0]], img1, warpMatrix=warpmat, motionType=cv2.MOTION_TRANSLATION, criteria=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 50, 1e-6), inputMask=None, gaussFiltSize=5)
                        except Exception as e:
                            warpmat = np.array([[1, 0, loc[0]], [0, 1, loc[1]]])
                        return warpmat[:,2] - np.array(radius)
                    elif fitting == 'local_quadratic':
                        subloc = subpixel_parabolic(res, loc)
                        return np.array(loc) + subloc - np.array(radius)
                    elif fitting == 'local_gaussian':
                        subloc = subpixel_gaussian(res, loc)
                        return np.array(loc) + subloc - np.array(radius)
                    elif fitting == 'quadratic':
                        try:
                            coeffs, residual = fit_quadratic_2d(res)
                            peak = find_peak_from_quadratic(coeffs)
                            return np.array(peak) - np.array(radius)
                        except:
                            pass
                    elif fitting == 'gaussian':
                        try:
                            coeffs, pcov = fit_gaussian_2d(res)
                            peak = (coeffs[1], coeffs[2])
                            return np.array(peak) - np.array(radius)
                        except:
                            pass
                return np.array(loc) - np.array(radius)

            radius = (radius[0]*2, radius[1]*2)

        raise ValueError("Failed to find translation")

    def phase_correlation(self, img1, img2):
        """Estimate shift using phase correlation"""
        # Use OpenCV's phaseCorrelate function
        shift, response = cv2.phaseCorrelate( np.float32(img2), np.float32(img1))
        return np.array(shift), response

    def estimate(self, image):
        self.progress += 1
        # if self.fitting == 'phase_correlation':
        #     shift, confidence = self.phase_correlation(self.src_x.data, image)
        #     self.shifts.append(self.Shift(shift + np.array([self.shifts[self.src_x.id].shift[0], self.shifts[self.src_y.id].shift[1]]), confidence))
        #     self.src_x = self.Image(self.progress, image)
        #     self.src_y = self.src_x
        #     return self.shifts[-1].shift, self.shifts[-1].confidence
        try:
            if self.src_x.id == self.src_y.id:
                shift = self._findTranslation(self.src_x.data, image, self.method, self.radius, self.retry, self.fitting)
            else:
                sx = self._findTranslation(self.src_x.data, image, self.method, self.radius, self.retry, self.fitting)
                sy = self._findTranslation(self.src_y.data, image, self.method, self.radius, self.retry, self.fitting)
                shift = np.array([sx[0], sy[1]])

            self.shifts.append(self.Shift(shift+np.array([self.shifts[self.src_x.id].shift[0], self.shifts[self.src_y.id].shift[1]]), 0.5))
            if shift[0] < 1:
                pass
            else:
                if shift[0] >= 1 and shift[0] < 3:
                    for i, dx in enumerate(np.linspace(0, shift[0], num=self.progress-self.src_x.id, endpoint=False)[1:]):
                        self.shifts[self.src_x.id + 1 + i].shift[0] = self.shifts[self.src_x.id].shift[0] + dx
                self.src_x = self.Image(self.progress, image)

            if shift[1] < 1:
                pass
            else:
                if shift[1] >= 1 and shift[1] < 3:
                    for i, dy in enumerate(np.linspace(0, shift[1], num=self.progress-self.src_y.id, endpoint=False)[1:]):
                        self.shifts[self.src_y.id + 1 + i].shift[1] = self.shifts[self.src_y.id].shift[1] + dy
                self.src_y = self.Image(self.progress, image)
        except Exception as e:
            self.src_x = self.Image( self.progress, image)
            self.src_y = self.src_x
            self.shifts.append(self.Shift(np.array([0,0]), 0.0))

        return self.shifts[-1].shift, self.shifts[-1].confidence

    def shift(self, i):
        return self.shifts[i].shift, self.shifts[i].confidence
        
    @property
    def num_shifts(self):
        return len(self.shifts)

    def save(self, filepath):
        """Save estimated shifts to a text file"""
        with open(filepath, 'w') as f:
            for shift in self.shifts:
                f.write(f"{shift.shift[0]:.2f},\t{shift.shift[1]:.2f},\t{shift.confidence:.2f}\n")

        return True

    def load(self, filepath):
        """Load estimated shifts from a text file"""
        with open(filepath, 'r') as f:
            lines = f.readlines()

        self.shifts = []
        for line in lines:
            items = line.strip().split(',')
            self.shifts.append(self.Shift(np.array([float(items[0]), float(items[1])]), float(items[2])))

        return True

        
import os
from scipy import ndimage

class Tester:
    def __init__(self, output_dir = 'test'):
        self.output_dir = output_dir
        # os.makedirs(self.output_dir, exist_ok=True)

    def _create_micro_pattern(self, size):
        """Create pattern with high spatial frequency for subpixel accuracy"""
        h, w = size
        
        # Multiple frequency patterns
        pattern = np.zeros((h, w), dtype=np.float32)
        
        # 1. High frequency sine pattern (for localization)
        x = np.linspace(0, 20*np.pi, w)
        y = np.linspace(0, 20*np.pi, h)
        X, Y = np.meshgrid(x, y)
        pattern += 0.3 * np.sin(X) * np.sin(Y)
        
        # 2. Medium frequency pattern
        pattern += 0.2 * np.sin(X/3) * np.cos(Y/3)
        
        # 3. Fine checkerboard
        checker_size = 4
        checker = ((np.indices((h, w)).sum(axis=0) // checker_size) % 2)
        pattern += 0.2 * checker.astype(np.float32)
        
        # 4. Random texture (high frequency noise)
        np.random.seed(42)
        texture = np.random.randn(h, w)
        texture = ndimage.gaussian_filter(texture, sigma=0.5)
        pattern += 0.1 * (texture - texture.min()) / (texture.max() - texture.min())
        
        # Normalize
        pattern = (pattern - pattern.min()) / (pattern.max() - pattern.min())
        
        return pattern
    
    def _fourier_shift(self, image, dx, dy):
        """Precise subpixel shift using Fourier transform"""
        # Fourier transform
        f = np.fft.fft2(image)
        
        # Create shift in frequency domain
        h, w = image.shape
        y_freq = np.fft.fftfreq(h).reshape(-1, 1)
        x_freq = np.fft.fftfreq(w).reshape(1, -1)
        
        # Phase shift
        phase_shift = np.exp(-2j * np.pi * (dy * y_freq + dx * x_freq))
        
        # Apply shift and inverse transform
        shifted_f = f * phase_shift
        shifted = np.fft.ifft2(shifted_f).real
        
        return shifted

    def create_subpixel_sequence(self, base_image=None, 
                                translations=None,
                                noise_levels=None,
                                image_size=(256, 256)):
        """
        Create sequence with known subpixel translations
        
        Args:
            translations: List of (dx, dy) subpixel translations
            noise_levels: List of noise levels for each frame
        """
        if base_image is None:
            base_image = self._create_micro_pattern(image_size)
        
        h, w = base_image.shape
        
        if translations is None:
            # Create a spiral of subpixel translations
            num_frames = 20
            translations = []
            for i in range(num_frames):
                angle = i * 0.3
                radius = 0.8
                dx = radius * np.cos(angle)  # Subpixel X
                dy = radius * np.sin(angle)  # Subpixel Y
                translations.append((dx, dy))
        
        if noise_levels is None:
            noise_levels = [0.01] * len(translations)
        
        images = [base_image]
        gt_translations = [(0.0, 0.0)]
        min_val, max_val = base_image.min(), base_image.max()
        
        for (dx, dy), noise_level in zip(translations, noise_levels):
            # Use Fourier shift theorem for accurate subpixel shifting
            shifted = self._fourier_shift(base_image, dx, dy)
            
            # Add noise
            if noise_level > 0:
                noise = np.random.normal(0, noise_level, shifted.shape)
                shifted = np.clip(shifted + noise, min_val, max_val)
            
            images.append(shifted.astype(np.float32))
            gt_translations.append((-dx, -dy))
        
        return images, gt_translations

    def test(self, images, gt_translations):
        """
        Test for pure translation estimation
        """
        num_frames = len(images)
        
        errors_x = []
        errors_y = []
        estimated_translations = []

        estimator = ShiftEstimator(search_radius=(10,10), method=cv2.TM_CCOEFF_NORMED, retry=3, fitting='ecc')
        estimator.begin(images[0])
        
        for i in range(1, num_frames):
            # Initialize warp matrix
            
            try:
                estimated_translations.append(estimator.estimate(images[i])[0])
                dx, dy = estimated_translations[-1]
                
                # Compute errors
                gt_dx, gt_dy = gt_translations[i]
                error_x = abs(dx - gt_dx)
                error_y = abs(dy - gt_dy)
                
                errors_x.append(error_x)
                errors_y.append(error_y)
                
                print(f"Frame {i:2d}: GT({gt_dx:6.3f}, {gt_dy:6.3f}) | "
                      f"Est({dx:6.3f}, {dy:6.3f}) | "
                      f"Err({error_x:6.3f}, {error_y:6.3f})")
                
            except cv2.error as e:
                print(f"Frame {i}: failed - {str(e)}")
                estimated_translations.append((0, 0))
                errors_x.append(float('inf'))
                errors_y.append(float('inf'))
        
        return estimated_translations, errors_x, errors_y

if __name__ == "__main__":
    tester = Tester()

    print("Scenario 1: very small shifts in x")
    small_shifts = [( 0.1*i, 1.3*i) for i in range(20)]

    baseimage = cv2.imread('/Users/xlingsky/Desktop/demo/FJ12/B4-000038.tif', cv2.IMREAD_UNCHANGED)

    images, gt_translations = tester.create_subpixel_sequence(base_image= baseimage, translations=small_shifts)
    
    # image0 = cv2.normalize(images[0], None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    # image1 = cv2.normalize(images[-1], None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    # cv2.imshow('First Frame', image0)
    # cv2.imshow('Last Frame', image1)
    # cv2.waitKey()
    # cv2.destroyAllWindows()
    
    est, err_x, err_y = tester.test(images, gt_translations)

    print(f'Max error in x: {max(err_x):.4f}, y: {max(err_y):.4f}\n')
    print(f'Mean error in x: {np.mean(err_x):.4f}, y: {np.mean(err_y):.4f}\n')

    