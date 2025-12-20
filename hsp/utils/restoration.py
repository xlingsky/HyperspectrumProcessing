from skimage import restoration, filters
import numpy as np

CONFIG_GAUSSIAN_BLUR = {
    'apply': True,
    'sigma': None,  # If None, will be computed from radius and truncate
    'radius': [3, 5],
    'truncate': 4.0
}
CONFIG_NON_LOCAL_MEANS = {
    'apply': True,
    'h': 0.115,
    'patch_size': 5,
    'patch_distance': 7,
    'fast_mode': True,
}

CONFIG_FRANGI = {
    'apply': True,
    'sigmas': range(1, 3),
}

def denoise(image, config):
    for method, params in config.items():
        if not params['apply']:
            continue
        if method == 'gaussian_blur':
            sigma = params['sigma'] or np.array(params['radius']) / params['truncate']
            image = filters.gaussian(image, sigma=sigma)
        elif method == 'non_local_means':
            image = restoration.denoise_nl_means(image, h=params['h'],
                                                 patch_size=params['patch_size'],
                                                 patch_distance=params['patch_distance'],
                                                 fast_mode=params['fast_mode'])
        elif method == 'frangi':
            image = filters.frangi(image, sigmas=params['sigmas'])
        else:
            raise ValueError(f'Unknown denoising method: {method}')

    return image


