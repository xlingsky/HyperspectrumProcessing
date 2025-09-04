# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import os
import sys
import errno
import datetime
import warnings
import tempfile
import subprocess
import numpy as np
import rasterio
import random
import cv2


# silent rasterio NotGeoreferencedWarning
warnings.filterwarnings("ignore",
                        category=rasterio.errors.NotGeoreferencedWarning)

# add the bin folder to system path
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
bin_dir = os.path.join(parent_dir, 'build')
os.environ['PATH'] = bin_dir + os.pathsep + os.environ['PATH']

def print_elapsed_time(since_first_call=False):
    """
    Print the elapsed time since the last call or since the first call.

    Args:
        since_first_call:
    """
    t2 = datetime.datetime.now()
    if since_first_call:
        print("Total elapsed time:", t2 - print_elapsed_time.t0)
    else:
        try:
            print("Elapsed time:", t2 - print_elapsed_time.t1)
        except AttributeError:
            print("Elapsed time:", t2 - print_elapsed_time.t0)
    print_elapsed_time.t1 = t2

def run(cmd, env=os.environ, timeout=None, shell=False, verbose=False):
    """
    Runs a shell command, and print it before running.

    Arguments:
        cmd: list of a command and its arguments, or as a fallback,
            a string to be passed to a shell that will be split into a list.
        env (optional, default value is os.environ): dictionary containing the
            environment variables
        timeout (optional, int): time in seconds after which the function will
            raise an error if the command hasn't returned

        TODO: remove the temporary `shell` argument once all commands use shell=False
        shell (bool): run the command in a subshell. Defaults to False.

    Both stdout and stderr of the shell in which the command is run are those
    of the parent process.
    """
    if verbose:
        print("\nRUN: %s" % cmd)
    t = datetime.datetime.now()
    if not isinstance(cmd, list) and not shell:
        cmd = cmd.split()
    subprocess.run(cmd, shell=shell, stdout=sys.stdout, stderr=sys.stderr,
                   env=env, timeout=timeout, check=True)
    if verbose:
        print(datetime.datetime.now() - t)


def matrix_translation(x, y):
    t = np.eye(3)
    t[0, 2] = x
    t[1, 2] = y
    return t


def rio_read_as_array_with_nans(im):
    """
    Read an image replacing gdal nodata value with np.nan

    Args:
        im: path to the input image file

    Returns:
        array: raster as numpy array
    """
    with rasterio.open(im, 'r') as src:
        array = src.read()
        nodata_values = src.nodatavals

    for band, nodata in zip(array, nodata_values):
        if nodata is not None:
            band[band == nodata] = np.nan

    return array.squeeze()


def rasterio_write(path, array, profile={}, tags={}):
    """
    Write a numpy array in a tiff or png file with rasterio.

    Args:
        path (str): path to the output tiff/png file
        array (numpy array): 2D or 3D array containing the image to write.
        profile (dict): rasterio profile (ie dictionary of metadata)
        tags (dict): dictionary with additional geotiff tags
    """
    # determine the driver based on the file extension
    extension = os.path.splitext(path)[1].lower()
    if extension in ['.tif', '.tiff']:
        driver = 'GTiff'
    elif extension in ['.png']:
        driver = 'png'
    else:
        raise NotImplementedError('format {} not supported'.format(extension))

    # read image size and number of bands
    array = np.atleast_3d(array)
    height, width, nbands = array.shape

    # define image metadata dict
    profile.update(driver=driver, count=nbands, width=width, height=height,
                   dtype=array.dtype)

    # write to file
    with rasterio.Env():
        with rasterio.open(path, 'w', **profile) as dst:
            dst.write(np.transpose(array, (2, 0, 1)))
            dst.update_tags(**tags)

def bounding_box2D(pts):
    """
    bounding box for the points pts
    """
    dim = len(pts[0])  # should be 2
    bb_min = [min([t[i] for t in pts]) for i in range(dim)]
    bb_max = [max([t[i] for t in pts]) for i in range(dim)]
    return bb_min[0], bb_min[1], bb_max[0] - bb_min[0], bb_max[1] - bb_min[1]


def crop_array(img, x, y, w, h, fill_value=0):
    """
    Crop an image represented as an array.

    Args:
        img (array): 2D input image
        x, y (ints): coordinate of the top-left corner of the crop
        w, h (ints): width and height of the crop
        fill_value (img.dtype): constant value used for filling the crop
            outside the input image domain

    Returns:
        array with the cropped image
    """
    crop = fill_value * np.ones((h, w), dtype=img.dtype)

    y0 = max(y, 0)
    y1 = min(y + h, img.shape[0])
    x0 = max(x, 0)
    x1 = min(x + w, img.shape[1])

    if y0 < y1 and x0 < x1:  # the requested crop overlaps the image
        crop[y0 - y:y1 - y, x0 - x:x1 - x] = img[y0:y1, x0:x1]

    return crop


def print_elapsed_time(since_first_call=False):
    """
    Print the elapsed time since the last call or since the first call.

    Args:
        since_first_call:
    """
    t2 = datetime.datetime.now()
    if since_first_call:
        print("Total elapsed time:", t2 - print_elapsed_time.t0)
    else:
        try:
            print("Elapsed time:", t2 - print_elapsed_time.t1)
        except AttributeError:
            print("Elapsed time:", t2 - print_elapsed_time.t0)
    print_elapsed_time.t1 = t2
    print()


def linear_stretching_and_quantization_8bit(img, p=1):
    """
    Simple 8-bit quantization with linear stretching.

    Args:
        img (np.array): image to quantize
        p (float): percentage of the darkest and brightest pixels to saturate,
            from 0 to 100.

    Returns:
        numpy array with the quantized uint8 image
    """
    a, b = np.nanpercentile(img, (p, 100 - p))
    return np.round(255 * (np.clip(img, a, b) - a) / (b - a)).astype(np.uint8)

def scan_existing_files(directory, callback):
    for file in os.listdir(directory):
        callback(directory, file)

def concatenate_images(image_paths, directory = None):
    images = []
    for path in image_paths:
        if directory is not None:
            path = os.path.join(directory, path)
        with rasterio.open(path, 'r') as src:
            img = src.read()
            images.append(img)
    return np.concatenate(images, axis=0)

def generate_random_color(seed=None):
    """
    Generate a random color.

    :param seed: An optional seed for the random number generator.
    :return: A tuple (r, g, b) representing the color.
    """
    if seed is not None:
        random.seed(seed)

    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)

    return (r, g, b)

def draw( image, seeds, color, text = None, kwargs = dict()):
    radius = kwargs.get('radius', 2)
    thickness = kwargs.get('thickness', 1)
    font_scale = kwargs.get('fontScale', 0.3)
    for seed in seeds:
        cv2.circle(image, (int(seed[0]), int(seed[1])), radius, color, thickness )

    if text is not None:
        if isinstance(text, str):
            text = [text]*len(seeds)
        for seed, t in zip(seeds, text):
            cv2.putText(image, t, (int(seed[0]), int(seed[1])-radius), cv2.FONT_HERSHEY_SIMPLEX, font_scale, color)

    if kwargs.get('showCoordinates', False) :
        for seed in seeds:
            cv2.putText(image, f'({seed[0]:.1f},{seed[1]:.1f})', (int(seed[0]), int(seed[1])+radius), cv2.FONT_HERSHEY_SIMPLEX, font_scale, color)

    if kwargs.get('linked', False):
        for i in range(len(seeds)-1):
            s1 = seeds[i]
            s2 = seeds[i+1]
            cv2.line(image, (int(s1[0]), int(s1[1])), (int(s2[0]), int(s2[1])), color , thickness )

    return image

def normalize_to_8bit_rgb(image):
    if image.dtype != np.float32:
        image = image.astype(np.float32)

    normalized = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX)
    u8 = normalized.astype(np.uint8)

    return cv2.cvtColor(u8, cv2.COLOR_GRAY2BGR)

def rasterio_read_as_rgb24(path, window=None):
    with rasterio.open(path, 'r') as src:
        img = src.read(window=rasterio.windows.Window(*window) if window is not None else None)[0,:,:]

        return normalize_to_8bit_rgb(img)

def get_image_shape(path):
    with rasterio.open(path, 'r') as src:
        return src.width, src.height, src.count