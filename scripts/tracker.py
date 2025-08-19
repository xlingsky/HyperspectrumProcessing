import cv2 as cv
from osgeo import gdal
import numpy as np
import random

def read_seeds_file(seeds_file):
    """
    Read a seeds file and return a list of seed points.
    
    Args:
    seeds_file (str): Path to the seeds file.
    
    Returns:
    list: A list of tuples representing seed points (x, y).
    """
    all_seeds = list()
    with open(seeds_file, 'r') as f:
        line = f.readline()
        num = int(line.strip())
        for i in range(num):
            line = f.readline()
            if not line.strip():
                break  
            frame = int(line.strip())
            line = f.readline()
            sz = int(line.strip())
            seeds = []
            for j in range(sz):
                x = [int(float(x)) for x in f.readline().strip().split()]
                seeds.append(x[:2])
            all_seeds.append((frame, seeds))
    return all_seeds

def draw_trajectory(image, trajectory, text, color=255):
    """
    Draw a trajectory on the image.

    Args:
    image (numpy.ndarray): The image on which to draw the trajectory.
    trajectory (list): A list of tuples representing the trajectory points (x, y).

    Returns:
    numpy.ndarray: The image with the trajectory drawn.
    """
    for i in range(len(trajectory) - 1):
        cv.line(image, tuple(trajectory[i]), tuple(trajectory[i + 1]), color)
    cv.putText(image, text, tuple(trajectory[0]), cv.FONT_HERSHEY_SIMPLEX, 0.5, color)
    return image


def draw_seeds(image, seeds, color=255):
    """
    Draw seed points on the image.

    Args:
    image (numpy.ndarray): The image on which to draw the seeds.
    seeds (list): A list of tuples representing seed points (x, y).

    Returns:
    numpy.ndarray: The image with seed points drawn.
    """
    for seed in seeds:
        cv.rectangle(image,
                     (int(seed[0]), int(seed[1])),
                     (int(seed[2]), int(seed[3])),
                     color)
    return image

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

def normalize_to_8bit_rgb(image):
    if image.dtype != np.float32:
        image = image.astype(np.float32)

    normalized = cv.normalize(image, None, 0, 1, cv.NORM_MINMAX)
    u8 = (normalized*255).astype(np.uint8)

    return cv.cvtColor(u8, cv.COLOR_GRAY2BGR)

if __name__ == "__main__":
    import sys
    import os
    if len(sys.argv) < 3:
        print("Usage: python tracker.py <seeds_file> frames")
        sys.exit(1)

    seeds_file = sys.argv[1]
    trajectories = read_seeds_file(seeds_file)

    image_file = sys.argv[2]
    src = gdal.Open(image_file)
    if src is None:
        print(f"Could not open image file: {image_file}")
        sys.exit(1)
    cols, rows = src.RasterXSize, src.RasterYSize

    dst = np.zeros((rows, cols, 3), dtype=np.uint8)

    for frame, trajectory in trajectories:
        draw_trajectory(dst, trajectory, str(frame)+"+"+str(len(trajectory)), color=generate_random_color(frame))

    output_file = os.path.splitext(image_file)[0] + "_trajectories.tif"
    cv.imwrite(output_file, dst)

    bandseeds = {i:[] for i in range(src.RasterCount)}
    for id, (frame, trajectory) in enumerate(trajectories):
        color = generate_random_color(frame)
        for i, pt in enumerate(trajectory):
            bandseeds[frame+i].append([pt, color, id])

    output_file = os.path.splitext(image_file)[0] + "_trajectories.mp4"
    video = cv.VideoWriter(output_file, cv.VideoWriter_fourcc(*'mp4v'), 1, (src.RasterXSize, src.RasterYSize))
    imagedata = src.ReadAsArray()
    for i in range(src.RasterCount):
        rgb = normalize_to_8bit_rgb(imagedata[i, :, :])
        seeds = bandseeds[i]

        if len(seeds) > 0:
            for pt, color, id in seeds:
                cv.circle(rgb, tuple(pt), 3, color, -1)
                cv.putText(rgb, str(id), tuple(pt), cv.FONT_HERSHEY_SIMPLEX, 0.5, color)

        video.write(rgb)

    video.release()