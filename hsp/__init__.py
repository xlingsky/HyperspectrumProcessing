#!/usr/bin/env python

import sys
import os.path
import json
import datetime
import multiprocessing
import numpy as np
import rasterio
from sklearn.linear_model import RANSACRegressor

import fire
from matplotlib import pyplot as plt

from hsp.config import cfg
from hsp import common
# from hsp import parallel
from hsp import initialization

def make_path_relative_to_file(path, f):
    return os.path.join(os.path.abspath(os.path.dirname(f)), path)

def convert_meanstd_to_axb(meanfile, stdfile, afile, bfile):
    mean = np.loadtxt(meanfile)
    std = np.loadtxt(stdfile)
    a,b = common.meanstd_to_axb(mean, std)
    np.savetxt(afile, a, fmt='%.6e')
    np.savetxt(bfile, b, fmt='%.6e')

def compute_meanstd(path, im, tile = 20):
    with rasterio.open(im,'r') as src:
        data = src.read()[0]
        mean,std = common.compute_hstripe_meanstd(data, tile)
        np.savetxt(path+'mean.txt', mean, fmt='%.2f')
        np.savetxt(path+'std.txt', std, fmt='%.2f')

def test(meanfile = '/Users/xlingsky/Desktop/DATA/GF5B_mean.txt', im = '/Users/xlingsky/Desktop/DATA/GF5B_AHSI_SW_20220930_223_069_L00000194508_SW.tif'):
    mean = np.loadtxt(meanfile)
    x = mean[:,0].reshape(-1,1)
    y = mean[:,1310].reshape(-1,1)
    reg = RANSACRegressor().fit(x, y)
    plt.figure()
    plt.plot(x,y,"b+")
    plt.plot(x[reg.inlier_mask_,0], y[reg.inlier_mask_,0], "ro")
    dy = reg.predict(x)-y

def main(user_cfg):
    """
    Launch the hsp pipeline with the parameters given in a json file.

    Args:
        user_cfg: user config dictionary
        start_from: the step to start from (default: 0)
    """
    common.print_elapsed_time.t0 = datetime.datetime.now()
    initialization.build_cfg(user_cfg)
    initialization.make_dirs()

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes'] is not None:
        nb_workers = cfg['max_processes']

def read_config_file(config_file):
    """
    Read a json configuration file and interpret relative paths.

    If any input or output path is a relative path, it is interpreted as
    relative to the config_file location (and not relative to the current
    working directory). Absolute paths are left unchanged.
    """
    with open(config_file, 'r') as f:
        user_cfg = json.load(f)

    # output paths
    if not os.path.isabs(user_cfg['dstdir']):
        user_cfg['dstdir'] = make_path_relative_to_file(user_cfg['dstdir'],
                                                         config_file)

    # input paths
    for img in user_cfg['data']:
        for d in ['img', 'rpc']:
            if d in img and isinstance(img[d], str) and not os.path.isabs(img[d]):
                img[d] = make_path_relative_to_file(img[d], config_file)

    for task in user_cfg['tasklist']:
        for t in ['template']:
            if t in task and not os.path.isabs(task[t]):
                task[t] = make_path_relative_to_file(task[t], config_file)

    return user_cfg

if __name__ == "__main__":
    argv = sys.argv
    # compute_meanstd(argv[1], argv[2])
    test(argv[1]+'mean.txt')