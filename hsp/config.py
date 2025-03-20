# This module contains a dictionary, cfg, containing all the parameters of the
# hsp pipeline. This dictionary is updated at runtime with parameters defined
# by the user in the config.json file. All the optional parameters (that the
# user is not forced to define in config.json) must be defined here, otherwise
# they won't have a default value.

cfg = {}

# path to output directory
cfg['dstdir'] = "hsp_output"

# path to directory where (many) temporary files will be stored
cfg['tmpdir'] = "hsp_tmp"

# temporary files are erased when hsp terminates. Switch to False to keep them
cfg['clean_tmp'] = True

# remove all generated files except from ply point clouds and tif raster dsm
cfg['clean_intermediate'] = False

# hsp processes the images tile by tile. The tiles are squares cropped from the
# reference image. The width and height of the tiles are given by this param, in pixels.
cfg['tile_size'] = 800

# margins used to increase the footprint of the rectified tiles, to
# account for poor disparity estimation close to the borders
cfg['horizontal_margin'] = 50  # for regularity and occlusions
cfg['vertical_margin'] = 10  # for regularity

# max number of processes launched in parallel. None means the number of available cores
cfg['max_processes'] = None

# max number of OMP threads used by programs compiled with openMP
cfg['omp_num_threads'] = 1

# timeout in seconds, after which a function that runs on a single tile is not
# waited for
cfg['timeout'] = 600

# debug mode (more verbose logs and intermediate results saved)
cfg['debug'] = False
