import os
import shutil
import argparse

import hsp


def main():
    """
    Command line parsing for hsp command line interface.
    """
    parser = argparse.ArgumentParser(description=('HSP: Hyper-spectral Pipeline'))
    parser.add_argument('config', metavar='config.json',
                        help=('path to a json file containing the paths to '
                              'input and output files and the task xml'
                              'specifies the algorithm parameters'))
    # parser.add_argument('--start_from', dest='start_from', type=int,
    #                     default=0, help="Restart the process from a given step in "
    #                                     "case of an interruption or to try different parameters.")
    args = parser.parse_args()

    user_cfg = hsp.read_config_file(args.config)

    hsp.main(user_cfg)

    # Backup input file for sanity check
    if not args.config.startswith(os.path.abspath(hsp.cfg['dstdir']+os.sep)):
        shutil.copy2(args.config,os.path.join(hsp.cfg['dstdir'],'config.json.orig'))
