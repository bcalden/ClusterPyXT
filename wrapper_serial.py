import pix2pix
import cluster
import argparse
#import multiprocessing as mp
import numpy as np
import time
import pypeline_io as io

def get_arguments():
    help_str = """
    pypeline help string
    """
    prog = 'pypeline'

    parser = argparse.ArgumentParser(description=help_str, prog=prog)
    parser.add_argument("--cluster_config_file", "-c", dest="cluster_config",
                        action="store", default='',
                        help="Points to the cluster config file")
    parser.add_argument("--continue", dest='cont',
                        action='store_true', default=False,
                        help='Continue an unfinished run')
    parser.add_argument("--resolution", "-r", dest='resolution',
                        action='store', default=2, type=int,
                        help='Generate a low, medium, or high resolution temperature map. Low = 1, Med = 2, High = 3. '
                             'High resolution is a fit for every pixel, medium pixels are 3x3, low pixels are 5x5.')

    args = parser.parse_args()


    return args


if __name__ == '__main__':
    args = get_arguments()
    if args.cluster_config is not None:
        clstr = cluster.read_cluster_data(args.cluster_config)
        if args.cont:
            print('Continuing spectral fits')
            regions = clstr.unfinished_regions_to_fit(args.resolution)
            original = len(clstr.scale_map_regions_to_fit(args.resolution))
            print('{reg} of {orig} regions left to fit.'.format(
                reg=len(regions),
                orig=original)
            )
        else:
            clstr.initialize_best_fits_file()
            clstr.initialize_worst_fits_file()
            clstr.initialize_all_fits_file()
            regions = clstr.scale_map_regions_to_fit(args.resolution)

        #region_lists = np.array_split(regions, args.num_cpus)
        start_time = time.time()

        for region in regions:
            print("Running pix2pix on region {region}".format(region=region))
            pix2pix.pix_to_pix(clstr, region)

        print("Time elapsed: {time:0.2f} s.".format(time=time.time()-start_time))