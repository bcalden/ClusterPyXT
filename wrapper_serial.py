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

    #max_cpu = mp.cpu_count()

    # logger.debug("Getting commandline arguments.")
    parser = argparse.ArgumentParser(description=help_str, prog=prog)
    parser.add_argument("--cluster_config_file", "-c", dest="cluster_config",
                        action="store", default='/Volumes/blackbox/A115/A115_pypeline_config.ini',
                        help="Points to the cluster config file")
    # parser.add_argument("--num_cpus", "-n", dest="num_cpus",
    #                     action="store", default=max_cpu,
    #                     help="Number of CPUs to use in the fitting process. "
    #                     "Default is max number of cpus available ({} for "
    #                     "this system.)".format(max_cpu))
    args = parser.parse_args()

    #args.num_cpus == int(args.num_cpus)

    return args


# def run_in_parallel(cluster_obj, region_list, start_time):
#     regions_processed = []
#     num_regions = len(region_list)
#     for region in region_list:
#         pix2pix.pix_to_pix(cluster_obj, region)
#         regions_processed.append(region)
#
#         if len(regions_processed) % 10 == 0:
#             current_time = time.time()
#             elapsed = (current_time - start_time)# in seconds
#
#             num_regions_processed = len(regions_processed)
#             seconds_per_region = elapsed/num_regions_processed
#             name = mp.current_process().name
#
#             completion_seconds = ((num_regions-num_regions_processed)*seconds_per_region) # in seconds
#             completion_time = time.strftime("%H Hours %M minutes %S seconds.",
#                                            time.gmtime(completion_seconds))
#
#             log_string = "\n{name}: Elapsed time (s): {elapsed}\n"\
#                          "{name}: seconds_per_region: {seconds_per_region}\n"\
#                          "{num_processed} of {total_regions} processed\n" \
#                          "Estimated time for completion: {completion}.\n".format(
#                 name=name,
#                 elapsed=elapsed,
#                 seconds_per_region=seconds_per_region,
#                 num_processed=num_regions_processed,
#                 total_regions=num_regions,
#                 completion=completion_time
#             )
#             print(log_string)
#             log_file = "{super_comp_dir}/{p_name}_time.txt".format(
#                 super_comp_dir=cluster_obj.super_comp_dir,
#                 p_name=mp.current_process().name
#             )
#             io.write_contents_to_file(log_string, log_file, False)
#
#
#
#     mp.Queue().put('Regions procssed: {}'.format(len(regions_processed)))


if __name__ == '__main__':
    args = get_arguments()
    if args.cluster_config is not None:
        clstr = cluster.read_cluster_data(args.cluster_config)
        clstr.initialize_best_fits_file()
        clstr.initialize_worst_fits_file()
        regions = clstr.regions_to_fit
        #region_lists = np.array_split(regions, args.num_cpus)
        start_time = time.time()

        for region in regions:
            print("Running pix2pix on region {region}".format(region=region))
            pix2pix.pix_to_pix(clstr, region)

        print("Time elapsed: {time:0.2f} s.".format(time=time.time()-start_time))