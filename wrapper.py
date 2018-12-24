import pix2pix
import cluster
import argparse
import multiprocessing as mp
import numpy as np
import time
import gc
import pypeline_io as io


def get_arguments():
    help_str = """
    pypeline help string
    """
    prog = 'pypeline'

    max_cpu = mp.cpu_count()

    # logger.debug("Getting commandline arguments.")
    parser = argparse.ArgumentParser(description=help_str, prog=prog)
    parser.add_argument("--cluster_config_file", "-c", dest="cluster_config",
                        action="store",default=None,
                        help="Points to the cluster config file")
    parser.add_argument("--parallel", "-p", dest="parallel",
                        action="store_true", default=False,
                        help='Run in parallel (default False)')
    parser.add_argument("--num_cpus", "-n", dest="num_cpus",
                        action="store", default=max_cpu, type=int,
                        help="Number of CPUs to use in the fitting process. "
                        "Default is max number of cpus available ({} for "
                        "this system.)".format(max_cpu))
    parser.add_argument("--continue", dest='cont',
                        action='store_true', default=False,
                        help='Continue an unfinished run')
    parser.add_argument("--resolution", "-r", dest='resolution',
                        action='store', default=2, type=int,
                        help='Generate a low, medium, or high resolution temperature map. Low = 1, Med = 2, High = 3. '
                        'High resolution is a fit for every pixel, medium pixels are 3x3, low pixels are 5x5.')
    args = parser.parse_args()

    return args


def run_in_parallel(cluster_obj: cluster.ClusterObj, region_list, start_time):
    regions_processed = []
    num_regions = len(region_list)
    time_of_last = 0
    for region in region_list:
        pix2pix.pix_to_pix(cluster_obj, region, mp.current_process().name)
        regions_processed.append(region)
        if len(regions_processed) % 10 == 0:
            current_time = time.time()
            time_since_last = current_time - time_of_last
            time_of_last = time.time()

            average_second_per_completion = time_since_last/10

            elapsed = current_time - start_time
            elapsed_str = time.strftime("%H Hours %M minutes %S seconds.",
                                            time.gmtime(elapsed))


            num_regions_processed = len(regions_processed)
            seconds_per_region = elapsed/num_regions_processed
            name = mp.current_process().name

            completion_seconds = ((num_regions-num_regions_processed)*seconds_per_region)  # in seconds
            completion_time = time.strftime("%H Hours %M minutes %S seconds.",
                                            time.gmtime(completion_seconds))

            log_string = "\n{color}*******************************************{reset}\n"\
                         "\n{name}: Elapsed time: {elapsed}\n"\
                         "{name}: Seconds per region (total average): {seconds_per_region}\n"\
                         "Seconds per region (last 10 regions fit): {recent_average}\n"\
                         "{num_processed} of {total_regions} processed\n" \
                         "Estimated time for completion: {completion}.\n" \
                         "\n{color}*******************************************{reset}\n".format(
                color=io.Colors.RED,
                reset=io.Colors.RESET,
                name=name,
                elapsed=elapsed_str,
                seconds_per_region=seconds_per_region,
                recent_average=average_second_per_completion,
                num_processed=num_regions_processed,
                total_regions=num_regions,
                completion=completion_time
            )
            print(log_string)
            log_file = "{super_comp_dir}/{p_name}_time.txt".format(
                super_comp_dir=cluster_obj.super_comp_dir,
                p_name=mp.current_process().name
            )
            io.write_contents_to_file(log_string, log_file, False)
            print("Manual garbage collection.")
            gc.collect()


    mp.Queue().put('Regions processed: {}'.format(len(regions_processed)))


if __name__ == '__main__':
    args = get_arguments()
    if args.cluster_config is not None:
        clstr = cluster.read_cluster_data(args.cluster_config)
        if args.cont:
            print("Continuing spectral fits")
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

        start_time = time.time()

        if args.parallel:
            region_lists = np.array_split(regions, args.num_cpus)

            pool = mp.Pool(processes=args.num_cpus)

            processes = [mp.Process(target=run_in_parallel, args=(clstr, x, start_time)) for x in region_lists]

            for process in processes:
                process.start()

            for process in processes:
                process.join()

            results = [mp.Queue().get() for p in processes]

            for result in results:
                print(result)
        else:
            num_regions = len(regions)
            counter = 0
            for region in regions:
                print("Running pix2pix on region {region}".format(region=region))
                pix2pix.pix_to_pix(clstr, region)
                counter += 1
                if counter % 10 == 0 or counter == num_regions:
                    print("{} of {} regions complete".format(counter, num_regions))
                    time_elapsed = time.strftime("%H hours %M minutes %S seconds.",
                                                 time.gmtime(time.time() - start_time))
                    print("Time elapsed: {}.".format(time_elapsed))