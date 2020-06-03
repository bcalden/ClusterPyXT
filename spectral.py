import cluster
import argparse
import multiprocessing as mp
import numpy as np
import time
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
                        action="store_true", default=True,
                        help='Run in parallel (default False)')
    parser.add_argument("--num_cpus", "-n", dest="num_cpus",
                        action="store", default=max_cpu, type=int,
                        help="Number of CPUs to use in the fitting process. "
                        "Default is max number of cpus available ({} for "
                        "this system.)".format(max_cpu))
    parser.add_argument("--continue", dest='cont',
                        action='store_true', default=True,
                        help='Continue an unfinished run')
    parser.add_argument("--resolution", "-r", dest='resolution',
                        action='store', default=2, type=int,
                        help='Generate a low, medium, or high resolution temperature map. Low = 1, Med = 2, High = 3. '
                        'High resolution is a fit for every pixel, medium pixels are 3x3, low pixels are 5x5.')
    parser.add_argument("--dev-mtpc", dest='dev_mtpc', type=int, action="store",
                        default=0, help="Max Task Per Child sent to mp.Pool()")
    args = parser.parse_args()

    return args


def print_stage_tmap_prep(cluster: cluster.ClusterObj):
    prep_str = """Now ready for spectral fitting. 

    If offloaded, copy the spectral fits file, {spectral_fits},
    back to the local machine. 

    Next, run 

        python acb.py --temperature_map --resolution 2 --cluster_config_file /path/to/cluster/A115/A115_pypeline_config.ini

    This will create the temperature map and allow for the creation of the pressure maps.""".format(
        spectral_fits=cluster.spec_fits_file
    )

    print(prep_str)


def print_iteration_string(start_time, iteration, total):
    i = iteration
    num_region_lists = total
    try:
        elapsed = time.time() - start_time
        time_elapsed_str = time.strftime("%H hours %M minutes %S seconds",
                                     time.gmtime(elapsed))

        time_per_iteration = elapsed/iteration
        time_per_iteration_str = time.strftime("%H hours %M minutes %S seconds",
                                           time.gmtime(time_per_iteration))
        time_remaining = time_per_iteration * (total-iteration)
        time_remaining_str = time.strftime("%H hours %M minutes %S seconds",
                                       time.gmtime(time_remaining))
        io.print_red("\n*********************************\n")
        io.print_red("Iteration {} of {}".format(i + 1, num_region_lists))
        io.print_red("Time elapsed: {elapsed}".format(elapsed=time_elapsed_str))
        io.print_red("Approximate time per iteration: {}".format(time_per_iteration_str))
        io.print_red("Approximate time remaining: {}".format(time_remaining_str))
        io.print_red("\n*********************************\n")
    except ZeroDivisionError:
        io.print_red("\n*********************************\n")
        io.print_red("Iteration {} of {}".format(i + 1, num_region_lists))
        io.print_red("\n*********************************\n")


def calculate_spectral_fits(clstr: cluster.ClusterObj, num_cpus=mp.cpu_count()):
    args = get_arguments()
    if args.cont \
    and io.file_exists(clstr.spec_fits_file) \
    and io.file_exists(clstr.bad_fits_file) \
    and io.file_exists(clstr.all_fits_file):
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
    print("Regions to fit: {reg}".format(reg=len(regions)))

    start_time = time.time()
    num_regions = len(regions)
    if num_cpus:
        args.num_cpus = num_cpus

    if args.parallel and args.num_cpus > 1:
        num_region_lists = (num_regions//args.num_cpus)+1
        region_lists = np.array_split(regions, num_region_lists)
        print("Starting {} iterations with ~{} fits per iteration.".format(len(region_lists), len(region_lists[0])))
        for i, small_region_list in enumerate(region_lists):
            print_iteration_string(start_time, i, num_region_lists)
            processes = [mp.Process(target=clstr.fit_region_number, args=(region,)) for region in small_region_list]

            for process in processes:
                process.start()

            for process in processes:
                process.join()
    else:
        num_regions = len(regions)
        counter = 0
        for region in regions:
            print("Fitting region number {region}".format(region=region))
            clstr.fit_region_number(region)
            counter += 1
            if counter % 10 == 0 or counter == num_regions:
                print("{} of {} regions complete".format(counter, num_regions))
                time_elapsed = time.strftime("%H hours %M minutes %S seconds.",
                                                time.gmtime(time.time() - start_time))
                print("Time elapsed: {}.".format(time_elapsed))
    
    print_stage_tmap_prep(clstr)



if __name__ == '__main__':
    args = get_arguments()
    if args.cluster_config is not None:
        clstr = cluster.load_cluster(args.cluster_config)
        calculate_spectral_fits(clstr)
    else:
        print("Check help string. spectral.py --help")
            
     
