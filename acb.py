import cluster
import pypeline_io as io
import numpy as np
import time
import argparse
import data_operations as do

def get_arguments():
    help_str = """
    This part of the pypeline creates all of the adaptive circular binned (acb) files
    needed to do the spectral fitting. This fitting should likely be offloaded onto a
    high performance computer.
    
    Sample call (with the ciao environment running):
    python acb.py --cluster_config_file /data_dir/A115/A115_pypeline_config.ini --resolution 2
    python acb.py --cluster_config_file /data/dir/A115/A115_pypeline_config.ini --temperature_map
    """
    prog = 'python acb.py'

    # logger.debug("Getting commandline arguments.")
    parser = argparse.ArgumentParser(description=help_str, prog=prog)
    parser.add_argument("--cluster_config_file", "-c", dest="cluster_config",
                        action="store", default=None,
                        help="Path to the cluster configuration file")
    # parser.add_argument("--parallel", "-p", dest="parallel",
    #                     action="store_true", default=False,
    #                     help='Run in parallel (default False)')
    parser.add_argument("--temperature_map", "-t", dest='temperature_map',
                        action="store_true", default=False,
                        help="Create a temperature map after the spectral fitting process.")
    parser.add_argument("--resolution", "-r", dest='resolution',
                        action='store', default=2, type=int,
                        help='Generate a low, medium, or high resolution temperature map. Low = 1, Med = 2, High = 3. '
                        'High resolution is a fit for every pixel, medium pixels are 3x3, low pixels are 5x5.')
    parser.add_argument('--make_fitting_commands', dest='commands', action='store_true', default=False)
    parser.add_argument('--eff_times_to_fits', dest='eff_times_fits', action='store_true', default=False)
    parser.add_argument('--make_pressure_map', dest='pressure', action='store_true', default=False)
    parser.add_argument('--shock_finder', dest='shock', action='store_true', default=False)
    args = parser.parse_args()

    return args, parser


def create_circle_regions(cluster):
    from astropy.io import fits

    start_time = time.time()

    scale_map_fits = fits.open(cluster.scale_map_file)
    mask_fits = fits.open(cluster.combined_mask)

    region_map = cluster.scale_map_region_index #region_map_fits[0].data

    scale_map = scale_map_fits[0].data
    mask = mask_fits[0].data
    bounds = scale_map.shape

    xvals = np.arange(bounds[1])
    yvals = np.arange(bounds[0])

    for observation in cluster.observations:
        print("Making circular fitting regions for observation {}".format(observation.id))
        image_fits = fits.open(observation.acisI_comb_img)
        image_header = image_fits[0].header

        cdelt1p = image_header['CDELT1P']
        cdelt2p = image_header['CDELT2P']
        crval1p = image_header['CRVAL1P']
        crval2p = image_header['CRVAL2P']
        crpix1p = image_header['CRPIX1P']
        crpix2p = image_header['CRPIX2P']

        radii = mask * scale_map * cdelt1p

        newx = ((xvals + 1 - crpix1p) * cdelt1p) + crval1p
        newy = ((yvals + 1 - crpix2p) * cdelt2p) + crval2p

        xx, yy = np.meshgrid(newx, newy)

        non_zero_indeces = np.nonzero(radii)
        nz_rad = radii[non_zero_indeces]
        nz_x = xx[non_zero_indeces]
        nz_y = yy[non_zero_indeces]
        obs_regions = region_map[non_zero_indeces]

        region_array = np.array([(i, j, k, l) for i, j, k, l in zip(nz_x, nz_y, nz_rad, obs_regions)])

        regions = [["circle({x},{y},{rad})".format(x=x[0], y=x[1], rad=x[2]), int(x[3])] for x in region_array]

        observation.scale_map_region_list = regions

    end_time = time.time()

    print("Time elapsed making regions for fit: {:0.2f} (s)".format(end_time-start_time))


def create_region_index_map(cluster):
    from astropy.io import fits

    mask_fits = fits.open(cluster.combined_mask)
    mask = mask_fits[0].data

    sz = mask.shape
    nx = sz[0]
    ny = sz[1]

    indexmap = np.zeros(sz)

    position = 0
    region_string = []

    for ci in range(nx):
        for cj in range(ny):
            if mask[ci,cj] == 1:
                if ci % 3 == 0 and cj % 3 == 0:  # makes it a lower resolution image than it needs to be.
                    region_string.append(str(position))
                indexmap[ci,cj] = position
                position += 1
    region_string = '\n'.join(region_string)
    with open(cluster.region_list, 'w') as f:
        f.write(region_string)

    region_file = mask_fits
    region_file[0].data = indexmap

    region_file[0].writeto(cluster.region_to_index, overwrite=True)


def create_scale_map_region_index(cluster: cluster.ClusterObj):
    from astropy.io import fits
    scale_map = cluster.scale_map
    scale_map_regions = np.zeros(scale_map.shape)
    sx = scale_map.shape[0]
    sy = scale_map.shape[1]
    region_num = 1
    for x in range(sx):
        for y in range(sy):
            if scale_map[x,y] != 0:
                scale_map_regions[x,y] = region_num
                region_num += 1
            else:
                scale_map_regions[x,y] = np.nan

    fits.writeto(cluster.scale_map_region_file,  # filename
                 scale_map_regions,  # data to write
                 cluster.scale_map_header,  # header so coordinate information is written
                 overwrite=True)  # self explanatory


def _update_completed_things(current, max_num, thing):
    io.clear_line()
    io.write("{current} out of {max} {thing} complete.                                                         ".format(
        current=current,
        max=max_num,
        thing=thing
    ))
    io.flush()

def _source_free_region(counter, current, max_num):
    io.clear_line()
    io.write("Encountered a source-free region -- recalculating...{counter} - {current}/{max} complete".format(
        counter=counter,
        current=current,
        max=max_num
    ))
    io.flush()

def _update_effective_exposure_time(current_region, number_regions, time_elapsed):
    io.clear_line()
    io.write("{current_region} of {num_regions} complete. Time elapsed: {time}".format(
        current_region=current_region,
        num_regions=number_regions,
        time=time_elapsed
    ))
    io.flush()

def create_scale_map(cluster):
    from astropy.io import fits

    target_sn = cluster.target_sn

    mask_fits = fits.open(cluster.combined_mask)

    cts_image = np.zeros(mask_fits[0].data.shape)
    back_rescale = np.zeros(mask_fits[0].data.shape)

    for obs in cluster.observations:
        obs_cts_image_fits = fits.open(obs.acisI_comb_img)
        obs_back_image_fits = fits.open(obs.backI_comb_img)

        cts_image += obs_cts_image_fits[0].data
        t_obs = obs_cts_image_fits[0].header['EXPOSURE']

        # back_rescale += obs_back_image_fits[0].data
        t_back = obs_back_image_fits[0].header['EXPOSURE']

        back_rescale += (t_obs/t_back)*obs_back_image_fits[0].data

    # sigma_back = np.sqrt(back_rescale)
    # sigma_cts = np.sqrt(cts_image)

    signal = cts_image - back_rescale

    signal[np.where(signal < 0)] = 0

    # noise = np.sqrt(sigma_cts**2 + sigma_back**2)
    #
    # scale_file = cluster.scale_map_file
    # sn_file = cluster.sn_map

    sz = signal.shape
    nx = sz[0]
    ny = sz[1]

    # radius_map = np.zeros(sz)

    scale_map = np.zeros(sz)
    sn_map = np.zeros(sz)

    pix_x = np.zeros(sz)
    pix_y = np.zeros(sz)

    for j in range(ny):
        for i in range(nx):
            pix_x[i, j] = float(i)
            pix_y[i, j] = float(j)

    bpix_x = 0
    bpix_y = 0
    epix_x = nx
    epix_y = ny

    max_radius = 100
    r = max_radius + 1  # +1 may be a remnant of IDL indexing
    #dr = 24.0
    min_dr = 0.125

    det = mask_fits[0].data

    num_pix = nx*ny

    ci=0
    counter=0
    start_time = time.time()

    # output_queue = mp.Queue()


    for cj in range(bpix_y,epix_y):
        #print("{} out of {} pixels complete.".format(cj*ci, num_pix))
        _update_completed_things(cj*ci, num_pix, "pixels")
        for ci in range(bpix_x, epix_x):
            if det[ci,cj] == 1:
                delta_x = ci-pix_x
                delta_y = cj-pix_y

                # an array where each value is the distance away from ci,cj
                radius = np.sqrt(delta_x**2 + delta_y**2)

                dr = 24.0
                hilo = 0
                niter = 0

                while (dr > min_dr) and (niter < 100):
                    indeces = np.where(radius <= r)
                    cts_map_total = np.sum(cts_image[indeces])

                    if cts_map_total == 0:
                        counter+=1
                        _source_free_region(counter, ci*cj, num_pix)
                        sn_val = 0
                        hilo = -1
                    else:
                        backmap_tot = np.sum(back_rescale[indeces])
                        signal_total = cts_map_total - backmap_tot
                        noise_total = np.sqrt(cts_map_total+backmap_tot)
                        sn_val = signal_total/noise_total
                    if sn_val < target_sn:
                        if r > max_radius:
                            r = max_radius + 1

                            niter = 110
                            # exit by setting niter=110.
                            # (niter=100 means niter hit max niter.
                            # niter=110 means radius hit max radius)

                            sn_map[ci,cj] = 0
                            scale_map[ci,cj] = 0
                        else:
                            if hilo == 1:
                                dr *= 0.5
                            r += dr
                            hilo = -1
                    else:
                        snmapval = sn_map[ci,cj]
                        if (sn_val < snmapval) or (snmapval == 0.0):
                            sn_map[ci,cj] = sn_val
                            scale_map[ci,cj] = r
                        if hilo == -1:
                            dr *= 0.5
                        r -= dr
                        hilo = 1

                    niter += 1
                # end while

    end_time = time.time()
    print("Time elapsed {:0.2f} seconds.".format(end_time - start_time))

    io.make_directory(cluster.acb_dir)

    obs_cts_image_fits[0].data = scale_map
    obs_cts_image_fits.writeto(cluster.scale_map_file, overwrite=True)

    obs_cts_image_fits[0].data = sn_map
    obs_cts_image_fits.writeto(cluster.sn_map, overwrite=True)


def prepare_efftime_circle(cluster):
    try:
        from ciao_contrib import runtool as rt
    except ImportError:
        print("Failed to import CIAO python scripts. ")
        raise

    for observation in cluster.observations:
        io.delete_if_exists(observation.effbtime)
        io.delete_if_exists(observation.effdtime)

        if not io.file_exists(observation.acisI_nosrc_combined_mask):
            print("Removing point sources from the observations combined mask file.")
            print("dmcopy infile='{}[exclude sky=region({})]' outfile={} clobber=True".format(
                observation.acisI_combined_mask,
                cluster.sources_file,
                observation.acisI_nosrc_combined_mask
            ))
            rt.dmcopy.punlearn()
            rt.dmcopy(
                infile="{fits_file}[exclude sky=region({source_file})]".format(
                    fits_file=observation.acisI_combined_mask,
                    source_file=cluster.sources_file
                ),
                outfile=observation.acisI_nosrc_combined_mask,
                clobber=True
            )
        else:
            print("{acis} already exists.".format(
                acis=observation.acisI_nosrc_combined_mask
            ))

        if not io.file_exists(observation.acisI_high_energy_combined_image):
            print("Creating high band (9.5-12 keV) source image cropped to combined region.")
            rt.dmcopy.punlearn()
            rt.dmcopy(
                infile="{fits_file}[sky=region({crop_file})]".format(
                    fits_file=observation.clean,
                    crop_file=cluster.master_crop_file
                ),
                outfile=observation.acisI_high_energy_temp_image,
                clobber=True
            )

            rt.dmcopy.punlearn()
            rt.dmcopy(
                infile="{fits_file}[EVENTS][bin sky=4][energy=9500:12000]".format(
                    fits_file=observation.acisI_high_energy_temp_image
                ),
                outfile=observation.acisI_high_energy_combined_image,
                option="image",
                clobber=True
            )
        else:
            print("{fits_file} already exists.".format(
                fits_file=observation.acisI_high_energy_combined_image
            ))

        io.delete_if_exists(observation.acisI_high_energy_temp_image)

        if not io.file_exists(observation.backI_high_energy_combined_image):
            print("Creating high band (9.5-12 keV) background image cropped to combined region.")
            rt.dmcopy.punlearn()
            rt.dmcopy(
                infile="{fits_file}[sky=region({crop_file})]".format(
                    fits_file=observation.back,
                    crop_file=cluster.master_crop_file
                ),
                outfile=observation.backI_high_energy_temp_image,
                clobber=True
            )

            rt.dmcopy.punlearn()
            rt.dmcopy(
                infile="{fits_file}[EVENTS][bin sky=4][energy=9500:12000]".format(
                    fits_file=observation.backI_high_energy_temp_image
                ),
                outfile=observation.backI_high_energy_combined_image,
                option="image",
                clobber=True
            )
        else:
            print("{fits_file} already exists.".format(
                fits_file=observation.backI_high_energy_combined_image
            ))

        io.delete_if_exists(observation.backI_high_energy_temp_image)


def calculate_effective_times(cluster: cluster.ClusterObj):
    from astropy.io import fits

    start_time = time.time()
    scale_map = cluster.scale_map
    number_of_regions = cluster.number_of_regions

    nx = scale_map.shape[0]
    ny = scale_map.shape[1]

    effective_data_times = np.zeros(scale_map.shape)
    effective_background_times = np.zeros(scale_map.shape)

    for observation in cluster.observations:
        print("Starting observation {obs}".format(obs=observation.id))
        high_energy_data = fits.open(observation.acisI_high_energy_combined_image)[0].data
        background = fits.open(observation.backI_high_energy_combined_image)[0].data

        sum_acis_high_energy = np.sum(high_energy_data)  # get the total counts in the high energy image
        sum_back_high_energy = np.sum(background)

        bg_to_data_ratio = sum_back_high_energy / sum_acis_high_energy

        source_subtracted_data = fits.open(observation.acisI_nosrc_combined_mask)[0].data
        exposure_time = fits.open(observation.acisI_high_energy_combined_image)[0].header['EXPOSURE']

        YY, XX = np.meshgrid(np.arange(ny), np.arange(nx))

        counter = 0
        print("Starting effective exposure time calculations...")
        for x in range(nx):
            for y in range(ny):
                if scale_map[x,y] >= 1:
                    radius = np.sqrt((x - XX)**2 + (y - YY)**2)
                    region = np.where(radius <= scale_map[x, y])

                    source_subtracted_area = np.sum(source_subtracted_data[region])
                    total_area = source_subtracted_data[region].size
                    fractional_area = source_subtracted_area / total_area

                    fractional_exposure_time = fractional_area * exposure_time

                    effective_data_times[x, y] = fractional_exposure_time
                    effective_background_times[x, y] = fractional_exposure_time * bg_to_data_ratio

                    counter += 1
                    if counter % 1000 == 0 or counter == number_of_regions or counter == 1:

                        time_elapsed = time.strftime("%H hours %M minutes %S seconds.",
                                                    time.gmtime(time.time()-start_time))
                        _update_effective_exposure_time(current_region=counter,
                                                        number_regions=number_of_regions,
                                                        time_elapsed=time_elapsed
                                                        )

        observation.effective_data_time = effective_data_times
        observation.effective_background_time = effective_background_times



def prepare_for_spec(cluster_obj: cluster.ClusterObj):
    try:
        import ciao
    except ImportError:
        print("Must be running CIAO before running prepare_for_spec.")
        raise
    io.make_directory(cluster_obj.super_comp_dir)
    #io.make_directory(cluster_obj.sherpa_save_dir)
    cluster_obj.initialize_best_fits_file()
    print("Preparing files for spectral analysis and copying to {super_comp_dir} for offloading computation.".format(
        super_comp_dir=cluster_obj.super_comp_dir
    ))
    io.copy(cluster_obj.configuration_filename, cluster_obj.super_comp_cluster_config)
    for observation in cluster_obj.observations:
        print("Copying files for {obsid}".format(obsid=observation.id))
        io.copy(observation.clean, cluster_obj.acisI_clean_obs(observation.id))
        io.copy(observation.back, cluster_obj.backI_clean_obs(observation.id))
        io.copy(observation.aux_response_file, observation.arf_sc)
        io.copy(observation.redistribution_matrix_file, observation.rmf_sc)
        io.copy(observation.acis_mask, observation.acis_mask_sc)

        exposure = ciao.get_exposure(observation.clean)

        io.write_contents_to_file(exposure, observation.exposure_time_file, binary=False)



def make_commands_lis(cluster: cluster.ClusterObj, resolution):
    from astropy.io import fits

    print("Creating {}".format(cluster.command_lis))

    offset = [None, 5, 3, 1][resolution]


    start_time = time.time()

    scalemap = cluster.scale_map
    # mask = fits.open(cluster.combined_mask)[0].data
    region_list = cluster.scale_map_regions_to_fit(resolution)
    #sz = scalemap.shape
    #sx = sz[0]
    #sy = sz[1]

    command_string = []

    pypeline_dir = io.get_user_input("Enter the directory containing the pix2pix.py portion of the pypeline on the remote machine: ")
    pix2pix_path = "{pypeline_dir}/pix2pix.py".format(pypeline_dir=pypeline_dir)

    data_dir = io.get_user_input("Enter the directory containing the cluster data on the remote machine:\n"
                                 "For example: /home/user/data/clustername/\n")

    for region in region_list:
        new_command = "python {pix2pix} {cluster_config} {region}".format(
            pix2pix=pix2pix_path,
            cluster_config="{data_dir}/{name}_pypeline_config.ini".format(
                data_dir=data_dir,
                name=cluster.name
            ),
            region=region
        )
        command_string.append(new_command)

    # for x in range(sx):
    #     for y in range(sy):
    #         if mask[i,j] == 1:
    #             image_radius = scalemap[i,j]
    #             if 0 < image_radius <= 100:
    #                 if i % offset == 0 and j % offset == 0:
    #                     new_command = "python {pix2pix} {cluster_config} {region}".format(
    #                         pix2pix=pix2pix_path,
    #                         cluster_config="{data_dir}/{name}_pypeline_config.ini".format(
    #                             data_dir=data_dir,
    #                             name=cluster.name
    #                         ),
    #                         region=position
    #                     )
    #                     command_string.append(new_command)
    #                     region_list.append("{}".format(position))
    #                 #end if
    #             # end if
    #             position += 1
    command_lis = "\n".join(command_string)
    region_string = '\n'.join([str(x) for x in region_list])
    io.write_contents_to_file(command_lis, cluster.command_lis, binary=False)
    io.write_contents_to_file(region_string, cluster.filtered_region_list, binary=False)
    end_time = time.time()
    print("Time elapsed: {time:0.2f} sec".format(time=(end_time-start_time)))


# def get_pixel_coordinates(cluster: cluster.ClusterObj):
#     region_map = fits.open(cluster.region_to_index)[0].data
#     coordinates = np.array(np.nonzero(region_map)).T
#     return coordinates


# def get_coords_for_region(cluster: cluster.ClusterObj, region_num: int, scale_map):
#     region_index_map = fits.open(cluster.scale_map_region_index)
#     sz = region_index_map.shape
#     for x in range(sz[0]):
#         for y in range(sz[1]):
#             if region_index_map[x,y] == region_num:
#                 return x,y
#     print("Fail!")
#     return 0,0

def make_temperature_map(cluster: cluster.ClusterObj, resolution, average=False):

    from astropy.io import fits

    #coordinates = get_pixel_coordinates(cluster)
    # indeces of this array are the region number minus 1
    # that is, region number 1 is coordinate array index 0
    # region 100 = coordinates[99]

    # high_res_offset = 0
    # med_res_offset = 1
    # low_res_offset = 2

    io.make_directory(cluster.output_dir)

    offset = [None, 2, 1, 0][resolution]

    mask_fits = fits.open(cluster.combined_mask)
    mask = mask_fits[0].data

    scale_map_regions = cluster.scale_map_region_index
    temps_with_errors = cluster.average_temperature_fits if average else cluster.temperature_fits

    temperature_map = np.zeros(mask.shape)
    temperature_error_map = np.zeros(mask.shape)
    temperature_fractional_error_map = np.zeros(mask.shape)

    regions = temps_with_errors['region']
    temperatures = temps_with_errors['temperature']
    temp_error_plus = temps_with_errors['temp_err_plus']
    temp_error_minus = temps_with_errors['temp_err_minus']
    for i, region in enumerate(regions):
        if i % 1000 == 0:
            _update_completed_things(i, len(regions), "regions")
        coordinates = cluster.coordinates_for_scale_map_region(region, scale_map_regions)
        #x,y = get_coords_for_region(cluster, regions[i])#coordinates[i-1]  # region 1 = index 0, reg 2 = index 1...
        x = int(coordinates[0])
        y = int(coordinates[1])
        low_x = x - offset
        high_x = x + offset + 1
        low_y = y - offset
        high_y = y + offset + 1

        temperature_map[low_x:high_x, low_y:high_y] = temperatures[i]
        temperature_error_map[low_x:high_x, low_y:high_y] = (np.abs(temp_error_plus[i] -
                                                             temp_error_minus[i]))/2
        temperature_fractional_error_map[low_x:high_x, low_y:high_y] = \
            (temperature_error_map[x,y]/temperature_map[x,y])*100

    _update_completed_things(i, len(regions), "regions")

    header = mask_fits[0].header
    # This header contains all coordinate information needed
    fits.writeto(cluster.temperature_map_filename,
                 temperature_map,
                 header,
                 overwrite=True)

    fits.writeto(cluster.temperature_error_map_filename,
                 temperature_error_map,
                 header,
                 overwrite=True)

    fits.writeto(cluster.temperature_fractional_error_map_filename,
                 temperature_fractional_error_map,
                 header,
                 overwrite=True)


def fitting_preparation(clstr, args=None):

    if args is None:
        resolution = 2
    else:
        resolution = args.resolution

    print("Creating the scale map.")
    create_scale_map(clstr)

    print("Creating the region index map.")
    create_scale_map_region_index(clstr)

    print("Preparing the high-energy images and backgrounds.")
    prepare_efftime_circle(clstr)

    print("Calculating effective times.")
    calculate_effective_times(clstr)

    print("Creating circular fitting regions.")
    create_circle_regions(clstr)

    print("Preparing for the spectral fits.")
    prepare_for_spec(clstr)

    print("Making the command list for use with mpiexec.")
    make_commands_lis(clstr, resolution)  # 3 for high_res, 2 for medium res, # 1 for low res

    print("Finished all of the preparation for the spectral fits. At this point you should copy\n"
          "the data folder over to a high performance computer. If speed and space aren't an issue,\n"
          "copy the entire cluster folder over. The only files really needed are the configuration\n"
          "file [cluster_pypeline_config.ini] and the acb folder (maintaining the original directory\n"
          "structure).")

    print("\nNext, you can either run wrapper.py directly or create a small bash file for scheduled\n"
          "running. If you need to create the bash file to call wrapper.py, make sure you setup the\n"
          "ciao environment in the bash file before calling python wrapper.py --configfile --etc.")

def eff_times_to_fits(clstr: cluster.ClusterObj):
    from astropy.io import fits
    for observation in clstr.observations:
        effbt = observation.effective_background_time
        effdt = observation.effective_data_time

        temp = fits.open(clstr.scale_map_file)

        temp[0].data = effbt
        print(io.get_path("writing {}/{}_{}_eff_bkg_time.fits".format(clstr.directory, clstr.name, observation.id)))
        temp.writeto(io.get_path('{}/{}_{}_eff_bkg_time.fits'.format(clstr.directory, clstr.name, observation.id)))

        temp[0].data = effdt
        print(io.get_path('writing {}/{}_{}_eff_data_time.fits'.format(clstr.directory, clstr.name, observation.id)))
        temp.writeto(io.get_path('{}/{}_{}_eff_data_time.fits'.format(clstr.directory, clstr.name, observation.id)))


def make_density_map(clstr: cluster.ClusterObj):
    from astropy.io import fits

    xray_sb_fits = fits.open(clstr.xray_surface_brightness_nosrc_filename)
    xray_sb_header = xray_sb_fits[0].header

    xray_sb = xray_sb_fits[0].data
    # xray surface brightness is proportional to density squared.
    # therefore, relative density is the square root of the xray surface brightness
    n = np.sqrt(xray_sb)
    fits.writeto(clstr.density_map_filename, n, header=xray_sb_header, overwrite=True)

    return n


def make_pressure_map(clstr: cluster.ClusterObj):

    from astropy.io import fits
    n = make_density_map(clstr)
    T = clstr.temperature_map
    header = clstr.temperature_map_header

    n, T = do.make_sizes_match(n, T)

    P = n*T

    norm_P = do.normalize_data(P)

    fits.writeto(clstr.pressure_map_filename, norm_P, header=header, overwrite=True)


if __name__ == '__main__':
    args, parser = get_arguments()
    if args.cluster_config is not None:
        if io.file_exists(args.cluster_config):
            clstr = cluster.read_cluster_data(args.cluster_config)
        else:
            clstr_config = cluster.get_cluster_config(args.cluster_config)
            if io.file_exists(clstr_config):
                clstr = cluster.read_cluster_data(clstr_config)
            else:
                print("Error finding cluster configuration file. Try passing the full path and filename.")
                exit(1)

        if args.commands:
            make_commands_lis(clstr, args.resolution)
        if args.eff_times_fits:
            eff_times_to_fits(clstr)
        elif args.temperature_map:
            print("Creating temperature map.")
            make_temperature_map(clstr, args.resolution)
        elif args.pressure:
            make_pressure_map(clstr)
        elif args.shock:
            import shockfinder
            if io.check_yes_no("This is likely not going to work. Continue?"):
                shockfinder.find_shock_in(clstr)
        else:
            fitting_preparation(clstr, args)
    else:
        parser.print_help()


