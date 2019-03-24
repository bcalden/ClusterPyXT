import pypeline_io as io
import os
import config
import subprocess
import numpy as np
import acb
import cluster
import sys
import data_operations as do

try:
    from ciao_contrib.cda.data import download_chandra_obsids
    import ciao_contrib.logger_wrapper as lw

    lw.initialize_logger("download", verbose=1)
    from ciao_contrib import runtool as rt
except ImportError:
    print("Failed to import CIAO python scripts. ")
    raise


def download_data(cluster):
    # refactor to make use of the observation class

    io.set_working_directory(cluster.directory)
    print("Downloading observation id(s) {}".format(cluster.observation_ids))

    success = []
    if not isinstance(cluster.observation_ids, list):
        cluster.observation_ids = [cluster.observation_ids]
    for observation in cluster.observation_ids:
        print("Downloading data from observation {}".format(observation))

        success.append(download_chandra_obsids([observation]))
        if success[-1]:
            print("Successfully downloaded data for observation {}.".format(observation))
            cluster.observation(observation).set_ccds()
            #observation.set_ccds()
        else:
            print("Failed trying to download data for observation {}.".format(observation))

    return not (False in success)


def prepare_to_merge_observations_from(cluster_obj):
    print("Preparing to merge the observations.")

    for observation in cluster_obj.observations:
        print("Preparing observation: {id}".format(id=observation.id))
        io.make_directory(observation.analysis_directory)
        io.make_directory(observation.reprocessing_directory)
        io.copytree(src=observation.primary_directory, dst=observation.analysis_directory)
        io.copytree(src=observation.secondary_directory, dst=observation.analysis_directory)

    return
    #
    # for observation in cluster_obj.observation_ids:
    #     observation_dir = "{}/{}".format(cluster_obj.directory, observation)
    #     print("\nObservation directory: {}".format(observation_dir))
    #     analysis_dir = "{}/analysis/".format(observation_dir)
    #     io.make_directory("{}/analysis/repro".format(observation_dir))
    #     for obs_type in ['primary', 'secondary']:
    #         dir_to_copy = "{}/{}".format(observation_dir, obs_type)
    #         print(dir_to_copy)
    #         io.copytree(dir_to_copy, analysis_dir)
    # return


def unzip_data_from(cluster_obj):
    print("Unzipping data for {}.".format(cluster_obj.name))
    for observation in cluster_obj.observation_list:
        print("Unzipping {}/{}.".format(cluster_obj.name, observation))
        current_analysis_dir = "{}/{}/analysis/".format(cluster_obj.directory, observation)
        io.set_working_directory(current_analysis_dir)
        all_files_in_directory = os.listdir()

        files_to_unzip = [gz_file for gz_file in all_files_in_directory if gz_file.endswith('.gz')]
        for gz_file in files_to_unzip:
            print("Unzipping {}/{}/{}.".format(cluster_obj.name, observation, gz_file))
            io.gz_unzip(gz_file)
    return


def dates_and_versions_match(acis_filename, background_filename):
    acis = {'date': io.get_date_from_filename(acis_filename),
            'version': io.get_version_from_filename(acis_filename)}

    background = {'date': io.get_version_from_filename(background_filename),
                  'version': io.get_version_from_filename(background_filename)}

    return (acis['date'] == background['date']) and (acis['version'] == background['version'])


def acis_process_events(gainfile, infile, outfile, clobber=False):
    rt.acis_process_events.punlearn()

    kwargs = {'infile': infile,
              'outfile': outfile,
              'acaofffile': None,
              'stop': None,
              'doevtgrade': False,
              'apply_cti': True,
              'apply_tgain': False,
              'calculate_pi': True,
              'pix_adj': None,
              'gainfile': gainfile,
              'clobber': clobber,
              'eventdef': '{s:ccd_id,s:node_id,i:expno,s:chip,s:tdet,f:det,f:sky,s:'
                          'phas,l:pha,l:pha_ro,f:energy,l:pi,s:fltgrade,s:grade,x:status}'}

    rt.acis_process_events(**kwargs)

    return None


def reprocess(cluster, observation, acis_gain, background_gain, acis_id):
    local_path = "{cluster_dir}/{observation}/analysis/".format(cluster_dir=cluster.directory,
                                                                observation=observation)
    outfile = "{dir}/back_newgain_ccd{acis_id}.fits".format(dir=local_path,
                                                            acis_id=acis_id)

    infile = "{dir}/back_ccd{acis_id}.fits".format(dir=local_path,
                                                   acis_id=acis_id)

    gainfile = "{ciao_dir}/CALDB/data/chandra/acis/det_gain/{acis_gain}".format(ciao_dir=config.ciao_directory(),
                                                                                acis_gain=acis_gain)
    print("Reprocessing {cluster}/{observation}/{acis_id}".format(cluster=cluster.name,
                                                                  observation=observation,
                                                                  acis_id=acis_id))

    acis_process_events(gainfile=gainfile, infile=infile, outfile=outfile)

    return outfile


def ciao_back(cluster, overwrite=False):
    print("Running ciao_back on {}.".format(cluster.name))

    for observation in cluster.observations:
        pcad_file = make_pcad_lis(cluster, observation.id)
        backI_lis = []
        backS_lis = []
        analysis_path = observation.analysis_directory
        filelist = io.read_contents_of_file(observation.ccd_merge_list).split('\n')
        pcad = io.read_contents_of_file(pcad_file)
        for acis_file in filelist:
            rt.acis_bkgrnd_lookup.punlearn()
            print("Finding background for {}".format(acis_file))
            path_to_background = rt.acis_bkgrnd_lookup(infile=acis_file)
            print("Found background at {}".format(path_to_background))
            acis_id = int(acis_file.split('/')[-1].split('.')[-2][-1])
            assert isinstance(acis_id, int), "acis_id = {}".format(acis_id)
            assert not isinstance(path_to_background, type(None)), "Cannot find background {}".format(acis_file)

            local_background_path = io.get_path("{}/back_ccd{}.fits".format(analysis_path, acis_id))
            try:
                if io.file_exists(local_background_path) and overwrite:
                    io.delete(local_background_path)
                io.copy(path_to_background, local_background_path)
            except OSError:
                print("Problem copying background file {}. Do you have the right permissions and a full CALDB?".format(
                    path_to_background))

                raise

            acis_gain = rt.dmkeypar(infile=acis_file,
                                    keyword="GAINFILE",
                                    echo=True)
            background_gain = rt.dmkeypar(infile=local_background_path,
                                          keyword="GAINFILE",
                                          echo=True)

            print("{}/{}/acis_ccd{}.fits gain: {}".format(cluster.name, observation.id, acis_id, acis_gain))
            print("{}/{}/back_ccd{}.fits gain: {}".format(cluster.name, observation.id, acis_id, background_gain))

            if dates_and_versions_match(acis_gain, background_gain):
                print("Date/version numbers don't match on the acis data and background. Reprocessing.")

                local_background_path = reprocess(cluster, observation.id, acis_gain, background_gain, acis_id)

            print("Reprojecting background")
            rt.reproject_events.punlearn()
            infile = local_background_path
            outfile = io.get_path("{local_path}/back_reproj_ccd{acis_id}.fits".format(local_path=analysis_path,
                                                                                      acis_id=acis_id))
            match = acis_file

            print(
                "Running:\n reproject_events(infile={infile}, outfile={outfile}, aspect={pcad}, match={match})".format(
                    infile=infile, outfile=outfile, pcad=pcad, match=match)
            )
            rt.reproject_events(infile=infile,
                                outfile=outfile,
                                aspect="{pcad_file}".format(pcad_file=pcad),
                                match=match,
                                random=0,
                                clobber=True)

            back_reproject = outfile
            datamode = rt.dmkeypar(infile=io.get_filename_matching(io.get_path("{}/acis*evt1*.fits".format(analysis_path))),
                                   keyword="DATAMODE")
            if datamode == "VFAINT":
                print("VFAINT Mode, resetting setting status bits")
                rt.dmcopy.punlearn()
                rt.dmcopy(infile="{}[status=0]".format(back_reproject),
                          outfile=outfile,
                          clobber=True)
            if acis_id <= 3:
                backI_lis.append(back_reproject)
            else:
                backS_lis.append(back_reproject)

        merged_back_list = backI_lis + backS_lis

        print("writing backI.lis and backS.lis")
        io.write_contents_to_file("\n".join(backI_lis), io.get_path("{}/backI.lis".format(analysis_path)),
                                  binary=False)
        io.write_contents_to_file("\n".join(backS_lis), io.get_path("{}/backS.lis".format(analysis_path)),
                                  binary=False)
        io.write_contents_to_file("\n".join(merged_back_list), observation.merged_back_lis, binary=False)

    return


def check_chandra_repro(observation):
    if observation.original_reprocessed_evt2_file_exists and observation.original_reprocessed_bad_pixel_file_exists:
        return True
    else:
        print("Error during reprocessing. CIAO command chandra_repro produced no results for {}.".format(observation.id))
        print("Likely an error during chandra_repro. Try deleting the cluster directory and")
        print("starting the cluster over. If problem persists, please post an issue on GitHub.")
        sys.exit(1)


def reprocess_cluster(cluster):
    print("Reprocessing {}".format(cluster.name))
    for observation in cluster.observations:
        print("Reprocessing {}/{}".format(cluster.name, observation.id))
        chandra_repro(indir=observation.directory, outdir=observation.reprocessing_directory)
        print("Checking to ensure {} was reprocessed.".format(observation.id))
        check_chandra_repro(observation)
        print("Reprocessing appears succesful. Copying reprocessed event files.")
        copy_reprocessed_event_files(observation)
    return


def ciao_merge_stack(stack_lis_file):
    rt.dmmerge.punlearn()
    infile = "@{}[subspace -expno]".format(stack_lis_file)
    outfile = io.change_extension(stack_lis_file, "fits")

    rt.dmmerge(infile=infile,
               outfile=outfile,
               outBlock="",
               columnList="",
               clobber=True)

    return outfile


def ciao_merge_background(cluster):
    for observation in cluster.observations:
        print("Merging background files from {}/{}".format(cluster.name, observation.id))

        merged = ciao_merge_stack(observation.merged_back_lis)

        print("Merged background written to {}".format(merged))
        ciao_hiE_sources(observation)


def chandra_repro(indir="./", outdir="./repro", set_ardlib=False, clobber=True):
    indir = io.get_path(indir)
    outdir = io.get_path(outdir)

    kwargs = {'indir': indir,
              'outdir': outdir,
              'set_ardlib': set_ardlib,
              'clobber': clobber,
              'verbose': 1}

    rt.chandra_repro.punlearn()

    rt.chandra_repro(**kwargs)

    return


def copy_reprocessed_event_files(observation):
    evt2_filename = observation.original_reprocessed_bad_pixel_filename
    bpix1_filename = observation.original_reprocessed_bad_pixel_filename

    io.copy(evt2_filename, observation.reprocessed_evt2_filename)
    io.copy(bpix1_filename, observation.reprocessed_bad_pixel_filename)

    print("Copied level 2 event files")
    return None


def ccd_sort(cluster):
    print("Running ccd_sort on {}.".format(cluster.name))
    for observation in cluster.observations:
        print("Working on {}/{}".format(cluster.name, observation.id))
        analysis_path = observation.analysis_directory
        os.chdir(analysis_path)
        evt1_filename = io.get_path("{}/{}".format(analysis_path,
                                                            io.get_filename_matching("acis*evt1.fits")[0]))
        evt2_filename = io.get_path("{}/{}".format(analysis_path,
                                                            io.get_filename_matching("evt2.fits")[0]))
        detname = rt.dmkeypar(infile=evt1_filename, keyword="DETNAM", echo=True)
        print("evt1 : {}\nevt2 : {}\ndetname : {}".format(evt1_filename,
                                                          evt2_filename,
                                                          detname))
        assert not isinstance(detname, type(None)), "detname returned nothing!"
        detnums = [int(x) for x in detname.split('-')[-1]]

        for acis_id in detnums:
            print("{cluster}/{observation}: Making level 2 event file for ACIS Chip id: {acis_id}".format(
                cluster=cluster.name,
                observation=observation.id,
                acis_id=acis_id))
            rt.dmcopy(infile="{evt2_file}[ccd_id={acis_id}]".format(evt2_file=evt2_filename,
                                                                    acis_id=acis_id),
                      outfile="acis_ccd{acis_id}.fits".format(acis_id=acis_id),
                      clobber=True)

        acisI_list = io.get_filename_matching("acis_ccd[0-3].fits")
        for i in range(len(acisI_list)):
            acisI_list[i] = io.get_path("{obs_analysis_dir}/{file}".format(obs_analysis_dir=observation.analysis_directory,
                                                                           file=acisI_list[i]))
        io.write_contents_to_file("\n".join(acisI_list), observation.ccd_merge_list, binary=False)
        merge_data_and_backgrounds(cluster, acisI_list)

    return


def merge_data_and_backgrounds(cluster, acis_list):
    rt.dmmerge.punlearn()

    merged_file = "acisI.fits"
    rt.dmmerge(infile="@acisI.lis[subspace -expno]",
               outfile=merged_file,
               clobber=True)

    detname = rt.dmkeypar(infile=io.get_filename_matching("acis*evt1.fits"),
                          keyword="DETNAM")
    # acisI3 = detname.find("3")
    # acisS3 = detname.find("7")

    rt.dmlist.punlearn()
    rt.dmlist(infile=merged_file,
              opt="header")

    return None


def actually_merge_observations_from(cluster):
    print("Merging observations from {}.".format(cluster.name))

    merged_directory = io.get_path('{}/merged_obs_evt2/'.format(cluster.directory))

    io.make_directory(merged_directory)

    os.chdir(merged_directory)
    merged_observations = []

    for observation in cluster.observation_list:
        evt2_file = "{}/{}/analysis/evt2.fits".format(cluster.directory, observation)
        merged_observations.append(evt2_file)

    merged_lis = "{}/merged_obs.lis".format(merged_directory)
    io.write_contents_to_file("\n".join(merged_observations), merged_lis, binary=False)
    outroot = io.get_path("{}/{}".format(cluster.directory,cluster.name))

    infile = "@{infile}[ccd_id=0:3]".format(infile=merged_lis) # for ACIS-I
    # infile = "@{infile}".format(infile=merged_lis) # for ACIS-I & ACIS-S

    xygrid = "1500:6500:4,1500:6500:4"

    if len(merged_observations) == 1:
        rt.fluximage.punlearn()
        rt.fluximage(infile=infile,
                     outroot=outroot,
                     xygrid=xygrid,
                     clobber=True)
        print("Only single observation, flux image created.")

    elif len(merged_observations) > 1:
        rt.merge_obs.punlearn()
        rt.merge_obs(infiles=infile,
                     outroot=outroot,
                     xygrid=xygrid,
                     clobber=True,
                     parallel=True,
                     nproc=12)


def ciao_hiE_sources(observation):
    #  I don't know that this function is actually worthwhile
    data = observation.data_filename
    #background = observation.back_filename

    #cts_min = '625'
    # print("Finding optimal binning fators for ACIS images")
    #acis_bin = 4
    rt.dmcopy.punlearn()
    infile = "{infile}".format(infile=data)
    outfile = "{}/img_acisI_fullE.fits".format(observation.analysis_directory)
    rt.dmcopy(infile=infile, outfile=outfile, clobber=True)


# final part of runpipe1
def merge_observations(cluster):
    prepare_to_merge_observations_from(cluster)
    unzip_data_from(cluster)
    reprocess_cluster(cluster)
    ccd_sort(cluster)
    ciao_back(cluster)
    ciao_merge_background(cluster)
    actually_merge_observations_from(cluster)
    return


def sources_file_exists(cluster):
    return os.path.isfile(cluster.sources_file)


def remove_sources_from_observation(observation):
    print("removing sources from {}".format(observation.id))

    # remove sources from foreground and background
    fore_or_back = [observation.data_filename, observation.back_filename]

    for i, type_of_obs in enumerate(fore_or_back):
        infile = "{type_of_obs}[exclude sky=region({sources})]".format(
            type_of_obs=type_of_obs,
            sources=observation.cluster.sources_file
        )
        outfile = [observation.acis_nosrc_filename, observation.background_nosrc_filename][i]
        clobber = True

        print("infile: {}".format(infile))
        print("outfile: {}".format(outfile))

        rt.dmcopy.punlearn()
        rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)
        if type_of_obs is observation.background_nosrc_filename:
            print("Copying background to {}".format(observation.back))
            io.copy(outfile, observation.back)


def remove_sources(observation):
    if sources_file_exists(observation.cluster):
        remove_sources_from_observation(observation)


def generate_light_curve(observation):

    # filter out high energy background flares
    obsid_analysis_dir = observation.analysis_directory
    data = observation.acis_nosrc_filename
    background = observation.background_nosrc_filename

    infile = "{}[energy=9000:12000]".format(data)
    outfile = "{}/acisI_hiE.fits".format(obsid_analysis_dir)
    clobber = True

    rt.dmcopy.punlearn()
    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    data_hiE = outfile
    infile = "{}[bin sky=8]".format(data_hiE)
    outfile = "{}/img_acisI_hiE.fits".format(obsid_analysis_dir)

    rt.dmcopy.punlearn()
    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    backbin = 259.28

    echo = True
    tstart = rt.dmkeypar(infile=data_hiE, keyword="TSTART", echo=echo)
    tstop = rt.dmkeypar(infile=data_hiE, keyword="TSTOP", echo=echo)

    print("Creating a lightcurve from the high energy events list with dmextract")

    rt.dmextract.punlearn()
    infile = "{}[bin time={}:{}:{}]".format(data_hiE, tstart, tstop, backbin)
    outfile = "{}/acisI_lcurve_hiE.lc".format(obsid_analysis_dir)

    print('Running dmextract infile={} outfile={} opt=ltc1 clobber=True'.format(infile, outfile))

    rt.dmextract(infile=infile,
                 outfile=outfile,
                 opt='ltc1', clobber=True)

    lcurve_hiE = outfile

    print("cleaning the lightcurve for {}, press enter to continue.".format(observation.id))

    rt.deflare.punlearn()

    outfile = "{}/acisI_gti_hiE.gti".format(obsid_analysis_dir)
    method = "clean"
    save = "{}/acisI_lcurve_hiE".format(obsid_analysis_dir)

    rt.deflare(infile=lcurve_hiE, outfile=outfile, method=method, save=save)

    gti_hiE = outfile

    print("Filtering the event list using GTI info from high energy flares.")

    infile = "{}[@{}]".format(data, gti_hiE)
    outfile = "{}/acisI_nosrc_hiEfilter.fits".format(obsid_analysis_dir)

    print("running: dmcopy infile={} outfile={} clobber={}".format(infile, outfile, clobber))

    rt.dmcopy.punlearn()
    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    data_nosrc_hiEfilter = outfile

    infile = "{}[bin sky=8]".format(data_nosrc_hiEfilter)
    outfile = "{}/img_acisI_nosrc_fullE.fits".format(obsid_analysis_dir)

    rt.dmcopy.punlearn()

    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)


def lightcurves_with_exclusion(cluster):
    for observation in cluster.observations:


        # data_nosrc_hiEfilter = "{}/acisI_nosrc_fullE.fits".format(obs_analysis_dir)

        data_nosrc_hiEfilter = "{}/acisI_nosrc_hiEfilter.fits".format(observation.analysis_directory)

        print("Creating the image with sources removed")

        data = observation.acis_nosrc_filename

        image_nosrc = "{}/img_acisI_nosrc_fullE.fits".format(observation.analysis_directory)

        if io.file_exists(observation.exclude_file):
            print("Removing sources from event file to be used in lightcurve")

            infile = "{}[exclude sky=region({})]".format(data_nosrc_hiEfilter, observation.exclude)
            outfile = "{}/acisI_lcurve.fits".format(observation.analysis_directory)
            clobber = True

            rt.dmcopy.punlearn()
            rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

            data_lcurve = "{}/acisI_lcurve.fits".format(observation.analysis_directory)
        else:
            yes_or_no = io.check_yes_no(
                "Are there sources to be excluded from observation {} while making the lightcurve? ".format(observation.id))

            if yes_or_no:  # yes_or_no == True
                print("Create the a region file with the region to be excluded and save it as {}".format(observation.exclude_file))
            else:
                data_lcurve = data_nosrc_hiEfilter

        backbin = 259.28

        echo = True
        tstart = rt.dmkeypar(infile=data_nosrc_hiEfilter, keyword="TSTART", echo=echo)
        tstop = rt.dmkeypar(infile=data_nosrc_hiEfilter, keyword="TSTOP", echo=echo)

        print("Creating lightcurve from the events list with dmextract")

        infile = "{}[bin time={}:{}:{}]".format(data_lcurve, tstart, tstop, backbin)
        outfile = "{}/acisI_lcurve.lc".format(observation.analysis_directory)
        opt = "ltc1"

        rt.dmextract.punlearn()
        rt.dmextract(infile=infile, outfile=outfile, opt=opt, clobber=clobber)

        lcurve = outfile

        print("Cleaning the lightcurve by removing flares with deflare. Press enter to continue.")

        rt.deflare.punlearn()
        infile = lcurve
        outfile = "{}/acisI_gti.gti".format(observation.analysis_directory)
        method = "clean"
        save = "{}/acisI_lcurve".format(observation.analysis_directory)

        rt.deflare(infile=infile, outfile=outfile, method=method, save=save)

        gti = outfile

        print("filtering the event list using GTI info just obtained.")

        infile = "{}[@{}]".format(data_nosrc_hiEfilter, gti)
        outfile = observation.clean
        clobber = True

        rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

        data_clean = outfile

        print("Don't forget to check the light curves!")


def sources_and_light_curves(cluster):
    for observation in cluster.observations:
        remove_sources(observation)
        generate_light_curve(observation)


def create_global_response_file_for(cluster, obsid, region_file):
    observation = cluster.observation(obsid)
    #min_counts = 525

    obs_analysis_dir = observation.analysis_directory
    global_response_dir = "{}/globalresponse/".format(obs_analysis_dir)
    io.make_directory(global_response_dir)

    clean = observation.clean
    back = observation.back

    pbk0 = io.get_filename_matching("{}/acis*pbk0*.fits".format(obs_analysis_dir))[0]
    bad_pixel_file = io.get_filename_matching("{}/bpix1_new.fits".format(obs_analysis_dir))[0]

    rt.ardlib.punlearn()

    rt.acis_set_ardlib(badpixfile=bad_pixel_file)

    mask_file = io.get_filename_matching("{}/*msk1.fits".format(obs_analysis_dir))

    make_pcad_lis(cluster, obsid)

    infile = "{}[sky=region({})]".format(clean, region_file)
    outroot = "{}/acisI_region_0".format(global_response_dir)
    weight = True
    correct_psf = False
    pcad = "@{}/pcad_asol1.lis".format(obs_analysis_dir)
    combine = False
    bkg_file = ""
    bkg_resp = False
    group_type = "NUM_CTS"
    binspec = 1
    clobber = True

    rt.specextract(infile=infile, outroot=outroot, weight=weight, correctpsf=correct_psf,
                   asp=pcad, combine=combine, mskfile=mask_file, bkgfile=bkg_file, bkgresp=bkg_resp,
                   badpixfile=bad_pixel_file, grouptype=group_type, binspec=binspec, clobber=clobber)

    infile = "{}[sky=region({})][bin pi]".format(back, region_file)
    outfile = "{}/acisI_back_region_0.pi".format(global_response_dir)
    clobber = True

    rt.dmextract.punlearn()
    print("Running: dmextract infile={}, outfile={}, clobber={}".format(infile, outfile, clobber))
    rt.dmextract(infile=infile, outfile=outfile, clobber=clobber)

    rt.dmhedit.punlearn()
    infile = "{}/acisI_region_0.pi".format(global_response_dir)
    filelist = ""
    operation = "add"
    key = "BACKFILE"
    value = outfile

    rt.dmhedit(infile=infile, filelist=filelist, operation=operation, key=key, value=value)

    observation = cluster.observation(obsid)

    aux_response_file = '{global_response_directory}/acisI_region_0.arf'.format(
        global_response_directory=observation.global_response_directory)

    redist_matrix_file = '{global_response_directory}/acisI_region_0.rmf'.format(
        global_response_directory=observation.global_response_directory)

    io.copy(aux_response_file, observation.aux_response_file)
    io.copy(redist_matrix_file, observation.redistribution_matrix_file)


def make_pcad_lis(cluster, obsid):
    analysis_dir = cluster.obs_analysis_directory(obsid)
    search_str = "{}/*asol1.fits".format(analysis_dir)
    pcad_files = io.get_filename_matching(search_str)
    pcad_list_string = "\n".join(pcad_files)
    pcad_filename = "{}/pcad_asol1.lis".format(analysis_dir)

    io.write_contents_to_file(pcad_list_string, pcad_filename, binary=False)

    return pcad_filename


def make_response_files(cluster):
    for obsid in cluster.observation_ids:
        print("Making response files for observation {}".format(obsid))
        obs_analysis_dir = cluster.obs_analysis_directory(obsid)
        region_file = io.get_path("{}/acisI_region_0.reg".format(obs_analysis_dir))

        if (not io.file_exists(region_file)) or (io.file_size(region_file) == 0):
            print("Region file {} does not exist.".format(region_file))
            print("When DS9 opens, draw a small circle that covers a piece of each ACIS-I chip (~20 pixels) and save it as:\n" \
                  "{}".format(region_file))
            print("Opening SAO DS9")
            io.write_contents_to_file("", region_file, False)
            ds9_arguments = "ds9 -regions system physical -regions shape circle -regions format ciao -zoom 0.5 " \
                            "-bin factor 4 {}/acisI_clean.fits".format(obs_analysis_dir)
            subprocess.run([ds9_arguments], shell=True)
        print('Creating global response file.')
        create_global_response_file_for(cluster, obsid, region_file)


def make_mask_file(observation):
    from astropy.io import fits
    print("Creating an image mask for {}.".format(observation.id))

    original_fits_filename = observation.acisI_comb_img

    mask = fits.open(original_fits_filename)

    print("{} shape: {}".format(original_fits_filename, mask[0].shape))
    mask[0].data = np.ones_like(mask[0].data)

    mask_filename = observation.temp_acis_comb_mask_filename

    mask.writeto(mask_filename, overwrite=True)

    rt.dmcopy.punlearn()
    # infile = "{mask_filename}[sky=region({fov_file})][opt full]".format( # for ACIS-I & ACIS-S
    infile = "{mask_filename}[sky=region({fov_file}[ccd_id=0:3])][opt full]".format(  # for ACIS-I
        mask_filename=mask_filename,
        fov_file=observation.fov_file
    )
    outfile = observation.acisI_combined_mask
    clobber = True

    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    print("Image mask created for {obsid} and saved as {filename}".format(
        obsid=observation.id, filename=outfile
    ))

    io.delete(observation.temp_acis_comb_mask_filename)


def make_cumulative_mask_file(cluster, observation):
    from astropy.io import fits
    cumulative_mask_filename = cluster.combined_mask

    current_obs_mask_filename = observation.acisI_combined_mask

    if not io.file_exists(cumulative_mask_filename):
        print("Cumulative mask file not found. Creating it.")
        cumulative_mask = fits.open(current_obs_mask_filename)
        cumulative_mask.writeto(cumulative_mask_filename)
    else:
        current_mask = fits.open(current_obs_mask_filename)
        cumulative_mask = fits.open(cumulative_mask_filename)

        print("Cumulative mask {} shape:{}".format(cumulative_mask_filename,
                                                   cumulative_mask[0].shape))
        print("current mask {} shape:{}".format(current_obs_mask_filename,
                                                   current_mask[0].shape))

        current_mask[0].data = current_mask[0].data + cumulative_mask[0].data

        current_mask[0].data[np.where(current_mask[0].data > 1)] = 1
        current_mask.writeto(cumulative_mask_filename, overwrite=True)


def make_acisI_and_back_for(observation, cluster):
    from astropy.io import fits

    rt.dmcopy.punlearn()
    rt.dmcopy(
        infile="{clean_file}[sky=region({mask})]".format(
            clean_file=observation.clean, mask=cluster.master_crop_file),
        outfile=cluster.temp_acisI_comb,
        clobber=True
    )

    shp = fits.open(cluster.temp_acisI_comb)[0].shape
    print(observation.clean)
    print("{} shape {}".format(cluster.temp_acisI_comb,
                               shp))

    rt.dmcopy.punlearn()
    rt.dmcopy(
        infile="{temp_acisI_combined}[bin sky=4][energy=700:8000]".format(
            temp_acisI_combined=cluster.temp_acisI_comb),
        outfile=observation.acisI_comb_img,
        clobber=True
    )

    shp = fits.open(observation.acisI_comb_img)[0].shape

    print("{} shape {}".format(observation.acisI_comb_img,
                               shp))

    rt.dmcopy.punlearn()
    rt.dmcopy(
        infile="{back_file}[sky=region({mask})]".format(
            back_file=observation.back,
            mask=cluster.master_crop_file),
        outfile=cluster.temp_backI_comb,
        clobber=True
    )

    rt.dmcopy.punlearn()
    rt.dmcopy(
        infile="{temp_backI_combined}[bin sky=4][energy=700:8000]".format(temp_backI_combined=cluster.temp_backI_comb),
        outfile=observation.backI_comb_img,
        clobber=True
    )

    io.delete(cluster.temp_acisI_comb)
    io.delete(cluster.temp_backI_comb)


def run_ds9_for_master_crop(cluster):
    print("Need to create a box region containing all parts of the image you want included.")
    print("Save this file as {master_crop}.".format(master_crop=cluster.master_crop_file))
    ds9_arguments = "ds9 -regions system physical -regions shape circle -regions format ciao -zoom 0.5 " \
                    "-bin factor 4 {cluster_dir}/{name}_broad_flux.img".format(cluster_dir=cluster.directory,
                                                                              name=cluster.name)
    subprocess.run([ds9_arguments], shell=True)


def runpipe5(cluster):
    print("runpipe5")
    from astropy.io import fits
    # This portion of the pypeline 
    combined_dir = cluster.combined_directory

    io.make_directory(combined_dir)

    while not io.file_exists(cluster.master_crop_file):
        print("Master crop file not found")
        run_ds9_for_master_crop(cluster)

    # the contents of this for loop should be refactored/replaced with the make_acisI_and_back function
    for observation in cluster.observations:
        infile = "{}[sky=region({})]".format(observation.clean, cluster.master_crop_file)
        outfile = cluster.temp_acisI_comb
        clobber = True

        rt.dmcopy.punlearn()
        rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

        print("{} shape: {}".format(outfile, fits.open(outfile)[0].shape))

        infile = "{}[bin sky=4][energy=700:8000]".format(cluster.temp_acisI_comb)
        outfile = observation.acisI_comb_img
        clobber = True

        print("ObsID: {}\t- Extracting just 0.7keV - 8keV.".format(observation.id))
        rt.dmcopy.punlearn()
        rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

        # background
        infile = "{}[sky=region({})]".format(observation.back, cluster.master_crop_file)
        outfile = cluster.temp_backI_comb
        clobber = True

        rt.dmcopy.punlearn()
        rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

        infile = "{}[bin sky=4][energy=700:8000]".format(cluster.temp_backI_comb)
        outfile = observation.backI_comb_img
        clobber = True

        rt.dmcopy.punlearn()
        rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

        io.delete(cluster.temp_acisI_comb)
        io.delete(cluster.temp_backI_comb)

        make_mask_file(observation)
        make_cumulative_mask_file(cluster, observation)

    create_combined_images(cluster)
    make_nosrc_cropped_xray_sb(cluster)

def create_combined_images(cluster):
    from astropy.io import fits
    print("Combining images")
    mask = fits.open("{combined_dir}/acisI_comb_mask.fits".format(
        combined_dir=cluster.combined_directory
    ))

    counts_image = np.zeros(mask[0].data.shape)
    back_rescale = np.zeros(mask[0].data.shape)

    t_obs = t_back = 0.0

    good_crop = False
    while not good_crop:
        completed_obs = 0
        for obsid in cluster.observation_ids:
            print("Working on observation id {obsid}.".format(obsid=obsid))

            obs_img = fits.open("{combined_dir}/acisI_comb_img-{obsid}.fits".format(
                combined_dir=cluster.combined_directory,
                obsid=obsid))

            print("obs_img = {combined_dir}/acisI_comb_img-{obsid}.fits".format(combined_dir=cluster.combined_directory,
                                                                                obsid=obsid))
            print("Counts_image shape = {}".format(counts_image.shape))
            print("obs_img.shape = {}".format(obs_img[0].shape))

            if counts_image.shape != obs_img[0].data.shape:
                print("Padding the observation in the periphery to match dimensions.")
                counts_image, obs_img[0].data = do.make_sizes_match(counts_image, obs_img[0].data)

                #
                # # make use of the make_sizes_match()
                # print(
                #     'The shapes don\'t match. Remake the {master_crop} in DS9 and rerun this step of the pypeline.'.format(
                #         master_crop=cluster.master_crop_file))
                # # raise
                # run_ds9_for_master_crop(cluster)
                # break

            counts_image += obs_img[0].data

            t_obs = obs_img[0].header['EXPOSURE']

            print("Type of t_obs is {}".format(type(t_obs)))

            back_img = fits.open("{combined_dir}/backI_comb_img-{obsid}.fits".format(
                combined_dir=cluster.combined_directory,
                obsid=obsid
            ))

            t_back = back_img[0].header['EXPOSURE']

            print("Type of t_back is {}".format(type(t_back)))
            print("Type of back_img[0].data = {}".format(back_img[0].data.dtype))
            back_rescale += (t_obs / t_back) * (back_img[0].data.astype(float))
            completed_obs += 1

        if completed_obs == len(cluster.observation_ids):
            good_crop = True


    signal = counts_image - back_rescale

    obs_img[0].data = signal
    obs_img.writeto(cluster.combined_signal, overwrite=True)

    obs_img[0].data = counts_image
    obs_img.writeto(cluster.counts_image, overwrite=True)

    obs_img[0].data = back_rescale
    obs_img.writeto(cluster.back_rescale, overwrite=True)

    print("Images combined!")


def extract_spec(observation, region_file, region_number, dtime, btime):
    infile = "{clean}[sky=region({region_file})][bin pi]".format(
        clean=observation.sc_clean,
        region_file=region_file
    )

    outfile = io.get_path(
        "{super_comp_dir}/{obsid}_{region_number}.pi".format(
            super_comp_dir=observation.cluster.super_comp_dir,
            obsid=observation.id,
            region_number=region_number
        ))

    rt.dmextract(infile=infile, outfile=outfile, clobber=True)

    infile = "{back}[sky=region({region_file})][bin pi]".format(
        back=observation.sc_back,
        region_file=region_file
    )

    outfile = io.get_path(
        "{super_comp_dir}/{obsid}_back_{region_number}.pi".format(
            super_comp_dir=observation.cluster.super_comp_dir,
            obsid=observation.id,
            region_number=region_number
        ))

    rt.dmextract(infile=infile, outfile=outfile, clobber=True)

    data_pi = "{super_comp_dir}/{obsid}_{region_number}.pi".format(
        super_comp_dir=observation.cluster.super_comp_dir,
        obsid=observation.id,
        region_number=region_number
    )

    back_pi = "{super_comp_dir}/{obsid}_back_{region_number}.pi".format(
        super_comp_dir=observation.cluster.super_comp_dir,
        obsid=observation.id,
        region_number=region_number
    )

    warf = "'{super_comp_dir}/{name}_{obsid}.arf'".format(
        super_comp_dir=observation.cluster.super_comp_dir,
        name=observation.cluster.name,
        obsid=observation.id
    )

    wrmf = "'{super_comp_dir}/{name}_{obsid}.rmf'".format(
        super_comp_dir=observation.cluster.super_comp_dir,
        name=observation.cluster.name,
        obsid=observation.id
    )

    # Put this background file into the 'grouped' data file for the region

    #rt.dmhedit(infile=data_pi, filelist="", operation="add", key="BACKFILE", value=back_pi)

    rt.dmhedit(infile=data_pi, filelist="", operation="add", key="EXPOSURE", value=dtime)
    rt.dmhedit(infile=data_pi, filelist="", operation="add", key="RESPFILE", value=wrmf)
    rt.dmhedit(infile=data_pi, filelist="", operation="add", key="ANCRFILE", value=warf)
    rt.dmhedit(infile=data_pi, filelist="", operation="add", key="BACKFILE", value=back_pi)
    rt.dmhedit(infile=back_pi, filelist="", operation="add", key="EXPOSURE", value=btime)

    io.append_to_file(observation.cluster.spec_lis(region_number), "{}\n".format(data_pi))

    return (data_pi, back_pi)


def spec_extract(observation, region_file, region_num, min_counts):
    infile = "{data_clean}[sky=region({region_file})]".format(
        data_clean=observation.cluster.acisI_clean_obs(observation.id),
        region_file=region_file
    )

    outroot = "{outdir}/acisI_region_{obsid}_{region_num}".format(
        outdir=observation.cluster.super_comp_dir,
        region_num=region_num,
        obsid=observation.id
    )

    print("Starting specextract on {}".format(infile))

    rt.specextract(infile=infile,
                   outroot=outroot,
                   weight='yes',
                   correct='no',
                   asp='@{}'.format(observation.pcad_asol),
                   combine='no',
                   mskfile=observation.acis_mask_sc,
                   #                   bkgfile=observation.back,
                   bkgresp="no",
                   badpixfile=observation.bad_pixel_file,
                   binspec=1,
                   clobber=True
                   )
    io.append_to_file(observation.cluster.spec_lis(region_num), "{}.pi\n".format(outroot))


def get_keyword_value(filename, keyword):
    value = rt.dmkeypar(infile=filename, keyword=keyword, echo=True)
    # print("{infile}['{keyword}'] = {value}".format(
    #     infile=filename,
    #     keyword=keyword,
    #     value=value
    # ))
    return value


def get_exposure(filename):
    return get_keyword_value(filename, "EXPOSURE")

def run_stage_1(cluster):
    download_data(cluster)
    merge_observations(cluster)


def start_from_last(cluster):
    complete_string = ['',
                       'Data downloaded. Now time to merge the observations',
                       'Merging complete!']

    # if cluster_filename is None:
    #     cluster_obj.initialize_cluster()
    #     kwargs = {'current_cluster_name': cluster_obj.name,
    #               'cluster_config_file': cluster_obj.configuration_filename}
    #     config.update_system_configuration(**kwargs)
    # else:
    #     cluster_obj = cluster.read_cluster_data(cluster_filename)

    continue_running = True
    while continue_running:
        pypeline_progress_index = int(cluster.last_step_completed)
        function_steps = [None,
                          download_data, # runpipe1 / stage 1
                          merge_observations, #runpipe1 / stage 1
                          sources_and_light_curves, #runpipe2 / stage 2
                          lightcurves_with_exclusion, #runpipe3 / stage 2
                          make_response_files, #runpipe4 / stage 3
                          runpipe5, #runpipe5
                          acb.fitting_preparation] #runpipe acb

        if pypeline_progress_index < len(function_steps):
            success = function_steps[pypeline_progress_index](cluster)

            if pypeline_progress_index == 1:
                if success:
                    print(complete_string[pypeline_progress_index])
                    cluster.last_step_completed = 2
                else:
                    print("Problem downloading the data. Try again")
                    continue_running = False
            elif pypeline_progress_index == 2:
                cluster.last_step_completed = 3
                print("Now time to select the sources. The broad_flux.img file in the main observation directory"
                      "is the file you need to create the source region file from. You can do this by opening "
                      "that file in SAO DS9 and selecting them by hand, or using an automated tool. Save the sources"
                      "file in the main observation directory as sources.reg and then rerun the pypeline.")


                #start DS9, load the broad flux FITS image, and a sources.reg file for the user to begin
                continue_running = False
            elif pypeline_progress_index == 3:
                print("Sources removed.")
                continue_running = False
                print("draw regions around sources to exclude when creating the lightcurve "
                      "and save as exclude.reg in the clusters main directory.")
                cluster.last_step_completed = 4
            elif pypeline_progress_index == 4:
                print("Lightcurve creation complete")
                continue_running = False
                cluster.last_step_completed = 5
            elif pypeline_progress_index == 5:
                cluster.last_step_completed = 6
                continue_running = True
            elif pypeline_progress_index == 6:
                cluster.last_step_completed = 7
                continue_running = False
            elif pypeline_progress_index == 7:
                cluster.last_step_completed = 8
                continue_running = False

        else:
            print("That functionality is not yet complete. Help by contributing on GitHub!")
            continue_running = False
    return

def initialize_cluster(name="", obsids=[], abundance=0.3, redshift=0.0, nH=0.0):
    clstr = cluster.ClusterObj(name=name, observation_ids=obsids, abundance=abundance,
                               redshift=redshift, hydrogen_column_density=nH,
                               data_directory=config.data_directory())
    print('Making initial cluster directory: {}'.format(clstr.directory))
    io.make_directory(clstr.directory)
    io.make_initial_directories(clstr)
    clstr.last_step_completed = 1
    print("Downloading cluster data.")
    download_data(clstr)
    clstr.last_step_completed = 2
    print("Merging observations.")
    merge_observations(clstr)
    clstr.last_step_completed = 3


def automated_cluster_init(batch_file):
    print("Automated cluster initialization using: {batch_file}".format(batch_file=batch_file))
    data_directory = config.data_directory()
    csv_clusters = io.get_cluster_info_from_csv(batch_file)
    for clstr in csv_clusters:
        cluster_obj = cluster.ClusterObj(name=clstr['name'],
                                         observation_ids=clstr['obsids'],
                                         data_directory=data_directory,
                                         abundance=clstr['abundance'],
                                         redshift=clstr['redshift'],
                                         hydrogen_column_density=clstr['hydrogen_column_density']
                                         )

        io.make_directory(cluster_obj.directory)
        cluster_obj.write_cluster_data()
        io.make_initial_directories(cluster_obj)
        cluster_obj.last_step_completed = 1
        download_data(cluster_obj)
        cluster_obj.last_step_completed = 2
        merge_observations(cluster_obj)
        cluster_obj.last_step_completed = 3


def copy_image_excluding_region(image_file, sources_file, output_file, overwrite=False):
    infile = "{file}[exclude sky=region({sources})]".format(
        file=image_file,
        sources=sources_file
    )
    outfile = output_file
    copy_image(infile, outfile, overwrite=overwrite)


def copy_image(infile, outfile, overwrite=False):
    rt.dmcopy.punlearn()
    rt.dmcopy(infile=infile, outfile=outfile, clobber=overwrite)

def copy_image_cropping_region(image_file, crop_region, output_file, overwrite=False):
    infile = "{file}[sky=region({crop_region})]".format(
        file=image_file,
        crop_region=crop_region
    )
    outfile = output_file
    copy_image(infile, outfile, overwrite=overwrite)


def make_cropped_xray_sb_image(clstr: cluster.ClusterObj):
    infile = "{file}[sky=region({master_crop})]".format(
        file=clstr.xray_surface_brightness_filename,
        master_crop=clstr.master_crop_file
    )
    outfile = clstr.xray_surface_brightness_cropped_filename
    copy_image(infile, outfile, overwrite=True)

def remove_sources_from_xray_surface_brightness(clstr: cluster.ClusterObj):
    copy_image_excluding_region(clstr.xray_surface_brightness_cropped_filename,
                                clstr.sources_file,
                                clstr.xray_surface_brightness_nosrc_filename,
                                overwrite=True)


def make_nosrc_cropped_xray_sb(clstr: cluster.ClusterObj):
    make_cropped_xray_sb_image(clstr)
    remove_sources_from_xray_surface_brightness(clstr)

