from enum import IntEnum
from errors import ClusterPyError
import pypeline_io as io
import os
import config
import subprocess
import numpy as np
import acb
import cluster
import sys
import data_operations as do
import time
import multiprocessing as mp
import astropy.io.fits as fits
import spectral
from tqdm import tqdm

try:
    from ciao_contrib.cda.data import download_chandra_obsids
    import ciao_contrib.logger_wrapper as lw

    lw.initialize_logger("download", verbose=1)
    from ciao_contrib import runtool as rt
except ImportError:
    print("Failed to import CIAO python scripts. Is CIAO running?")
    sys.exit(1)


class Stage(IntEnum):
    zero = 0
    one = 1
    two = 2
    three = 3
    four = 4
    five = 5
    tmap = 6


def download_obsid(obsid):
    # print(f'Downloading {obsid}')
    return download_chandra_obsids([obsid])

def download_data(cluster):
    io.set_working_directory(cluster.directory)

    obsids = [int(obsid) if obsid != '' else None for obsid in cluster.observation_ids]
    
    num_cpus = mp.cpu_count() // 2
    if num_cpus > 5:
        num_streams = 3
    else:
        num_streams = num_cpus // 2

    with mp.Pool(num_streams) as pool:
        results = pool.map(download_obsid, obsids)

    _ = [cluster.observation(obsid).set_ccds() for obsid in obsids]
    return results


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

    gainfile = "{ciao_dir}/CALDB/data/chandra/acis/det_gain/{acis_gain}".format(ciao_dir=config.sys_config.ciao_directory,
                                                                                acis_gain=acis_gain)
    print("Reprocessing {cluster}/{observation}/{acis_id}".format(cluster=cluster.name,
                                                                  observation=observation,
                                                                  acis_id=acis_id))

    acis_process_events(gainfile=gainfile, infile=infile, outfile=outfile)

    return outfile


def ciao_back(cluster, overwrite=False):
    print("Running ciao_back on {}.".format(cluster.name))

    for observation in cluster.observations:
        pcad_file = make_pcad_lis(observation)
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
                io.copy(path_to_background, local_background_path, replace=True)
            except OSError:
                print("Problem copying background file {}. Do you have the right permissions and a full CALDB?".format(
                    path_to_background))

                raise
            try:
                rt.dmkeypar.punlearn()
                print(f'Running dmkeypar {acis_file} "GAINFILE" echo=True')
                acis_gain = rt.dmkeypar(infile=acis_file,
                                        keyword="GAINFILE",
                                        echo=True)
                rt.dmkeypar.punlearn()
                print(f'Running dmkeypar {local_background_path} "GAINFILE" echo=True')
                background_gain = rt.dmkeypar(infile=local_background_path,
                                            keyword="GAINFILE",
                                            echo=True)
            except ValueError:
                
                print("Error getting parameter file in CIAO. Please close ClusterPyXT and re-try the stage. If the problem persists, please file a bug report on https://github.com/bcalden/ClusterPyXT with the following error message:")
                raise
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
                                aspect="@{pcad_file}".format(pcad_file=pcad_file),
                                match=match,
                                random=0,
                                clobber=True)

            back_reproject = outfile
            datamode = rt.dmkeypar(infile=observation.level_1_event_filename,
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


def reprocess_cluster_multiobs(cluster: cluster.ClusterObj):
    print("Reprocessing {}.".format(cluster.name))
    result = chandra_repro_multi(cluster)
    print(result)

    return


def reprocess_cluster(cluster):
    print("Reprocessing {}. This may take a a little while (potentially 10s of minutes)".format(cluster.name))

    for observation in cluster.observations:
        print("Reprocessing {}/{}".format(cluster.name, observation.id))
        result = chandra_repro(observation)
        print(result)
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


def chandra_repro(observation: cluster.Observation):

    rt.chandra_repro.punlearn()
    os.chdir(observation.analysis_directory)
    output = rt.chandra_repro(indir=observation.directory,
                              outdir=observation.reprocessing_directory,
                              set_ardlib=False,
                              clobber=True,
                              verbose=1)

    return output

def chandra_repro_multi(cluster: cluster.ClusterObj):
    rt.chandra_repro.punlearn()
    os.chdir(cluster.directory)
    obsids = ",".join(cluster.observation_ids)
    output = rt.chandra_repro(indir=obsids, outdir="", set_ardlib=False, clobber=True, verbose=1)

    return output

def ccd_sort(cluster):
    print("Running ccd_sort on {}.".format(cluster.name))
    for observation in cluster.observations:
        print("Working on {}/{}".format(cluster.name, observation.id))
        #analysis_path = observation.analysis_directory
        #os.chdir(analysis_path)
        evt1_filename = observation.level_1_event_filename
        evt2_filename = observation.reprocessed_evt2_filename
        detname = rt.dmkeypar(infile=evt1_filename, keyword="DETNAM", echo=True)
        print("evt1 : {}\nevt2 : {}\ndetname : {}".format(evt1_filename,
                                                          evt2_filename,
                                                          detname))
        assert not isinstance(detname, type(None)), "detnam keyword not in level 1 event file: {}".format(
            observation.level_1_event_filename
        )
        detnums = [int(x) for x in detname.split('-')[-1]]

        io.make_directory(observation.analysis_directory)

        for acis_id in detnums:
            print("{cluster}/{observation}: Making level 2 event file for ACIS Chip id: {acis_id}".format(
                cluster=cluster.name,
                observation=observation.id,
                acis_id=acis_id))
            try:
                rt.dmcopy(infile=observation.reprocessed_evt2_for_ccd(acis_id),
                          outfile=observation.acis_ccd(acis_id),
                          clobber=True)
            except OSError as oserr:
                print("Error generating event files for each CCD.")
                print("Observation: {}\t CCD: {}".format(observation.id, acis_id))
                print("File: {}".format(observation.reprocessed_evt2_for_ccd(acis_id)))
                if not io.file_sizes_match(observation.reprocessed_evt2_filename, observation.original_reprocessed_evt2_filename):
                    print("File sizes don't match for  {} and {}.".format(observation.reprocessed_evt2_filename,
                                                                          observation.original_reprocessed_evt2_filename))
                    print("These should be the same file. Try manually copying {og} to {new} and retrying.".format(
                        og=observation.original_reprocessed_evt2_filename,
                        new=observation.reprocessed_evt2_filename
                    ))
                print("Retry last pipeline step. If problem persists, please post an issue to GitHub.")
                raise
                #sys.exit(1)

        os.chdir(observation.analysis_directory)
        if observation.acis_type == 0:  # ACIS-I
            acis_list = io.get_filename_matching("acis_ccd[0-3].fits")
        elif observation.acis_type == 1:  # ACIS-S
            acis_list = io.get_filename_matching("acis_ccd[4-8].fits")
        for i in range(len(acis_list)):
            acis_list[i] = io.get_path("{obs_analysis_dir}/{file}".format(obs_analysis_dir=observation.analysis_directory,
                                                                           file=acis_list[i]))
        io.write_contents_to_file("\n".join(acis_list), observation.ccd_merge_list, binary=False)
        merge_data_and_backgrounds(cluster, acis_list)

    return


def merge_data_and_backgrounds(cluster, acis_list):
    # merges the
    rt.dmmerge.punlearn()

    merged_file = "acisI.fits" # needs to be renamed to reflect ACIS-I & S
    rt.dmmerge(infile="@acisI.lis[subspace -expno]",
               outfile=merged_file,
               clobber=True)

    #detname = rt.dmkeypar(infile=io.get_filename_matching("acis*evt1.fits"),
    #                      keyword="DETNAM")
    # acisI3 = detname.find("3")
    # acisS3 = detname.find("7")

    rt.dmlist.punlearn()
    rt.dmlist(infile=merged_file,
              opt="header")

    return None


def actually_merge_observations_from(cluster):
    print("Merging observations from {}.".format(cluster.name))

    merged_directory = cluster.merged_directory

    io.make_directory(merged_directory)

    os.chdir(merged_directory)
    merged_observations = []

    for observation in cluster.observations:
        merged_observations.append(observation.ccd_filtered_reprocessed_evt2_filename)

    merged_lis = "{}/merged_obs.lis".format(merged_directory)
    io.write_contents_to_file("\n".join(merged_observations), merged_lis, binary=False)
    outroot = io.get_path("{}/{}/".format(cluster.directory, cluster.name))

    infile = "@{infile}".format(infile=merged_lis) # for ACIS-I & ACIS-S

    xygrid = "-3000:10000:4,-3000:10000:4"

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
                     binsize=4,
                     #xygrid=xygrid,
                     clobber=True,
                     parallel=True,
                     nproc=12)


def make_point_spread_function_map(observation, ecf=0.3, energy=1.4):
    rt.mkpsfmap.punlearn()
    rt.mkpsfmap(infile=observation.broad_threshold_image_filename,
                outfile=observation.point_spread_function_map_filename,
                energy=energy,
                ecf=ecf,
                clobber=True)


def wav_detect(observation):
    rt.wavdetect.punlearn()
    rt.wavdetect(infile=observation.broad_threshold_image_filename,
                 outfile=observation.source_map_filename,
                 scellfile=observation.source_cell_map_filename,
                 imagefile=observation.source_image_filename,
                 defnbkgfile=observation.normalized_background_without_sources_filename,
                 scales="1.0 2.0 4.0 8.0 16.0",
                 psffile=observation.point_spread_function_map_filename,
                 regfile=observation.source_region_filename,
                 clobber=True)


def merge_source_files(cluster: cluster.ClusterObj):
    region_files = [observation.source_region_filename for observation in cluster.observations]
    io.merge_region_files(region_files, cluster.sources_file)


def find_sources(cluster: cluster.ClusterObj, ecf=0.1, energy=1.4):
    for observation in cluster.observations:
        print("Finding sources in {}".format(observation.id))
        make_point_spread_function_map(observation, ecf=ecf, energy=energy)
        wav_detect(observation)
    merge_source_files(cluster)


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

    reprocess_cluster_multiobs(cluster) #multiobs test
    ccd_sort(cluster)
    ciao_back(cluster)
    ciao_merge_background(cluster)
    actually_merge_observations_from(cluster)
    return


def sources_file_exists(cluster):
    return os.path.isfile(cluster.sources_file)


def remove_sources_from_observation(observation):
    # print("removing sources from {}".format(observation.id))

    # remove sources from foreground and background
    fore_or_back = [observation.data_filename, observation.back_filename]

    for i, type_of_obs in enumerate(fore_or_back):
        infile = "{type_of_obs}[exclude sky=region({sources})]".format(
            type_of_obs=type_of_obs,
            #sources=observation.cluster.sources_file
            sources=observation.source_region_filename
        )
        outfile = [observation.acis_nosrc_filename, observation.background_nosrc_filename][i]
        clobber = True

        # print("infile: {}".format(infile))
        # print("outfile: {}".format(outfile))

        rt.dmcopy.punlearn()
        rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber, verbose=0)
        if type_of_obs is observation.background_nosrc_filename:
            # print("Copying background to {}".format(observation.back))
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

    # print("Creating a lightcurve from the high energy events list with dmextract")

    rt.dmextract.punlearn()
    infile = "{}[bin time={}:{}:{}]".format(data_hiE, tstart, tstop, backbin)
    outfile = "{}/acisI_lcurve_hiE.lc".format(obsid_analysis_dir)

    # print('Running dmextract infile={} outfile={} opt=ltc1 clobber=True'.format(infile, outfile))

    rt.dmextract(infile=infile,
                 outfile=outfile,
                 opt='ltc1', clobber=True)

    lcurve_hiE = outfile

    # print("cleaning the lightcurve for {}, press enter to continue.".format(observation.id))

    rt.deflare.punlearn()

    outfile = "{}/acisI_gti_hiE.gti".format(obsid_analysis_dir)
    method = "clean"
    save = "{}/acisI_lcurve_hiE".format(obsid_analysis_dir)

    rt.deflare(infile=lcurve_hiE, outfile=outfile, method=method, save=save)

    gti_hiE = outfile

    # print("Filtering the event list using GTI info from high energy flares.")

    infile = "{}[@{}]".format(data, gti_hiE)
    outfile = "{}/acisI_nosrc_hiEfilter.fits".format(obsid_analysis_dir)

    # print("running: dmcopy infile={} outfile={} clobber={}".format(infile, outfile, clobber))

    rt.dmcopy.punlearn()
    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    data_nosrc_hiEfilter = outfile

    infile = "{}[bin sky=8]".format(data_nosrc_hiEfilter)
    outfile = "{}/img_acisI_nosrc_fullE.fits".format(obsid_analysis_dir)

    rt.dmcopy.punlearn()

    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber, verbose=0)


def lightcurve_with_exclusion_for(observation):
    data_nosrc_hiEfilter = "{}/acisI_nosrc_hiEfilter.fits".format(observation.analysis_directory)

    # print("Creating the image with sources removed")

    data = observation.acis_nosrc_filename

    image_nosrc = "{}/img_acisI_nosrc_fullE.fits".format(observation.analysis_directory)

    if io.file_exists(observation.exclude_file):
        # print("Removing sources from event file to be used in lightcurve")

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

    infile = "{}[bin time={}:{}:{}]".format(data_lcurve, tstart, tstop, backbin)
    outfile = "{}/acisI_lcurve.lc".format(observation.analysis_directory)
    opt = "ltc1"

    rt.dmextract.punlearn()
    rt.dmextract(infile=infile, outfile=outfile, opt=opt, clobber=clobber)

    lcurve = outfile

    rt.deflare.punlearn()
    infile = lcurve
    outfile = "{}/acisI_gti.gti".format(observation.analysis_directory)
    method = "clean"
    save = "{}/acisI_lcurve".format(observation.analysis_directory)

    rt.deflare(infile=infile, outfile=outfile, method=method, save=save)

    gti = outfile

    infile = "{}[@{}]".format(data_nosrc_hiEfilter, gti)
    outfile = observation.clean
    clobber = True
    
    rt.dmcopy.punlearn()
    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    data_clean = outfile

# def lightcurves_with_exclusion(cluster:cluster.ClusterObj, args):
#     num_obs = len(cluster.observations)
#     num_runs = (num_obs // args.num_cpus) + 1
#     obs_lists = np.array_split(cluster.observation_ids, num_runs)
#     for obs_list in tqdm(obs_lists):
#         processes = [mp.Process(target=lightcurve_with_exclusion_for, args=(cluster.observation(obsid),)) for obsid in obs_list]

#         for process in processes:
#             process.start()
#         for process in processes:
#             process.join()


def lightcurves_with_exclusion(cluster:cluster.ClusterObj, args):
    for observation in tqdm(cluster.observations, desc='Finishing light curves', unit='observation', total=len(cluster.observations)):
        lightcurve_with_exclusion_for(observation)
    

def sources_and_light_curves(cluster):
    print("Source removal:")
    for observation in cluster.observations:
        print("Removing sources for {obsid}".format(obsid=observation.id))
        remove_sources(observation)

    print("Light curves:")
    for observation in cluster.observations:
        print("Generating light curves for {obsid}".format(obsid=observation.id))
        generate_light_curve(observation)

def remove_sources_in_parallel(cluster, args):
    # Remove point sources for each observation in parallel.
    try:
        with mp.Pool(args.num_cpus) as pool:
            _ = list(tqdm(pool.imap(remove_sources, cluster.observations), 
                        total=len(cluster.observations), 
                        desc='Removing point sources', 
                        unit='observations')
                        )
    except:
        print(f"Error in removing sources in parallel. CPU count:{args.num_cpus}")
        print(f"Trying again with single core.")

        for observation in tqdm(cluster.observations, 
                            total=len(cluster.observations), 
                            desc='Removing point sources', 
                            unit='observation'):
            remove_sources(observation)


def generate_light_curves(cluster, args):
    # Generate light curves for each obsertion. Each observation gets its own process
    for observation in tqdm(cluster.observations, total=len(cluster.observations), desc='Generating light curves', unit='observation'):
        generate_light_curve(observation)

    ## For Parallel operation below - sometimes crashes (50/50)
    # num_obs = len(cluster.observations)
    # num_runs = (num_obs // args.num_cpus) + 1
    # obs_lists = np.array_split(cluster.observation_ids, num_runs)
    # for obs_list in tqdm(obs_lists, desc='Generating light curves in batches', unit='batch'):
    #     processes = [mp.Process(target=generate_light_curve, args=(cluster.observation(obsid), )) for obsid in obs_list]
    #     for process in processes:
    #         process.start()
    #     for process in processes:
    #         process.join()


def create_global_response_file_for(observation: cluster.Observation):

    #min_counts = 525

    obs_analysis_dir = observation.analysis_directory
    global_response_dir = observation.global_response_directory
    io.make_directory(global_response_dir)

    bad_pixel_file = observation.reprocessed_bad_pixel_filename
    clean = observation.clean
    back = observation.back

    rt.ardlib.punlearn()

    rt.acis_set_ardlib(badpixfile=bad_pixel_file)

    mask_file = observation.mask_file

    make_pcad_lis(observation)

    infile = "{}[sky=region({})]".format(clean, observation.response_file_region_covering_ccds)
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

    rt.specextract.punlearn()

    start_time = time.time()
    # print("Running specextract on {}".format(observation.id))
    # print("Size of region0: {}".format(observation.acisI_region_0_size))
    rt.specextract(infile=infile, outroot=outroot, weight=weight, correctpsf=correct_psf,
                   asp=pcad, combine=combine, mskfile=mask_file, bkgfile=bkg_file, bkgresp=bkg_resp,
                   badpixfile=bad_pixel_file, grouptype=group_type, binspec=binspec, clobber=clobber)
    elapsed = time.time() - start_time
    # print("Elapsed time: {:0.2f} sec(s)".format(elapsed))

    infile = "{}[sky=region({})][bin pi]".format(back, observation.response_file_region_covering_ccds)
    outfile = "{}/acisI_back_region_0.pi".format(global_response_dir)
    clobber = True

    rt.dmextract.punlearn()
    # print("Running dmextract")
    #print("Running: dmextract infile={}, outfile={}, clobber={}".format(infile, outfile, clobber))
    start_time = time.time()
    rt.dmextract(infile=infile, outfile=outfile, clobber=clobber)
    elapsed = time.time() - start_time
    # print("Elapsed time: {:0.2f} sec(s)".format(elapsed))

    rt.dmhedit.punlearn()
    infile = "{}/acisI_region_0.pi".format(global_response_dir)
    filelist = ""
    operation = "add"
    key = "BACKFILE"
    value = outfile

    rt.dmhedit(infile=infile, filelist=filelist, operation=operation, key=key, value=value)

    aux_response_file = '{global_response_directory}/acisI_region_0.arf'.format(
        global_response_directory=observation.global_response_directory)

    redist_matrix_file = '{global_response_directory}/acisI_region_0.rmf'.format(
        global_response_directory=observation.global_response_directory)

    io.copy(aux_response_file, observation.aux_response_file)
    io.copy(redist_matrix_file, observation.redistribution_matrix_file)


def make_pcad_lis(observation: cluster.Observation):
    search_str = "{}/*asol1.fits".format(observation.reprocessing_directory)
    pcad_files = io.get_filename_matching(search_str)
    pcad_list_string = "\n".join(pcad_files)
    pcad_filename = "{}/pcad_asol1.lis".format(observation.analysis_directory)

    io.write_contents_to_file(pcad_list_string, pcad_filename, binary=False)

    return pcad_filename


def do_function_on_observations_in_parallel(cluster: cluster.ClusterObj,
                                            function=None,
                                            num_cpus=1):
    observation_lists = cluster.parallel_observation_lists(num_cpus)
    num_observations = len(cluster.observations)

    start_time = time.time()

    for observation_list in observation_lists:
        print("Working on {} of {} observations.".format(len(observation_list), num_observations))
        processes = [mp.Process(target=function, args=(observation,)) for observation in observation_list]

        for process in processes:
            process.start()

        for process in processes:
            process.join()

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Elapsed time: {:2f} seconds".format(elapsed_time))


def make_response_files_in_parallel(cluster: cluster.ClusterObj, args):
    with mp.Pool(args.num_cpus) as pool:
        _ = list(tqdm(pool.imap(create_global_response_file_for, cluster.observations), desc='Creating global response files', total=len(cluster.observations), unit='observation'))


def make_response_files(cluster):
    for observation in cluster.observations:
        region_file = observation.response_file_region_covering_ccds
        print("Checking for response region file (region file covering at least a piece of all ACIS-I CCDs) for {}".format(observation.id))
        if (not io.file_exists(region_file)) or (io.file_size(region_file) == 0):
            print("Region file {} does not exist.".format(region_file))
            print("When DS9 opens, draw a small circle that covers a piece of each ACIS-I chip (~20 pixels) and save it as:\n" \
                  "{}".format(region_file))
            print("Opening SAO DS9")
            io.write_contents_to_file("", region_file, False)
            ds9_arguments = "ds9 -regions system physical -regions shape circle -regions format ciao -zoom 0.5 " \
                            "-bin factor 4 {clean_obs}".format(clean_obs=observation.clean)
            subprocess.run([ds9_arguments], shell=True)
        print('Creating global response file.')
        create_global_response_file_for(observation)


def make_mask_file(observation: cluster.Observation):
    from astropy.io import fits
    # print("Creating an image mask for {}.".format(observation.id))

    original_fits_filename = observation.acisI_comb_img

    mask = fits.open(original_fits_filename)

    # print("{} shape: {}".format(original_fits_filename, mask[0].shape))
    mask[0].data = np.ones_like(mask[0].data)

    mask_filename = observation.temp_acis_comb_mask_filename

    mask.writeto(mask_filename, overwrite=True)

    rt.dmcopy.punlearn()
    # need to check type of observation, S or I, and then generate a new string with the approriate ccd filtering

    if observation.acis_type == 0:  # ACIS-I
        ccd_filter = "0:3"
    else:
        ccd_filter = "4:9"

    #infile = "{mask_filename}[sky=region({fov_file})][opt full]".format( # for ACIS-I & ACIS-S
    infile = "{mask_filename}[sky=region({fov_file}[ccd_id={ccd_filter}])][opt full]".format(
        mask_filename=mask_filename,
        fov_file=observation.fov_file,
        ccd_filter=ccd_filter
    )
    outfile = observation.acisI_combined_mask_file
    clobber = True

    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    # print("Image mask created for {obsid} and saved as {filename}".format(
    #     obsid=observation.id, filename=outfile
    # ))

    io.delete(observation.temp_acis_comb_mask_filename)


def make_cumulative_mask_file(cluster, observation):
    cumulative_mask_filename = cluster.combined_mask

    current_obs_mask_filename = observation.acisI_combined_mask_file

    if not io.file_exists(cumulative_mask_filename):
        # print("Cumulative mask file not found. Creating it.")
        cumulative_mask = fits.open(current_obs_mask_filename)
        cumulative_mask.writeto(cumulative_mask_filename)
    else:
        current_mask = fits.open(current_obs_mask_filename)
        cumulative_mask = fits.open(cumulative_mask_filename)

        # print("Cumulative mask {} shape:{}".format(cumulative_mask_filename,
        #                                            cumulative_mask[0].shape))
        # print("current mask {} shape:{}".format(current_obs_mask_filename,
        #                                            current_mask[0].shape))
        try:
            cumulative_mask[0].data = current_mask[0].data + cumulative_mask[0].data
        except ValueError as err:
            # print("Shapes don't match, reprojecting image.")
            rt.reproject_image(infile=observation.acisI_combined_mask_file,
                               matchfile=cluster.combined_mask,
                               outfile=observation.temp_acis_comb_mask_filename)
            io.move(observation.temp_acis_comb_mask_filename, observation.acisI_combined_mask_file)

            # print("Combining masks")
            rt.dmimgcalc(infile=observation.acisI_combined_mask_file,
                         infile2=cluster.combined_mask,
                         operation='add',
                         outfile=cluster.combined_mask,
                         clobber=True)
            #cumulative_mask[0].data += fits.open(observation.temp_acis_comb_mask_filename)[0].data

            cumulative_mask[0].data = fits.open(cluster.combined_mask)[0].data
        cumulative_mask[0].data[np.where(cumulative_mask[0].data > 1)] = 1
        cumulative_mask.writeto(cumulative_mask_filename, overwrite=True)


def reproject(infile=None, matchfile=None, outfile=None, overwrite=False):
    rt.reproject_image.punlearn()
    rt.reproject_image(infile=infile, matchfile=matchfile, outfile=outfile, clobber=overwrite)


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


def stage_4_parallel(cluster: cluster.ClusterObj):
    print("Making observation masks.")
    do_function_on_observations_in_parallel(cluster, function=make_masks_for)

    print("Making the cumulative mask file.")
    make_cumulative_mask(cluster)


def make_cumulative_mask(cluster: cluster.ClusterObj):
    cumulative_mask_filename = cluster.combined_mask
    cumulative_mask = np.zeros(fits.open(cumulative_mask_filename)[0].data.shape)

    for observation in cluster.observations:
        current_obs_mask_filename = observation.acisI_combined_mask_file
        current_mask = observation.acisI_combined_mask
        try:
            cumulative_mask += current_mask
        except ValueError as err:
            print("Shapes don't match, reprojecting image.")
            rt.reproject_image(infile=observation.acisI_combined_mask_file,
                               matchfile=cluster.combined_mask,
                               outfile=observation.temp_acis_comb_mask_filename)
            io.move(observation.temp_acis_comb_mask_filename, observation.acisI_combined_mask_file)
            current_mask = fits.open(observation.acisI_combined_mask_file)

            print("Combining masks")
            cumulative_mask += current_mask

    # write any val > 1 = 1

    if not io.file_exists(cumulative_mask_filename):
        print("Cumulative mask file not found. Creating it.")
        cumulative_mask = fits.open(current_obs_mask_filename)
        cumulative_mask.writeto(cumulative_mask_filename)
    else:
        current_mask = fits.open(current_obs_mask_filename)
        cumulative_mask = fits.open(cumulative_mask_filename)


        cumulative_mask[0].data[np.where(cumulative_mask[0].data > 1)] = 1
        cumulative_mask.writeto(cumulative_mask_filename, overwrite=True)

def make_energy_filtered_image(observation: cluster.Observation):  # Energies in eV
    rt.dmcopy.punlearn()
    rt.dmcopy(infile=observation.cropped_clean_infile_string,
              outfile=observation.temp_acis_comb_filename,
              clobber=True)

    print("ObsID: {}\t- Extracting just 0.7keV - 8keV.".format(observation.id))
    rt.dmcopy.punlearn()
    rt.dmcopy(infile=observation.temporary_acis_combined_energy_filtered_infile_string,
              outfile=observation.acisI_comb_img,
              clobber=True)

    io.delete(observation.temp_acis_comb_filename)


def make_energy_filtered_background(observation: cluster.Observation):
    rt.dmcopy.punlearn()
    rt.dmcopy(infile=observation.cropped_background_infile_string,
              outfile=observation.temp_back_comb_filename,
              clobber=True)

    rt.dmcopy.punlearn()
    rt.dmcopy(infile=observation.temporary_back_combined_energy_filtered_infile_string,
              outfile=observation.backI_comb_img,
              clobber=True)

    io.delete(observation.temp_back_comb_filename)

def make_masks_for(observation: cluster.Observation):
    print("Making masks for: {}".format(observation.id))
    make_energy_filtered_image(observation)
    make_energy_filtered_background(observation)
    make_mask_file(observation)


def make_acisI_and_back(observation:cluster.Observation):
    cluster = observation.cluster
    infile = "{}[sky=region({})]".format(observation.clean, cluster.master_crop_file)
    outfile = cluster.temp_acisI_comb
    clobber = True

    rt.dmcopy.punlearn()
    rt.dmcopy(infile=infile, outfile=outfile, clobber=clobber)

    # print("{} shape: {}".format(outfile, fits.open(outfile)[0].shape))

    infile = "{}[bin sky=4][energy=700:8000]".format(cluster.temp_acisI_comb)
    outfile = observation.acisI_comb_img
    clobber = True

    # print("ObsID: {}\t- Extracting just 0.7keV - 8keV.".format(observation.id))
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

def stage_4(cluster: cluster.ClusterObj, args):
    combined_dir = cluster.combined_directory

    io.make_directory(combined_dir)

    while not io.file_exists(cluster.master_crop_file):
        print("Master crop file not found")
        run_ds9_for_master_crop(cluster)

    # the contents of this for loop should be refactored/replaced with the make_acisI_and_back function
    # do_function_on_observations_in_parallel(cluster, make_acisI_and_back, args.num_cpus)
    for observation in tqdm(cluster.observations, desc='Cropping observations', unit='observation', total=len(cluster.observations)):
        make_acisI_and_back(observation)
        
    create_combined_images(cluster)
    make_nosrc_cropped_xray_sb(cluster)


def  create_combined_images(cluster):
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
            observation = cluster.observation(obsid)
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
                io.copy(observation.acisI_combined_image_filename, observation.temp_acis_comb_filename)
                rt.reproject_image(infile=observation.temp_acis_comb_filename,
                                   matchfile=cluster.combined_mask,
                                   outfile=observation.acisI_combined_image_filename,
                                   clobber=True)
                obs_img = fits.open(observation.acisI_combined_image_filename)

            counts_image += obs_img[0].data

            t_obs = obs_img[0].header['EXPOSURE']

            print("Type of t_obs is {}".format(type(t_obs)))

            back_img = fits.open(observation.backI_comb_img)

            t_back = back_img[0].header['EXPOSURE']

            print("Type of t_back is {}".format(type(t_back)))
            print("Type of back_img[0].data = {}".format(back_img[0].data.dtype))
            try:
                back_rescale += (t_obs / t_back) * (back_img[0].data.astype(float))
            except ValueError as err:
                print("Shapes don't match. Reprojecting.")
                rt.reproject_image(infile=observation.backI_comb_img,
                                   matchfile=cluster.combined_mask,
                                   outfile=observation.backI_comb_temp_img,
                                   clobber=True)

                io.move(observation.backI_comb_temp_img, observation.backI_comb_img)

                back_img = fits.open(observation.backI_comb_img)

                back_rescale += (t_obs / t_back) * (back_img[0].data.astype(float))

            completed_obs += 1

        if completed_obs == len(cluster.observation_ids):
            good_crop = True


    signal = counts_image - back_rescale

    obs_img[0].data = signal
    obs_img.writeto(cluster.combined_signal, overwrite=True)

    obs_img[0].data = counts_image
    obs_img.writeto(cluster.counts_image_filename, overwrite=True)

    obs_img[0].data = back_rescale
    obs_img.writeto(cluster.back_rescale_filename, overwrite=True)

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
                   bkgfile=observation.back,
                   bkgresp="no",
                   badpixfile=observation.reprocessed_bad_pixel_filename,
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
    try:
        download_data(cluster)
    except TimeoutError as te:
        io.print_red("Download failed due to a timeout. Appears the connection was dropped. Please re-run stage 1.")
        io.print_red("Error: {te.strerror}")
        raise
    merge_observations(cluster)


def finish_stage_1(cluster: cluster.ClusterObj):
    print_stage_2_prep(cluster)


def print_stage_2_prep(cluster: cluster.ClusterObj):
    prep_msg = """Data downloaded and the observations are merged into a surface brightness map {sb_map_filename}. 
    Now it is time to filter out point sources and high energy flares. To do so, first open the surface brightness
    map and create regions around sources you want excluded from the data analysis. These are typically foreground
    point sources one does not want to consider when analyzing the cluster. Save these regions as a DS9 region file
    named {sources_file}. 

    Additionally, you need to create a region file containing any regions you wanted excluded from the deflaring process.
    This would include areas such as the peak of cluster emission as these regions may contain high energy events we want
    to consider in this analysis. Save this region file as {exclude_file}. 

    After both files are saved, you can continue ClusterPyXT on {cluster_name}""".format(
        sb_map_filename=cluster.xray_surface_brightness_filename,
        sources_file=cluster.sources_file,
        exclude_file=cluster.exclude_file,
        cluster_name=cluster.name
    )

    print(prep_msg)


def check_for_required_stage_2_files(cluster: cluster.ClusterObj):
    if io.file_exists(cluster.sources_file) and io.file_exists(cluster.exclude_file):
        return True
    else:
        io.print_red("Error: Missing {sources} and/or {exclude}".format(
            sources=cluster.sources_file,
            exclude=cluster.exclude_file
        ))
        print_stage_2_prep(cluster)

        #sys.exit(ClusterPyError.sources_or_exclude_not_found)


def run_stage_2(cluster):
    print("Starting Stage 2: {}".format(cluster.name))
    check_for_required_stage_2_files(cluster)
    sources_and_light_curves(cluster)
    make_nosrc_xray_sb(cluster)
    lightcurves_with_exclusion(cluster)

    return


def run_stage_2_parallel(cluster, args):
    
    check_for_required_stage_2_files(cluster)
    remove_sources_in_parallel(cluster,args)
    generate_light_curves(cluster, args)
    make_nosrc_xray_sb(cluster)
    lightcurves_with_exclusion(cluster, args)

    return


def finish_stage_2(cluster: cluster.ClusterObj):
    finish_str = """Stage 2 complete -  Point sources removed -> {xray_sb_nosrc}.
                                        High energy events filtered.""".format(
        xray_sb_nosrc=cluster.xray_surface_brightness_nosrc_filename
    )

    print(finish_str)
    print_stage_3_prep(cluster)


def print_stage_3_prep(cluster: cluster.ClusterObj):
    observation = cluster.observations[0]
    prep_str = """Next is stage 3. This stage extracts the RMF and ARF files. Before continuing the pipeline
    on {cluster_name}, you need to create a region file for each observation. Each observation
    will need its own region file named acisI_region_0.reg and saved in the respective analysis
    directory (e.g. {region_file}).
    
    To create this file, open the respective acisI_clean.fits file (e.g. {acisI_clean}) and draw
    a small circle region containing some of each of the ACIS-I CCD's. This region does not need
    to contain ALL of the chips, just a piece of each. It can be ~20 pixels (bigger circle=longer
    runtime). 
    
    After the region files for each observation are created, continue running ClusterPyXT on {cluster_name}""".format(
        xray_sb_nosrc=cluster.xray_surface_brightness_nosrc_filename,
        cluster_name=cluster.name,
        region_file=observation.response_file_region_covering_ccds,
        acisI_clean=observation.clean
    )

    print(prep_str)


def print_stage_3_file_message(cluster: cluster.ClusterObj):
    print("Need the response region file. Please check readme stage 3 for details.")


def check_for_required_stage_3_files(cluster: cluster.ClusterObj):
    all_files = True
    for observation in cluster.observations:
        if not io.file_exists(observation.response_file_region_covering_ccds):
            io.print_red("Cannot find: {}".format(observation.response_file_region_covering_ccds))
            all_files = False

    return all_files


def run_stage_3(cluster: cluster.ClusterObj, args):
    check_for_required_stage_3_files(cluster)
    make_response_files_in_parallel(cluster, args)


def finish_stage_3(cluster: cluster.ClusterObj):
    finish_str = """"""

    print(finish_str)
    print_stage_4_prep(cluster)


def print_stage_4_prep(cluster: cluster.ClusterObj):
    prep_str = """Now you need to create a region file enclosing the region you would like to crop
    the final analysis to. To do so, open the surface brightness file {xray_sb_file}
    and create a box region containing all parts of the image you want included in the analysis. 
    
    Save this file as: {master_crop}
    
    After this region file is created, continue running ClusterPyXT on {cluster_name}.
    
    Note: Due to processing complexities, during this stage you may encounter an error
     where two observations have slightly different dimensions (usually ~1) and the pipeline cannot
     combine them. This is due to the region splitting pixels and in some observations that pixel may
     be counted where in others it is not. If this happens, draw a new crop region, save it, and re-run.""".format(
        xray_sb_file=cluster.xray_surface_brightness_nosrc_filename,
        master_crop=cluster.master_crop_file,
        cluster_name=cluster.name
    )

    print(prep_str)


def run_stage_4(cluster: cluster.ClusterObj, args):
    stage_4(cluster, args)


def finish_stage_4(cluster: cluster.ClusterObj):
    finish_str = """Data filtered and cropped."""
    print(finish_str)

    print_stage_5_prep(cluster)


def print_stage_5_prep(cluster: cluster.ClusterObj):
    prep_str = """You are now ready for Stage 5. This stage only requires all previous stages to be completed. 
    Stage 5 calculates the adaptive circular bins, generates the scale map, and calculates exposure corrections. 
    It can take a long time (~10s of hours). 

    After stage 5 is complete, you are ready for spectral fitting.

    Please continue running ClusterPyXT on {cluster_name}.""".format(cluster_name=cluster.name)

    print(prep_str)


def run_stage_5(cluster: cluster.ClusterObj, args=None, num_cpus=mp.cpu_count()):
    acb.fitting_preparation(cluster, args, num_cpus=num_cpus)
    cluster.last_step_completed = Stage.five.value


def finish_stage_5(cluster: cluster.ClusterObj):
    finish_str = """Scale map created, adaptive circular bins generated, and various other files generated needed 
    to correct exposures."""

    print(finish_str)
    print_stage_spectral_fits_prep(cluster)


def print_stage_spectral_fits_prep(cluster: cluster.ClusterObj):
    prep_str = """Now ready for spectral fitting. This can be offloaded to a remote machine if necessary.
    If offloading, copy the cluster configuration file, {cluster_config},
    and the acb directory to the remote machine. Update the configuration file to reflect the appropriate
    path of your data.

    Next, with CIAO running, simply run:

        python spectral.py --parallel --resolution 2 -c {name}

    The resolution parameter can be set to either 1 - low resolution, 2 - medium resolution, or 3 - high resolution.
    If the parallel flag indicates to run in parallel. If the number of cpus is not set (--num_cpus), ClusterPyXT uses
    the total number of cores on your machine. 

    If you must restart the fitting for any reason, simply add the --continue flag in order to not redo any of the fits.""".format(
        cluster_config=cluster.configuration_filename,
        name=cluster.name
    )
    print(prep_str)


def run_stage_spectral_fits(cluster: cluster.ClusterObj, num_cpus=mp.cpu_count()):
    spectral.calculate_spectral_fits(cluster, num_cpus)

def finish_stage_spectral_fits(cluster: cluster.ClusterObj):
    print_stage_tmap_prep(cluster)


def print_stage_tmap_prep(cluster: cluster.ClusterObj):
    prep_str = """Now ready for spectral fitting. 

    If offloaded, copy the spectral fits file, {spectral_fits},
    back to the local machine. 

    Next, run 

        python acb.py --temperature_map --resolution 2 --cluster_config_file /path/to/cluster/A115/A115_pypeline_config.ini

    This will create the temperature map and allow for the creation of the pressure maps.""".format(
        spectral_fits=cluster.spec_fits_file
    )

def run_stage_tmap(cluster: cluster.ClusterObj):
    pass

def finish_stage_tmap(cluster: cluster.ClusterObj):
    pass


def start_from_last(cluster: cluster.ClusterObj, args=None):
    print("Continuing {}".format(cluster.name))

    last_stage_completed = int(cluster.last_step_completed)

    print("Last step completed: {}".format(last_stage_completed))

    if last_stage_completed == Stage.zero:
        run_stage_1(cluster)
        cluster.last_step_completed = Stage.one.value
        finish_stage_1(cluster)
        return

    elif last_stage_completed == Stage.one:
        run_stage_2_parallel(cluster, args)
        cluster.last_step_completed = Stage.two.value
        finish_stage_2(cluster)
        return

    elif last_stage_completed == Stage.two:
        run_stage_3(cluster, args.num_cpus)
        cluster.last_step_completed = Stage.three.value
        finish_stage_3(cluster)
        return

    elif last_stage_completed == Stage.three:
        run_stage_4(cluster, args)
        cluster.last_step_completed = Stage.four.value
        finish_stage_4(cluster)
        return

    elif last_stage_completed == Stage.four:
        run_stage_5(cluster, args=args)
        cluster.last_step_completed = Stage.five.value
        finish_stage_5(cluster)
        return

    elif last_stage_completed == Stage.five:
        run_stage_spectral_fits(cluster)
        cluster.last_step_completed = Stage.tmap.value
        finish_stage_spectral_fits(cluster)
        ### To be implemented
        # run_stage_tmap(cluster)
        # cluster.last_step_completed = Stage.tmap.value
        # finish_stage_tmap(cluster)
        return

    else:
        print("Error loading configuration file.")

    return


def initialize_cluster(name="", obsids=[], abundance=0.3, redshift=0.0, nH=0.0):
    clstr = cluster.ClusterObj(name=name, observation_ids=obsids, abundance=abundance,
                               redshift=redshift, hydrogen_column_density=nH,
                               data_directory=config.sys_config.data_directory)
    print('Making initial cluster directory: {}'.format(clstr.directory))
    io.make_directory(clstr.directory)
    io.make_initial_directories(clstr)
    #run_stage_1(clstr)


def automated_cluster_init(batch_file):
    print("Automated cluster initialization using: {batch_file}".format(batch_file=batch_file))
    data_directory = config.sys_config.data_directory
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


def remove_sources_from_cropped_xray_surface_brightness(clstr: cluster.ClusterObj):
    copy_image_excluding_region(image_file=clstr.xray_surface_brightness_cropped_filename,
                                sources_file=clstr.sources_file,
                                output_file=clstr.xray_surface_brightness_nosrc_cropped_filename,
                                overwrite=True)


def make_nosrc_xray_sb(clstr: cluster.ClusterObj):
    copy_image_excluding_region(image_file=clstr.xray_surface_brightness_filename,
                                sources_file=clstr.sources_file,
                                output_file=clstr.xray_surface_brightness_nosrc_filename,
                                overwrite=True)


def make_nosrc_cropped_xray_sb(clstr: cluster.ClusterObj):
    make_cropped_xray_sb_image(clstr)
    remove_sources_from_cropped_xray_surface_brightness(clstr)

