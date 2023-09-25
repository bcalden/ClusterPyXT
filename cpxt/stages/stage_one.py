"""
File: cpxt/stages/stage_one.py
Author: Brian Alden
Date: 30 Mar 2022
Description: This file contains the code to run stage one of `ClusterPyXT`. 
             Specifically, it downloads all the observation data entered for
             the particular cluster being processed. Next, we merge the 
             observations to create a merged surface brightness maps of all 
             observations. To finish, we find and merge the background files
             provided by CIAO in the CALDB. 
"""

# Internal imports
from cpxt.stages.stage_two import print_stage_2_prep
from cpxt.chandra.observation import Observation
from cpxt.core.stages import Stage, MessageType
from cpxt.chandra.cluster import Cluster
from cpxt.core.config import CPXTConfig
from cpxt.chandra.ccd import CCD
from cpxt.chandra import ciao
from cpxt.core import io

# External imports
from ciao_contrib import runtool as rt
import multiprocessing as mp
from pathlib import Path
from tqdm import tqdm
import shutil as sh
import logging
import os


logger = logging.getLogger(__name__)

################################################################################
################################################################################
##                                                                            ##
##              ███████╗████████╗ █████╗  ██████╗ ███████╗     ██╗            ##
##              ██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝    ███║            ##
##              ███████╗   ██║   ███████║██║  ███╗█████╗      ╚██║            ##
##              ╚════██║   ██║   ██╔══██║██║   ██║██╔══╝       ██║            ##
##              ███████║   ██║   ██║  ██║╚██████╔╝███████╗     ██║            ##
##              ╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝     ╚═╝            ##
##                                                                            ##
################################################################################
################################################################################

def run_on(cluster:Cluster) -> None:
    """Runs stage one of the pipeline on the provided cluster. First downloads 
    the observation data, reprocesses the observations to get them the same 
    full field of view, and generates the level 2 event file for each ACIS CCD. 
        
    Next we find and merges all the backgrounds found for this cluster in 
    the CALDB. 

    Finally, we actually merge the observations into one big X-ray surface 
    brightness map.

    Parameters
    ----------
    cluster : ClusterObj
        The ClusterPyXT ClusterObject of the particular cluster 
        you want to run stage one on.

    Returns
    -------
    None
    """
    download_xray_data_for(cluster)
    reprocess_observations_for(cluster)
    make_level_2_event_files_for(cluster)
    merge_foregrounds_for(cluster)
    reprocess_backgrounds_for(cluster)
    merge_observations_for(cluster)
    print_stage_1_comp(cluster)
    print_stage_2_prep(cluster)



################################################################################
################################################################################
##                                                                            ##
##   ██████╗  ██████╗ ██╗    ██╗███╗   ██╗██╗      ██████╗  █████╗ ██████╗    ##
##   ██╔══██╗██╔═══██╗██║    ██║████╗  ██║██║     ██╔═══██╗██╔══██╗██╔══██╗   ##
##   ██║  ██║██║   ██║██║ █╗ ██║██╔██╗ ██║██║     ██║   ██║███████║██║  ██║   ##
##   ██║  ██║██║   ██║██║███╗██║██║╚██╗██║██║     ██║   ██║██╔══██║██║  ██║   ##
##   ██████╔╝╚██████╔╝╚███╔███╔╝██║ ╚████║███████╗╚██████╔╝██║  ██║██████╔╝   ##
##   ╚═════╝  ╚═════╝  ╚══╝╚══╝ ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚═════╝    ##
##                                                                            ##
################################################################################
################################################################################

def download_xray_data_for(cluster:Cluster) -> bool:
    """Downloads the observation data for the provided cluster and the 
    relevant observation ids. 

    CIAO Commands executed: download_chandra_obsid

    Parameters
    ----------
    cluster : ClusterObj
        The ClusterPyXT ClusterObject of the particular cluster you want 
        to download data for.

    Returns
    -------
    bool
        True or false depending on successful data download
    """
    # Remove the primary and secondary observation directories to prevent
    # problems when we go to reprocess. 
    for observation in cluster.observations:
        try:
            sh.rmtree(observation.primary_dir)
            logger.info(f"Removing {observation.id} primary data")
        except FileNotFoundError:
            pass
        try:
            sh.rmtree(observation.secondary_dir)
            logger.info(f"Removing {observation.id} secondary data")
        except FileNotFoundError:
            pass
        
    # We must change directory to the clusters directory otherwise the 
    # observation downloads will go to whatever directory you are currently in.
    os.chdir(cluster.data_directory)
    
    logger.debug(f"Beginning to download {cluster.obsids}")

    results = ciao.download_observations(
        obsids=cluster.obsids, 
        num_streams=CPXTConfig().num_downloads
        )
    
    # Iterate through the list of results. If one is false, this step failed.
    for result in results:
        if not result:
            logger.error(f"Not all observations downloaded for {cluster.name}")
            return False    
    
    logger.info(f"Successfully downloaded observation data for {cluster.name}")
    
    return True

################################################################################
################################################################################
##                                                                            ##
## ██████╗ ███████╗██████╗ ██████╗  ██████╗  ██████╗███████╗███████╗███████╗  ##
## ██╔══██╗██╔════╝██╔══██╗██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔════╝██╔════╝  ##
## ██████╔╝█████╗  ██████╔╝██████╔╝██║   ██║██║     █████╗  ███████╗███████╗  ##
## ██╔══██╗██╔══╝  ██╔═══╝ ██╔══██╗██║   ██║██║     ██╔══╝  ╚════██║╚════██║  ##
## ██║  ██║███████╗██║     ██║  ██║╚██████╔╝╚██████╗███████╗███████║███████║  ##
## ╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝  ╚═╝ ╚═════╝  ╚═════╝╚══════╝╚══════╝╚══════╝  ##
##                                                                            ##
################################################################################
################################################################################

def reprocess_observations_for(cluster:Cluster) -> None:
    """Reprocesses the cluster foregrounds in order to expand each of their 
    initial sizes so that that each image size is now large enough to accomadate 
    the fields of all observations present. This is done so image coordinates on
    one observation match image coordinates on another observation. 

    CIAO commands run
    -----------------
    chandra_repro

    Parameters
    ----------
    cluster : ClusterObj
        The ClusterPyXT ClusterObject you are merging observations for.

    Returns
    -------
    bool
        True or false depending on successful merging of observations
    """
    logging.info(f"Reprocessing {cluster.obsids} for {cluster.name}")
    
    # chandra_repro is the first CIAO command we need to run to get all of the
    # observations the same size. It is easiest to call this command with the 
    # clusters working directory as the current working directory. First we 
    # will get the current directory to change back at the end, then we will
    # change the directory to the current clusters.
    start_directory = os.getcwd()
    os.chdir(cluster.data_directory)

    with mp.Pool(CPXTConfig().num_cores) as pool:
        results = list(tqdm(
            pool.imap(reprocess_obs, cluster.observations),
            desc="Reprocessing observations",
            total=len(cluster.observations),
            unit="observation"
            ))    
    # Change back to the starting directory
    os.chdir(start_directory)
    logging.info(f"Reprocessing {cluster.name} complete.")
    

def reprocess_obs(observation: Observation) -> str:
    """
    Takes a `ClusterPyXT` `Observation` object and calls `chanrda_repro` on it
    to reprocess the observation's datafiles. This checks each observation and 
    caluclates the field of view necessary to fit each observation onto it. Then
    it reprocesses each observation with this new field of view so each 
    observations pixel coordinates are on the same scale / point of reference.
    
    Parameters
    ----------
    observation : Observation
        The particular observation we are reprocessing.
    
    Returns
    -------
    str
        The results from the `CIAO` `chandra_repro` call.
    """
    keyword_args = {
        'indir': [f"{observation.id}"], 
        'set_ardlib': False, # can mess things up when processing multiple obs.
        'clobber': True,
        }
    try:
        output = ciao.run_command(rt.chandra_repro, **keyword_args)# type:ignore
    except:
        raise

    return str(output)


def make_level_2_event_files_for(cluster: Cluster) -> None:
    """
    From an early stage in processing, it can be beneficial to use a level 2 
    events file that is filtered to each specific CCD. This function takes a 
    `cpxt` cluster, and iteratively generates those files per observation. All
    of this is wrapped in a `tqdm` progress bar to help update the user as to 
    the status of the operation.
    
    Parameters
    ----------
    cluster : Cluster
        The cluster we want to make level2 event files for. 
    
    Returns
    -------
    None
    """
    logger.info(f"Making level 2 event files for {cluster.name}")
    for observation in tqdm(cluster.observations,
                            total=len(cluster.obsids),
                            desc="Generating Level 2 Event Files",
                            leave=True):
        # the leave=True flag helps keep this bar present while each of the
        # `tqdm` progress bars generated from the below function call
        #  will disappear when completed as their progress bar has the argument
        # leave=False.
        generate_level_2_event_files_for_obs(observation)
    

def generate_level_2_event_files_for_obs(observation: Observation) -> None:
    """
    For each observation, for various reasons it is beneficial to filter out
    each particular CCDs level 2 event data into its own file. This function 
    takes an observation, and generates a new level 2 event file filtered for
    each CCD. It is saved in that observations analysis directory as:
        /path/to/cluster/observation/analysis/acis_ccd#.fits
    With the corresponding CCD number replacing the #. 

    The function uses the `multiprocessing.Pool` library with a context manager
    to ensure the pool is closed properly and prevent issues. The actually call,
    `pool.imap` is wrapped in a `tqdm` progress bar to help convey the job 
    status to the user. 
    
    Parameters
    ----------
    observation : Observation
        The particular `Observation` object we want to split into its subsequent
        CCD images. 
    
    Returns
    -------
    None
    """
    logger.info(f"Generating level 2 event files for {observation}")

    try:
        observation.analysis_dir.mkdir(exist_ok=True)
    except:
        raise

    with mp.Pool(CPXTConfig().num_cores) as pool:
        results = list(tqdm(
            pool.imap(generate_level_2_event_files_for_ccd, observation.ccds),
            desc=f"Gen. {observation.id} evt2 for each CCD",
            total=len(observation.ccds),
            unit="CCD",
            leave=False 
            # The `leave=False` argument makes this bar disappear when complete.
            # We want this functionality as this function is called from within
            # a loop multiple times. We don't want to litter the screen with 
            # this functions progress bars. 
            ))
    
    
def generate_level_2_event_files_for_ccd(ccd: CCD) -> None:
    """
    This function is called in conjunction with the above function, 
    `generate_level_2_event_files_for_obs`. When that function is called on an
    observation, it calls this function for each of that observations CCDs. This
    function is the one that generates the new level 2 events file for that 
    particular CCD. It is saved in the following location:
        /path/to/cluster/observation/analysis/acis_ccd#.fits
    With the corresponding CCD number replacing the #.
    
    Files Generated
    ---------------
    ./cluster/obsid/analysis/acic_ccd{ccd}.fits
    
    Parameters
    ----------
    ccd : CCD
        The CCD we want to filter out and generate its own level 2 event file.
    
    Returns
    -------
    None
    """
    keyword_args = {
        "infile": ccd.observation.reprocessed_evt2_file,
        "outfile": ccd.acis_ccd_evt2_file,
        "clobber": True    
        }
    try:
        result = ciao.run_command(rt.dmcopy, **keyword_args)        #type:ignore
    except:
        raise

################################################################################
################################################################################    
##                                                                            ##
##   ███████╗ ██████╗ ██████╗ ███████╗                                        ## 
##   ██╔════╝██╔═══██╗██╔══██╗██╔════╝                                        ##
##   █████╗  ██║   ██║██████╔╝█████╗█████╗                                    ##
##   ██╔══╝  ██║   ██║██╔══██╗██╔══╝╚════╝                                    ##
##   ██║     ╚██████╔╝██║  ██║███████╗                                        ##
##   ╚═╝      ╚═════╝ ╚═╝  ╚═╝╚══════╝                                        ##
##                                                                            ##
##            ██████╗ ██████╗  ██████╗ ██╗   ██╗███╗   ██╗██████╗ ███████╗    ##
##           ██╔════╝ ██╔══██╗██╔═══██╗██║   ██║████╗  ██║██╔══██╗██╔════╝    ##
##           ██║  ███╗██████╔╝██║   ██║██║   ██║██╔██╗ ██║██║  ██║███████╗    ##
##           ██║   ██║██╔══██╗██║   ██║██║   ██║██║╚██╗██║██║  ██║╚════██║    ##
##           ╚██████╔╝██║  ██║╚██████╔╝╚██████╔╝██║ ╚████║██████╔╝███████║    ##
##            ╚═════╝ ╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝╚═════╝ ╚══════╝    ##
##                                                                            ##
################################################################################
################################################################################

def merge_foregrounds_for(cluster: Cluster) -> None:
    """
    This function merges the foregrounds using the CCD's filtered out for each
    particular observation. The filtering process is handled by the 
    `Observation` class and subsequently the `Cluster` class. We take a list of
    the filtered CCD counts images (level 2 event files) and merge them to a 
    combined mutlt-observation image (assuming obsids > 1). This is done using
    the CIAO commands listed below. 
    
    CIAO commands
    -------------
    dmmerge - https://cxc.harvard.edu/ciao/ahelp/dmmerge.html
        Merges the level 2 event files across all observations filtered to the 
        CCDs we care about.

    Files Generated
    ---------------
    ./cluster/

    Parameters
    ----------
    cluster : Cluster
        The cluster we are merging the individual CCD level 2 event files for. 

    Returns
    -------
    None    
    """

    # Need to write the ACIS CCD filenames to a file so we can call dmmerge once
    # just passing that file instead of building a long string we are likely to
    # mess up and clutter log files with. 
    long_string_of_filenames = "\n".join(cluster.acis_ccd_evt2_file_list)
    io.write_contents_to_file(long_string_of_filenames, 
                              str(cluster.acis_ccd_evt2_filenames))

    keyword_args = {
        "infile": f"@{cluster.acis_ccd_evt2_filenames}",
        "outfile": cluster.merged_acis_events,
        "clobber": True
        }
    ciao.run_command(rt.dmmerge, **keyword_args)                    #type:ignore

################################################################################
################################################################################    
##                                                                            ##
##  ██████╗  █████╗  ██████╗██╗  ██╗                                          ## 
##  ██╔══██╗██╔══██╗██╔════╝██║ ██╔╝                                          ##
##  ██████╔╝███████║██║     █████╔╝█████╗                                     ##
##  ██╔══██╗██╔══██║██║     ██╔═██╗╚════╝                                     ##
##  ██████╔╝██║  ██║╚██████╗██║  ██╗                                          ##
##  ╚═════╝ ╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝                                          ##
##                                                                            ##
##            ██████╗ ██████╗  ██████╗ ██╗   ██╗███╗   ██╗██████╗ ███████╗    ##
##           ██╔════╝ ██╔══██╗██╔═══██╗██║   ██║████╗  ██║██╔══██╗██╔════╝    ##
##           ██║  ███╗██████╔╝██║   ██║██║   ██║██╔██╗ ██║██║  ██║███████╗    ##
##           ██║   ██║██╔══██╗██║   ██║██║   ██║██║╚██╗██║██║  ██║╚════██║    ##
##           ╚██████╔╝██║  ██║╚██████╔╝╚██████╔╝██║ ╚████║██████╔╝███████║    ##
##            ╚═════╝ ╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝╚═════╝ ╚══════╝    ##
##                                                                            ##
################################################################################
################################################################################

def reprocess_backgrounds_for(cluster: Cluster) -> None:
    """
    The goal of this function is to iteratively go through the cluster's 
    observations, find backgrounds for each ACIS CCD image created above, and
    reprocess them to get them on the same field of view. This is done in
    parallel per observation per ccd.
    
    Parameters
    ----------
    cluster : Cluster
        The `Cluster` object we are working on.
    
    Returns
    -------
    None
    """
    logger.info(f"Finding and reprocessing backgrounds for {cluster.name}")

    # This is just a fancy call to `reprocess_obs_backgrounds_for` for each
    # of the clusters observations. The `with` context manager takes care of
    # opening and closing the multiprocessing pool. `tqdm` wraps the call in
    # order to display a progress bar to the user. 
    with mp.Pool(CPXTConfig().num_cores) as pool:
        results = list(tqdm(
            pool.imap(reprocess_obs_backgrounds_for, cluster.observations),
            total=len(cluster.observations),
            desc=f"Finding/reprocessing backgrounds for {cluster.name}",
            unit="Observation"
            ))


def reprocess_obs_backgrounds_for(observation: Observation) -> None:
    """
    This function goes through the observation and for each CCD present it calls
    `find_ccd_background_for(ccd)` and then reprocesses that background file to
    get it into the same field of view as we did with all of the CCD images. 
    This is done so that the x,y pixel coordinates are the same across all 
    images used in this project. 
    
    Parameters
    ----------
    observation : Observation
        The `cpxt` `Observation` class we are working on. 
    
    Returns
    -------
    None
    """
    logger.debug(f"Finding/reprocessing background for {observation.id}")
    for ccd in observation.ccds:
        caldb_background = find_ccd_background_for(ccd)
        logger.debug(f"Found background for {ccd} at {caldb_background}")
        
        # Copy the background to the observations analysis directory
        io.copy(caldb_background, ccd.background)

        keyword_args = {
            "infile": ccd.acis_ccd_evt2_file,
            "keyword": "GAINFILE",
            "echo":True
        }
        ccd_gain_file = \
            ciao.run_command(rt.dmkeypar, **keyword_args)         # type: ignore
        
        # Now get the background `GAINFILE`
        keyword_args['infile'] = ccd.background
        ccd_back_gain = \
            ciao.run_command(rt.dmkeypar, **keyword_args)         # type: ignore
        
        if not io.dates_and_versions_match(ccd_gain_file, ccd_back_gain):
            gain_path = caldb_background.split('/')[:-2] # CALDB path
            gain_path.append('det_gain')
            gain_path.append(ccd_gain_file)
            gain_file = '/'.join(gain_path)
        
            keyword_args = {
                "infile": caldb_background,
                "outfile": str(ccd.background), # convert from PathLib to str
                # "acaofffile": None,
                # "stop": None,
                "doevtgrade": False,
                "apply_cti": True,
                "apply_tgain": False,
                "calculate_pi": True,
                "pix_adj": 'EDSER',
                "gainfile": gain_file,
                "clobber": True,
                "eventdef": "{s:ccd_id,s:node_id,i:expno,s:chip,s:tdet,f:det,"
                            "f:sky,s:phas,l:pha,l:pha_ro,f:energy,l:pi,"
                            "s:fltgrade,s:grade,x:status}"
            }
            ciao.run_command(rt.acis_process_events, **keyword_args)#type:ignore



def find_ccd_background_for(ccd: CCD) -> str:
    """
    For each ACIS CCD image generated above, we are going to search for the 
    appropriate background from the CALDB. This is done by iterating through 
    every `CCD` from every `Observation` of the `Cluster`. If
    
    CIAO Calls
    ----------
    acis_bkgrnd_lookup
        Checks the CALDB for the background file for the given input file.
        This is the reason the full CALDB is required during installation.

    Parameters
    ----------
    ccd : CCD
        The CCD we are querying the CALDB for a background file. 
    
    Returns
    -------
    str
        A `str` containing the location of the background file for this
        observaiton.
    """
    keyword_args = {
        "infile": ccd.acis_ccd_evt2_file
    }
    background = \
        ciao.run_command(rt.acis_bkgrnd_lookup, **keyword_args)     #type:ignore
    
    # This command will return a list of background files for every ccd on the 
    
    split_background = background.split('\n')
    if len(split_background) > 1:
        files = [Path(filename) for filename in split_background]
        for file in files:
            if file.name.startswith(f"acis{ccd.id}"):
                return str(file)
    return background

################################################################################
################################################################################    
##                                                                            ##
##  ███╗   ███╗███████╗██████╗  ██████╗ ███████╗     ██████╗ ██████╗ ███████╗ ##
##  ████╗ ████║██╔════╝██╔══██╗██╔════╝ ██╔════╝    ██╔═══██╗██╔══██╗██╔════╝ ##
##  ██╔████╔██║█████╗  ██████╔╝██║  ███╗█████╗      ██║   ██║██████╔╝███████╗ ##
##  ██║╚██╔╝██║██╔══╝  ██╔══██╗██║   ██║██╔══╝      ██║   ██║██╔══██╗╚════██║ ##
##  ██║ ╚═╝ ██║███████╗██║  ██║╚██████╔╝███████╗    ╚██████╔╝██████╔╝███████║ ##
##  ╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝     ╚═════╝ ╚═════╝ ╚══════╝ ##
################################################################################
################################################################################

def merge_observations_for(cluster: Cluster) -> None:
    """
    Merge each of the observations into a combined flux image. 
    
    CIAO Calls
    ----------
    fluximage
        Called when the cluster is made single observation.
    
    merge_obs
        Called when the cluster is made up of multiple observations.

    Parameters
    ----------
    cluster : Cluster
        The current `Cluster` object we are working on.
    
    Returns
    -------
        None
    """
    
    logger.info(f'Merging observations from {cluster}')
    print(f'Merging observations from {cluster}')
    # create the file & directory for list of observations to merge
    merged_file = Path(f'{cluster.merged_directory}/merged_obs.lis')
    merged_file.parent.mkdir(parents=True, exist_ok=True)

    obs_to_merge = [str(obs.reprocessed_evt2_file_ccd_filter) \
                                                for obs in cluster.observations]
                                        
    io.write_contents_to_file('\n'.join(obs_to_merge), str(merged_file))

    kwargs = {
            'infile': f"@{str(merged_file)}", # @/path/to/merged_obs.lis
            'outroot': Path(f'{cluster.data_directory}/{cluster.name}'),
            'binsize': 4,
            'clobber': True,
            'cleanup': True
        }

    if len(cluster.observations) == 1:
        command = rt.fluximage  # type: ignore
    else:
        command = rt.merge_obs    # type: ignore
        kwargs['parallel'] = True
        kwargs['nproc'] = CPXTConfig().num_cores
    
    ciao.run_command(command, **kwargs)
    logger.info(f"Done merging observations for {cluster.name}")
    print("Done merging observations.")


def print_stage_1_prep() -> None:
    """This function prints the preparation message for stage 1.

    Parameters
    ----------
    None.
    
    Returns
    -------
    None
    """
    
    message = io.load_message(Stage.one, MessageType.preparation)
    print(message)
    logging.info(message)


def print_stage_1_comp(cluster: Cluster) -> None:
    """This function prints the completion message for stage 1.

    Parameters
    ----------
    cluster : Cluster
        The ClusterPyXT Cluster object of the cluster you are working on.
    
    Returns
    -------
    None
    """
    message_args = {
        "sb_map_filename": cluster.xray_sb_map_filename,
    }

    message = io.load_message(Stage.one, MessageType.completion, **message_args)
    print(message)
    logging.info(message)