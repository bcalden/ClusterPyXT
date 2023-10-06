"""
File: cpxt/stages/stage_three.py
Author: Brian Alden
Date: 20 Sep 2023
Description: This file contains the code to run stage three of `ClusterPyXT`. 
             
"""

# Internal imports
from cpxt.stages.stage_four import print_stage_4_prep
from cpxt.chandra.observation import Observation
from cpxt.core.stages import Stage, MessageType
from cpxt.chandra.cluster import Cluster
from cpxt.core.config import CPXTConfig
from cpxt.chandra.ccd import CCD
from cpxt.chandra import ciao
from cpxt.core.util import run_application
from cpxt.core import io

# External imports
from ciao_contrib import runtool as rt
import multiprocessing as mp
from pathlib import Path
from tqdm import tqdm
import logging


logger = logging.getLogger(__name__)

################################################################################
################################################################################
##                                                                            ##
##          ███████╗████████╗ █████╗  ██████╗ ███████╗    ██████╗             ##
##          ██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝    ╚════██╗            ##
##          ███████╗   ██║   ███████║██║  ███╗█████╗       █████╔╝            ##
##          ╚════██║   ██║   ██╔══██║██║   ██║██╔══╝       ╚═══██╗            ##
##          ███████║   ██║   ██║  ██║╚██████╔╝███████╗    ██████╔╝            ##
##          ╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝    ╚═════╝             ##
##                                                                            ##
################################################################################
################################################################################


def run_on(cluster:Cluster) -> None:
    """Runs stage three of the pipeline on the provided cluster. First 

    Parameters
    ----------
    cluster : ClusterObj
        The ClusterPyXT ClusterObject of the particular cluster 
        you want to run stage three on.

    Returns
    -------
    None
    """
    if has_required_files(cluster):
        extract_response_files_in_parallel_for(cluster)
        generate_rmf_files_in_parallel_for(cluster)
        generate_arf_files_in_parallel_for(cluster)
        print_stage_3_comp(cluster)
        print_stage_4_prep(cluster)
    else: 
        print_stage_3_prep(cluster)
        make_region_files(cluster)


def has_required_files(cluster: Cluster) -> bool:
    """Checks if the required files for stage three are present in the 
    cluster's main output directory. The required file is a region file named
    acis(i/s)_region_0.reg. This file contains a small region (approximately 
    20-40 arcseconds) covering a small piece of each of the ACIS CCDs. This 
    region is used to characterize the gain across all the different chips.

    Parameters
    ----------
    cluster : ClusterObj
        The ClusterPyXT ClusterObject of the particular cluster 
        you want to run stage three on.

    Returns
    -------
    bool
        True if the required files are present, False if not.
    """
    for observation in cluster.observations:
        if not io.file_exists(observation.gain_region_file):
            return False
    return True

def make_region_files(cluster: Cluster) -> None:
    for observation in tqdm(cluster.observations, desc="Creating region files"):
        region_file = observation.gain_region_file
        if (not io.file_exists(region_file)) or \
             (io.file_size(region_file) == 0):
            # print("Region file {} does not exist.".format(region_file))
            # print("When DS9 opens, draw a small circle that covers a piece of each ACIS-I chip (~20 pixels) and save it as:\n" \
            #       "{}".format(region_file))
            # print("Opening SAO DS9")
            io.write_contents_to_file("", region_file, False)
            ds9_arguments = [f"ds9 -regions system physical -regions shape " \
                             f"circle -regions format ciao -zoom 0.5 -bin " \
                             f"factor 4 -scale log " \
                             f"{observation.clean_data_filename}"]
            
            run_application("", ds9_arguments, shell=True)

def extract_response_files_in_parallel_for(cluster: Cluster) -> None:
    """Extract response files in parallel for the given cluster.
    
    Parameters
    ----------
    cluster : Cluster
        The current cluster you are working on
    
    Returns
    -------
    None
    """
    logger.info(f"Extracting response files for {cluster}.")
    observations = cluster.observations
    with mp.Pool(processes=len(observations)) as pool:
        results = list(tqdm(pool.imap(extract_response_files_for, observations),
                            total=len(observations),
                            unit="observation",
                            desc="Extracting response files"
                            ))
    
    if not all(results):
        logger.error(f"Failed to extract response files for {cluster}.")
        raise RuntimeError(f"Failed to extract response files for {cluster}.")

def extract_response_files_for(observation: Observation) -> bool:
    """Extract response files for the given cluster.
    
    Parameters
    ----------
    cluster : Cluster
        The current cluster you are working on
    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    logger.info(f"Extracting response files for {observation}.")
    try:
        obs_analysis_dir = observation.analysis_dir
        global_response_dir = observation.global_response_dir
        global_response_dir.mkdir(parents=True, exist_ok=True)
        
        bad_pixel_file = str(observation.reprocessed_bad_pixel_filename)
        clean = observation.clean_data_filename
        
        ardlib_args = {
            'badpixfile': bad_pixel_file
        }
        ciao.run_command(rt.acis_set_ardlib, **ardlib_args)                # type: ignore
        
        mask_file = str(observation.mask_file)
        make_pcad_lis(observation)
        
        specextract_args = {
            'infile': f"{clean}[sky=region({observation.gain_region_file})]",
            'outroot': f"{observation.gain_region_arf_file}"[:-4],
            'weight': True,
            'correctpsf': False,
            'asp': f"@{obs_analysis_dir}/pcad_asol1.lis",
            'combine': False,
            'mskfile': mask_file,
            'bkgfile': "",
            'bkgresp': False,
            'badpixfile': bad_pixel_file,
            'grouptype': "NUM_CTS",
            'binspec': 1,
            'clobber': True
        }
        ciao.run_command(rt.specextract, **specextract_args)      # type: ignore
        
        return True
    except:
        raise
        return False

def generate_rmf_files_in_parallel_for(cluster: Cluster):
    """Generate RMF files in parallel for the given cluster.
    
    Parameters
    ----------
    cluster : Cluster
        The current cluster you are working on
    
    Returns
    -------
    None
    """
    logger.info(f"Generating RMF files for {cluster}.")
    observations = cluster.observations
    with mp.Pool(processes=len(observations)) as pool:
        results = list(tqdm(pool.imap(generate_rmf_files_for, observations),
                            total=len(observations),
                            unit="observation",
                            desc="Generating RMF files"
                            ))
    
    if not all(results):
        logger.error(f"Failed to generate RMF files in parallel for {cluster}.")
    for observation in tqdm(cluster.observations, desc="Generating RMF files"):
        if not generate_rmf_files_for(observation):
            logger.error(f"Failed to generate RMF files for {cluster}.")
            raise RuntimeError(f"Failed to generate RMF files for {cluster}.")

def generate_rmf_files_for(observation: Observation) -> bool:
    """Generate RMF files for the given cluster.
    
    Parameters
    ----------
    cluster : Cluster
        The current cluster you are working on
    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    logger.info(f"Generating RMF files for {observation}.")
    try:
        global_response_dir = observation.global_response_dir
        global_response_dir.mkdir(parents=True, exist_ok=True)
        back = observation.background_filename
        
        dmextract_args = {
            'infile': f"{back}[sky=region(" \
                      f"{observation.gain_region_file})][bin pi]",
            'outfile': str(observation.background_gain_region_pi_filename),
            'clobber': True
        }
        ciao.run_command(rt.dmextract, **dmextract_args)          # type: ignore
        
        return True
    except:
        return False
    
def generate_arf_files_in_parallel_for(cluster: Cluster):
    """Generate ARF files in parallel for the given cluster.
    
    Parameters
    ----------
    cluster : Cluster
        The current cluster you are working on
    
    Returns
    -------
    None
    """
    logger.info(f"Generating ARF files for {cluster}.")
    observations = cluster.observations
    with mp.Pool(processes=len(observations)) as pool:
        results = list(tqdm(pool.imap(generate_arf_files_for, observations),
                            total=len(observations),
                            unit="observation",
                            desc="Generating ARF files"
                            ))
    
    if not all(results):
        logger.error(f"Failed to generate ARF files in parallel for {cluster}.")
    for observation in tqdm(cluster.observations, desc="Generating ARF files"):
        if not generate_arf_files_for(observation):
            logger.error(f"Failed to generate ARF files for {cluster}.")
            raise RuntimeError(f"Failed to generate ARF files for {cluster}.")

def generate_arf_files_for(observation: Observation) -> bool:
    """Generate ARF files for the given cluster.
    
    Parameters
    ----------
    observation : Observation
        The current observation you are working on
    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    logger.info(f"Generating ARF files for {observation}.")
    try:
        global_response_dir = observation.global_response_dir
        
        dmhedit_args = {
            'infile': observation.gain_region_pi_filename,
            'filelist': "",
            'operation': "add",
            'key': "BACKFILE",
            'value': observation.background_gain_region_pi_filename
        }
        ciao.run_command(rt.dmhedit, **dmhedit_args)              # type: ignore
        
        io.copy(observation.gain_region_arf_file,
                observation.aux_response_filename)

        io.copy(observation.gain_region_rmf_file,
                observation.redist_matrix_filename)
                
        return True
    except:
        raise
        return False


def make_pcad_lis(observation: Observation):
    search_str = "{}/*asol1.fits".format(observation.reprocessed_dir)
    pcad_files = [str(s) for s in io.get_filenames_matching(search_str)]
    pcad_list_string = "\n".join(pcad_files)
    pcad_filename = "{}/pcad_asol1.lis".format(observation.analysis_dir)

    io.write_contents_to_file(pcad_list_string, pcad_filename, binary=False)

    return pcad_filename


def print_stage_3_prep(cluster: Cluster) -> None:
    """This function prints the instructions for preparing the cluster for stage
    three of the pipeline. This includes the creation of the region file that
    will be used to extract the gain for each of the CCDs.

    Parameters
    ----------
    cluster : Cluster
        The ClusterPyXT Cluster object of the cluster you are working on
    
    Returns
    -------
    None
    """
    message_args = {
        "cluster_name": cluster.name,
        "region_file": cluster.observations[0].gain_region_file,
        "clean": cluster.observations[0].clean_data_filename
    }

    message = io.load_message(Stage.three, 
                              MessageType.preparation, 
                              **message_args)
    print(message)
    logging.info(message)

def print_stage_3_comp(cluster: Cluster):
    """This function prints the instructions for completing stage three of the
    pipeline. This includes the creation of the region file that will be used to
    extract the gain for each of the CCDs.

    Parameters
    ----------
    cluster : Cluster
        The ClusterPyXT Cluster object of the cluster you are working on
    
    Returns
    -------
    None
    """
    
    message = io.load_message(Stage.three, 
                              MessageType.completion)
    print(message)
    logging.info(message)
