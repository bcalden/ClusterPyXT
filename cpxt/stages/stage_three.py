"""
File: cpxt/stages/stage_three.py
Author: Brian Alden
Date: 20 Sep 2023
Description: This file contains the code to run stage three of `ClusterPyXT`. 
             
"""

# Internal imports
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
        extract_response_files_for(cluster)
        generate_rmf_files_for(cluster)
        generate_arf_files_for(cluster)
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
    return False

def make_region_files(cluster: Cluster) -> None:
    for observation in tqdm(cluster.observations, desc="Creating region files"):
        region_file = observation.region_file
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
        print('Creating global response file.')

def extract_response_files_for(cluster:Cluster) -> bool:
    """Function description
    
    Parameters
    ----------
    cluster : cl.ClusterObj
        The current cluster you are working on

    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    return False

def generate_rmf_files_for(cluster:Cluster) -> bool:
    """Function description
    
    Parameters
    ----------
    cluster : cl.ClusterObj
        The current cluster you are working on
    
    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    return False

def generate_arf_files_for(cluster:Cluster) -> bool:
    """Function description
    
    Parameters
    ----------
    cluster : cl.ClusterObj
        The current cluster you are working on

    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    return False

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
        "region_file": cluster.observations[0].region_file,
        "clean": cluster.observations[0].clean_data_filename
    }

    message = io.load_message(Stage.three, 
                              MessageType.preparation, 
                              **message_args)
    print(message)
    logging.info(message)