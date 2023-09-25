"""
File: cpxt/stages/stage_two.py
Author: Brian Alden
Date: 20 Sep 2023
Description: This file contains the code to run stage two of `ClusterPyXT`. 
             Specifically, it removes point sources and flares from the images.
             To do so, it requires a sources.reg and exclude.reg files. The 
             sources.reg file contains the ds9 regions of the point sources to
             be removed. The exclude.reg file contains the ds9 region(s) of the 
             peak of cluster emission to be excluded from the deflaring process.
"""

# Internal imports
from ast import keyword
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
##          ███████╗████████╗ █████╗  ██████╗ ███████╗    ██████╗             ##
##          ██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝    ╚════██╗            ##
##          ███████╗   ██║   ███████║██║  ███╗█████╗       █████╔╝            ##
##          ╚════██║   ██║   ██╔══██║██║   ██║██╔══╝      ██╔═══╝             ##
##          ███████║   ██║   ██║  ██║╚██████╔╝███████╗    ███████╗            ##
##          ╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝    ╚══════╝            ##
##                                                                            ##
################################################################################
################################################################################

def run_on(cluster:Cluster) -> None:
    """Runs stage two of the pipeline on the provided cluster. First 

    Parameters
    ----------
    cluster : ClusterObj
        The ClusterPyXT ClusterObject of the particular cluster 
        you want to run stage one on.

    Returns
    -------
    None
    """
    if check_for_required_files(cluster):
        remove_point_sources_for(cluster)
        generate_light_curves(cluster)
        filter_high_energy_flares_for(cluster)
    # remove_sources_in_parallel(cluster,args)
    # generate_light_curves(cluster, args)
    # make_nosrc_xray_sb(cluster)
    # lightcurves_with_exclusion(cluster, args)

def check_for_required_files(cluster:Cluster) -> bool:
    """This function checks to make sure that the cluster has the required files
    for stage two to run. Specifically, it checks for the sources.reg and the
    exclude.reg files. If either of these files are missing, it will print an
    error message and call the `print_stage_2_prep` function.
    

    Parameters
    ----------
        cluster (Cluster): The CPXT Cluster object of the cluster you are
                           working on

    Returns
    -------
        bool: True if the cluster has the required files, False otherwise
    """
    if io.file_exists(cluster.sources_file) and io.file_exists(cluster.exclude_file):
        return True
    else:
        io.print_red("Error: Missing {sources} and/or {exclude}".format(
            sources=cluster.sources_file,
            exclude=cluster.exclude_file
        ))
        print_stage_2_prep(cluster)
    
    return False


def remove_point_sources_for(cluster:Cluster) -> bool:
    """This function attempts to remove the point sources identified in the
    sources.reg file from the cluster's observations. It does so by calling
    the `remove_sources_from_obs` function on each observation in the cluster.
    If the `remove_sources_from_obs` function fails in parallel, it will try 
    again in serial.
    
    Parameters
    ----------
    cluster : cl.ClusterObj
        The current cluster you are working on

    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    try:
        with mp.Pool(CPXTConfig().num_cores) as pool:
            _ = list(tqdm(pool.imap(remove_sources_from_obs, 
                                    cluster.observations),
                          total=len(cluster.observations),
                          desc='Removing point sources',
                          unit='observations'))
    except Exception as e:  # Update to catch specific exceptions if possible.
        print(f"Error in removing sources in parallel.")
        print(f"CPU count:{CPXTConfig().num_cores}\nException: {e}\n")
        print("Trying again with single core.")
        
        # Remove sources in serial
        for observation in tqdm(cluster.observations,
                                total=len(cluster.observations),
                                desc='Removing point sources',
                                unit='observation'):
            remove_sources_from_obs(observation)
    return True

def remove_sources_from_obs(observation):
    logging.info(f"Removing sources from {observation.id}")

    # We are excluding sources in both foreground and background files
    # and copying the background to observation.back
    operations = [
        (observation.data_filename, observation.acis_nosrc_filename),
        (observation.back_filename, observation.background_nosrc_filename)
    ]

    sources_file = observation.sources_filename

    for original_file, no_src_file in operations:
        infile = f"{original_file}[exclude sky=region({sources_file})]"
        
        logging.debug(f"infile: {infile}")
        logging.debug(f"outfile: {no_src_file}")

        # Using dmcopy with the exclude sky=region() syntax to remove sources
        try:
            ciao.run_command(rt.dmcopy, 
                             infile=infile, 
                             outfile=no_src_file, 
                             clobber=True, 
                             verbose=0)
        except Exception as e:
            logging.error(f"Failed to run dmcopy with infile: {infile}," \
                          f"outfile: {no_src_file}. Error: {e}")
            raise
    
    # Copy the background file to the observation.back file
    logging.info(f"Copying background to {observation.back}")
    io.copy(observation.background_nosrc_filename, observation.back)


def generate_light_curves(cluster:Cluster) -> bool:
    """
    __description of function__
    
    Parameters
    ----------
    param1 : type
        description of param1
    
    Returns
    -------
    type : description of return
    """
    return False

def filter_high_energy_flares_for(cluster:Cluster) -> bool:
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

def print_stage_2_prep(cluster: Cluster) -> None:
    """This function prints the instructions for preparing the cluster for stage
    two of the pipeline. Specifically, it tells the user to open the surface
    brightness map and create regions around the point sources they want to
    exclude from the data analysis. It also tells the user to create a region
    file containing any regions they want excluded from the deflaring process.

    Parameters
    ----------
    cluster : Cluster
        The ClusterPyXT Cluster object of the cluster you are working on
    
    Returns
    -------
    None
    """
    message_args = {
        "sources_file": cluster.sources_file,
        "exclude_file": cluster.exclude_file,
        "cluster_name": cluster.name
    }

    message = io.load_message(Stage.two, 
                              MessageType.preparation, 
                              **message_args)
    print(message)
    logging.info(message)