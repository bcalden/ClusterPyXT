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
    if has_required_files(cluster):
        remove_point_sources_for(cluster)
        make_nosrc_xray_sb_image(cluster)
        generate_light_curves(cluster)
        filter_high_energy_flares_for(cluster)
        light_curves_with_exclusion(cluster) # rename function
    # remove_sources_in_parallel(cluster,args)
    # generate_light_curves(cluster, args)
    # make_nosrc_xray_sb(cluster)
    # lightcurves_with_exclusion(cluster, args)
    else: 
        print_stage_2_prep(cluster)


def has_required_files(cluster:Cluster) -> bool:
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
        logging.info(f"Found {cluster.sources_file} and {cluster.exclude_file}")
        return True
    else:
        io.print_red(f"Error: Missing {cluster.sources_file} and/or " \
                     f"{cluster.exclude_file}")
        logging.info(f"Missing {cluster.sources_file} and/or " \
                     f"{cluster.exclude_file}.")
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
        keyword_args = {
            "infile": infile,
            "outfile": no_src_file,
            "clobber": True,
            "verbose": 0
        }
        try:
            ciao.run_command(rt.dmcopy, **keyword_args)           # type: ignore
        except Exception as e:
            logging.error(f"Failed to run dmcopy with infile: {infile}," \
                          f"outfile: {no_src_file}. Error: {e}")
            raise
    
    # Copy the background file to the observation.back file
    logging.info(f"Copying background to {observation.back}")
    io.copy(observation.background_nosrc_filename, observation.back)


def make_nosrc_xray_sb_image(cluster:Cluster) -> bool:
    """This function creates a surface brightness map of the cluster with the 
    point sources identified in sources.reg removed. It achieves this by calling
    the `dmcopy` command from CIAO and excluding the region file `sources.reg`
    
    Parameters
    ----------
    cluster : cl.ClusterObj
        The current cluster you are working on.

    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    excluding = f"exclude sky=region({cluster.sources_file})"
    keyword_args = {
        "infile": f"{cluster.xray_sb_map_filename}[{excluding}]",
        "outfile": cluster.xray_sb_nosrc_map_filename,
        "clobber": True,
    }
    ciao.run_command(rt.dmcopy, **keyword_args)                   # type: ignore
    return True

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
    for observation in tqdm(cluster.observations, 
                            total=len(cluster.observations), 
                            desc='Generating light curves', 
                            unit='observation'):
        generate_light_curve(observation)
    return True

def generate_light_curve(observation: Observation) -> bool:
    """
    Generates a light curve for a given observation by performing a series of 
    data processing steps. These steps include filtering out high energy 
    background flares, binning data, extracting light curves, cleaning light 
    curves, and filtering event lists using GTI (Good Time Interval) information 
    derived from high energy flares. The results of each step are saved as files 
    in the observation's analysis directory.
    
    Parameters
    ----------
    observation : cpxt.Observation
        The observation object containing information and file paths related to 
        the observation. 
    
    Returns
    -------
    bool
        True if the light curve is successfully generated, False otherwise.

    Notes
    -----
    This function relies on CIAO tools for processing the observation data and 
    generates several intermediate and final output files in the observation's 
    analysis directory. Ensure that the required CIAO tools are available in the 
    environment, and the observation object is correctly initialized with the 
    necessary file paths.
    """
    obsid_analysis_dir = observation.analysis_dir
    data = observation.acis_nosrc_filename

    # Filter to just high energy events from 9-12 keV
    keyword_args = {
        "infile": f"{observation.acis_nosrc_filename}[energy=9000:12000]",
        "outfile": observation.high_energy_data_file,
        "clobber": True
    }
    ciao.run_command(rt.dmcopy, **keyword_args)                   # type: ignore

    # Bin the event data
    keyword_args = {
        "infile": f"{observation.high_energy_data_file}[bin sky=8]",
        "outfile": observation.binned_high_energy_data_file,
        "clobber": True
    }
    ciao.run_command(rt.dmcopy, **keyword_args)                   # type: ignore

    # Get start and stop times 
    keyword_args = {
        "infile": observation.high_energy_data_file,
        "keyword": "TSTART",
        "echo": True
    }
    tstart = ciao.run_command(rt.dmkeypar, **keyword_args)        # type: ignore

    keyword_args = {
        "infile": observation.high_energy_data_file,
        "keyword": "TSTOP",
        "echo": True
    }
    tstop = ciao.run_command(rt.dmkeypar, **keyword_args)         # type: ignore
    
    # Extract lightcurve
    bin_string = f"bin time={tstart}:{tstop}:259.28"
    keyword_args = {
        "infile": f"{observation.high_energy_data_file}[{bin_string}]",
        "outfile": observation.high_energy_light_curve_file,
        "opt": "ltc1",
        "clobber": True
    }
    ciao.run_command(rt.dmextract, **keyword_args)                # type: ignore
    
    # Clean the lightcurve
    keyword_args = {
        "infile": observation.high_energy_light_curve_file,
        "outfile": observation.good_time_interval_light_curve_file,
        "method": "clean",
        "save": observation.good_time_interval_light_curve_file.with_suffix("")
    }
    ciao.run_command(rt.deflare, **keyword_args)                  # type: ignore
    
    return True


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
    for observation in tqdm(cluster.observations, 
                            total=len(cluster.observations), 
                            desc='Filtering High Energy Flares', 
                            unit='observation'):
        filter_high_energy_flares_from(observation)
    return True

def filter_high_energy_flares_from(observation:Observation) -> bool:
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
    
    # Filter the event list using GTI info from high energy flares
    gti_file = observation.good_time_interval_light_curve_file
    keyword_args = {
        "infile": f"{observation.acis_nosrc_filename}[@{gti_file}]",
        "outfile": observation.acis_nosrc_high_energy_filename,
        "clobber": True
    }
    ciao.run_command(rt.dmcopy, **keyword_args)                   # type: ignore
    
    # Final binning
    keyword_args = {
        "infile": f"{observation.acis_nosrc_high_energy_filename}[bin sky=8]",
        "outfile": observation.binned_full_energy_nosrc_data_file,
        "clobber": True,
        "verbose": 0
    }
    ciao.run_command(rt.dmcopy, **keyword_args)                   # type: ignore
    
    return True


def light_curves_with_exclusion(cluster:Cluster) -> bool:
    for observation in tqdm(cluster.observations, 
                            desc='Finishing light curves', 
                            unit='observation', 
                            total=len(cluster.observations)):
        light_curve_with_exclusion_for(observation)
    return True

def light_curve_with_exclusion_for(observation:Observation) -> bool:
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