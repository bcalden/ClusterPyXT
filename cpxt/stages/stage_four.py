"""
File: cpxt/stages/stage_four.py
Author: Brian Alden
Date: 20 Sep 2023
Description: This file contains the code to run stage four of `ClusterPyXT`. 
             
"""

# Internal imports
from cpxt.chandra.observation import Observation
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
import logging


logger = logging.getLogger(__name__)



################################################################################
################################################################################
##                                                                            ##
##          ███████╗████████╗ █████╗  ██████╗ ███████╗    ██╗  ██╗            ##
##          ██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝    ██║  ██║            ##
##          ███████╗   ██║   ███████║██║  ███╗█████╗      ███████║            ##
##          ╚════██║   ██║   ██╔══██║██║   ██║██╔══╝      ╚════██║            ##
##          ███████║   ██║   ██║  ██║╚██████╔╝███████╗         ██║            ##
##          ╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝         ╚═╝            ##
##                                                                            ##
################################################################################
################################################################################

def crop_using_master_crop(cluster:Cluster):
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


def filter_data_to(cluster:Cluster, low_energy=0.7, high_energy=8.0):
    """Filters all used event data to the passed energy levels.
    Default is 0.7-8.0 keV. Alternative values can be used (2-5 is also common)
    These values are chosen based on the Quantum Efficiency of the ACIS CCDs.
    More info can be found in chapter 6 of the proposers guide. 
    https://cxc.cfa.harvard.edu/proposer/POG/html/chap6.html#tth_sEc6.5

    
    Parameters
    ----------
    cluster : cl.ClusterObj
        The current `cluster` you are working on
    
    low_energy : float, default=0.7
        Filter out all events below `low_energy` (in keV).

    high_energy : float, default=8.0
        Filter out all events above `high_energy` (in keV).

    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    return False
