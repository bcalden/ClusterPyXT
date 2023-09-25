"""
File: cpxt/stages/stage_five.py
Author: Brian Alden
Date: 20 Sep 2023
Description: This file contains the code to run stage five of `ClusterPyXT`. 
             
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
##          ███████╗████████╗ █████╗  ██████╗ ███████╗    ███████╗            ##
##          ██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝    ██╔════╝            ##
##          ███████╗   ██║   ███████║██║  ███╗█████╗      ███████╗            ##
##          ╚════██║   ██║   ██╔══██║██║   ██║██╔══╝      ╚════██║            ##
##          ███████║   ██║   ██║  ██║╚██████╔╝███████╗    ███████║            ##
##          ╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝    ╚══════╝            ##
##                                                                            ##
################################################################################
################################################################################

def calculate_acb_bins_for(cluster:Cluster) -> bool:
    """Function description
    
    Parameters
    ----------
    cluster : Cluster
        The current cluster you are working on

    
    Returns
    -------
    bool
        True or false depending on successful completion
    """
    return False

def generate_scale_map_for(cluster:Cluster) -> bool:
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

def calculate_exposure_corrections_for(cluster:Cluster) -> bool:
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