"""
File: cpxt/core/arguments.py 
Author: Brian Alden
Date: 28 Mar 2022
Description: This file deals with arguments passed on the command line to `cpxt`. 
"""
# Internal imports
from cpxt.core.config import reset_config, reset_log

# External imports
from argparse import ArgumentParser
import sys
import logging
import os

def get_cli_args():
    """
    This function uses the `ArgumentParser` class from `argparse`. It allows
    for a simple command line interface (CLI) into ClusterPyXT. 
    
    Parameters
    ----------
    N/A
    
    Returns
    -------
    dict
        A python dictionary containing the arguments passed using the CLI
    """
    
    logging.debug("Getting the arguments passed at the command line.")

    # load the CLI description from a file.
    with open(f"cpxt/core/text_files/cpxt_cli_desc.txt", "r") as f:
        description = f.read()
    
    program_name = "ClusterPyXT"

    parser = ArgumentParser(prog=program_name, description=description)

    parser.add_argument("--initialize_cluster", "-i",
        dest="init_cluster",
        action="store_true", 
        help="Initializes the cluster and adds it to ClusterPyXT.",
        default=False)
    
    parser.add_argument("-hydrogen-column-density", "-nH",
        dest="nH",
        action="store", 
        help="The hydrogen column density for the cluster in units of /cm^2. "\
            "Input data using sci. notation. E.g. 5.2e19 for 5.2x10^19 /cm^2"\
            "Values can be obtained at the following url:"\
            "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl",
        default=None)
    
    parser.add_argument("--redshift", "-z",
        dest="redshift", 
        action="store",
        help="The redshift of the cluster. Values can be found at: "\
             "https://ned.ipac.caltech.edu/", 
        default=None)
    
    parser.add_argument("--cluster_name", "-n", 
        dest="name",
        action="store",
        help="The name of the cluster you are going to process."\
            "This value will also be used as a filename prefix.",
        default=None)
    
    parser.add_argument("--obsids", "-o",
        dest="obsids", 
        action="store", 
        nargs="+",
        help="The observation ids from CXO Chaser you want to use. If "\
            "multiple OBSIDs are to be used, enter them one after another"\
            "with a space in between. E.g. -o 3233 13458 13459 15578 15581."\
            " OBSIDs can be queried using Chaser at: "\
            "https://cda.harvard.edu/chaser/",
        default=None)
    
    parser.add_argument("--abundance", "-a",
        dest="abundance", 
        help="The metallicity or solar abundance of the galaxy cluster. If you "\
            "a value for this can be found in the literature for the "\
            "cluster of concern. A fair estimate for most clusters is 0.2-0.3.",
        action="store", 
        default=None)
    
    parser.add_argument("--reset-config", "-r",
        dest="reset_config",
        help="Resets the `cpxt` configuration file. The next run of `cpxt`"\
            " will prompt you for the data storage directory again.",
        action="store_true",
        default=False)

    parser.add_argument("--reset-log", "-l",
        dest="reset_log",
        help="Removes the current log file and beings logging to a new one.",
        action="store_true",
        default=False)

    arguments = parser.parse_args()

    return arguments

def process_cli_args() -> None:
    """
    This function first calls `get_cli_args` (see above) to get the arguments
    passed at the command line. Then it processes them to see which step/stage
    of `ClusterPyXT` to run. 
    
    Parameters
    ----------
    N/A
    
    Returns
    -------
    N/A
    """
    logging.debug("Processing commandline arguments")
    arguments = get_cli_args()

    if arguments.reset_config:
        reset_config()

    if arguments.reset_log:
        reset_log()
    
    if arguments.reset_config or arguments.reset_log:
        # The files have been reset at this point. We need to query the user
        # to find a new directory to store data in and create a new config file
        # exit the program so the user can rerun and create these.
        return
    