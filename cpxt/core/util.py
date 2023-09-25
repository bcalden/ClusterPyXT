"""
File: cpxt/core/util.py
Author: Brian Alden
Date: 20 Sep 2023
Description: This file contains utility code for ClusterPyXT. 
"""

# Internal imports

# External imports
from configparser import ConfigParser
from pathlib import Path
from enum import Enum
import subprocess
import logging
import sys
import os

logger = logging.getLogger(__name__)


def run_application(application:str , arguments:list(str)) -> None:
    """
    This function runs the application provided with the arguments provided.

    Parameters:
    -----------
    application (str): The path of the application to run
    
    arguments (list(str)): A list of strings containing the arguments to be
                           supplied to the application.
    
    Returns
    -------
    None
    """
      
    # Check to see if the application exists
    if not Path(application).exists():
        # The application does not exist
        logger.critical(f"Application {application} does not exist.")
        sys.exit(-1)
    else:
        # The application exists
        logger.info(f"Running {application} with arguments {arguments}")

        command = [application] + arguments
        subprocess.run(command)