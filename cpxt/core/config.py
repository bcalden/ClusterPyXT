"""
File: cpxt/core/config.py
Author: Brian Alden
Date: 28 Mar 2022
Description: This file contains the code for managing the `ClusterPyXT` 
             configuration. 
"""

# Internal imports

# External imports
from configparser import ConfigParser
from pathlib import Path
from enum import Enum
import logging
import sys
import os

logger = logging.getLogger(__name__)

class CFGKey(str, Enum):
    """
    This is the enumeration coding the keys we use for our configuration file. 
    We use the enumeration as to allow for easier referencing of configuration
    dictionary keys. 
    """
    cpxt = 'cpxt'
    data_directory = 'data_directory'
    astro_query = 'astro_query'
    num_downloads = 'num_downloads'
    num_cores = 'num_cores'
    log_level = 'log_level'


class CPXTConfig:
    """
    The CPXTConfig class, and ultimately its instatiation, captures and queries
    information about the runtime environment and where to store data. 
    
    Attributes
    ----------
    config_filename : str
        The `ClusterPyXT` configuration filename.
    log_filename : str
        The `ClusterPyXT` log file.
    ciao_version : str
        The CIAO version reported in the shell environment.
    CALDB_directory : str
        The CALDB directory reported in the shell environment.
    
    Methods
    -------
    load_cpxt_config : dict
        CPXT Configuration dictionary
    """

    config_directory = Path(f"{Path.home()}/.config/cpxt/")
    config_directory.mkdir(parents=True, exist_ok=True)
    config_filename = f"{config_directory}/cpxt.cfg"
    log_filename = f"{config_directory}/cpxt.log"
    ciao_version = os.getenv("CONDA_DEFAULT_ENV")
    CALDB_directory = os.getenv("CALDB")

    def __init__(self):
        """
        CPXTConfig intialization function. Trys to load the configuration
        from, otherwise it assumes this is the first run and tries to create
        the configuration file.          
        
        Parameters
        ----------
        None

        Returns
        -------
        CPXTConfig
            An instatiation of the CPXTConfig class.
        """

        try:
            self._config = self.load_cpxt_config()
        except KeyError:
            logger.critical("Was unable to find the ClusterPyXT config file "
                           f"at {self.config_filename}. Try restarting."
                           )
            raise
        
        self.data_directory = self._config[CFGKey.data_directory.value]

        try:
            self.num_cores = int(self._config[CFGKey.num_cores.value])
        except KeyError:
            self.num_cores = 1
        
        try:
            self.astro_query = bool(self._config[CFGKey.astro_query.value])
        except KeyError:
            self.astro_query = False

        try:
            self.num_downloads = int(self._config[CFGKey.num_downloads.value])
        except KeyError:
            self.num_downloads= 1


    def load_cpxt_config(self) -> dict:
        """
        Trys to load the `ClusterPyXT` configuration file. 
        
        Parameters
        ----------
        N/A
        
        Returns
        -------
        dict
            A dictionary containing the `ClusterPyXT` configuration parameters.
        """
        
        config_parser = ConfigParser()
        config_parser.read(CPXTConfig.config_filename)
        
        try:
            config_dict = dict(config_parser[CFGKey.cpxt.value])  
            logger.debug(f"Loaded config file, {CPXTConfig.config_filename}")
        except KeyError:
            logger.error(f"Unable to process config file. Recreating.")
            config_dict = {}
            raise

        return config_dict


def write_cpxt_config(data_directory: str="", 
                      astro_query: bool=True,
                      num_downloads: int=1,
                      num_cores: int=1,
                      log_level: int=1) -> None:
    """
    This function writes the configuration settings to the `cpxt` config
    file. The file can be viewed and edited in a text editor. 
    
    Parameters
    ----------
    data_directory : str
        The location you want `ClusterPyXT` to store its data. 
    astro_query : bool
        Should `cpxt` reach out to NED and HEASARC to automatically get cluster
        properties. 
        NED - https://ned.ipac.caltech.edu/
        NASA HEASARC - https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl
    num_downloads : int
        The number of simultaneous downloads we allow at once (1-5). 
    num_cores : int
        The number of cores/threads to use for parallal processing.
    
    Returns
    -------
    None
    """

    try:
        os.makedirs(CPXTConfig.config_directory)
    except FileExistsError:
        pass

    config_parser= ConfigParser()
    config = {
        CFGKey.data_directory.value: data_directory,
        CFGKey.astro_query.value: astro_query,
        CFGKey.num_cores.value: num_cores,
        CFGKey.num_downloads.value: num_downloads,
        CFGKey.log_level.value: log_level
    }
    config_parser[CFGKey.cpxt.value] = config                      #type: ignore
    with open(CPXTConfig.config_filename, 'w') as f:
        config_parser.write(f)


def reset_config():
    """
    Removes the configuration file if it exists. 
    
    Files Deleted
    -------------
    ~/.config/cpxt/cpxt.cfg

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """
    try:
        os.remove(CPXTConfig.config_filename)
        logger.info(f"Config file, {CPXTConfig.config_filename} reset.")
    except FileNotFoundError:
        logger.error(
            f"No configuration file found at {CPXTConfig.config_filename}."
            )
        return


def reset_log():
    """
    This function deletes the log file if it exists. 
    
    Files Generated
    ---------------
    ~/.config/cpxt/cpxt.log

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    try:
        os.remove(CPXTConfig.log_filename)
    except:
        logger.error(f"No log file found at {CPXTConfig.log_filename}.")
