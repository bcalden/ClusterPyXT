"""
File: cpxt/chandra/ciao.py
Author: Brian Alden
Date: 29 Mar 2022
Description: This file contains all the code used to interact with CIAO 
"""

# Internal imports
# None

# External imports
from ciao_contrib.cda.data import download_chandra_obsids
from ciao_contrib.runtool import CIAOToolParFile
from ciao_contrib import runtool as rt
import multiprocessing as mp
from pathlib import Path
from tqdm import tqdm
import logging
import os

logger = logging.getLogger(__name__)

def download_obsid(obsid: int):
    """
    Uses the `download_chandra_obsids` function from CIAO to download the passed
    obsid. This function downloads to the current working directory. 
    
    Parameters
    ----------
    obsid : int
        The Chandra OBSID we are asking CIAO to download.  
    
    Returns
    -------
    str
        The result output of the CIAO command `download_chandra_obsids`.
    """
    logger.info(f"Downloading data for observation id {obsid}")
    return download_chandra_obsids([obsid])


def download_observations(obsids: list, num_streams: int=1):
    """
    Function description goes here.
    
    Parameters
    ----------
    obsids : list
        A list of integer obsids to download from the Chandra X-ray Observatory.
    num_streams : int
        The number of simultaneous download streams to use. An `int` from 1-5
    
    Returns
    -------
    list
        `download_chandra_obsids` returns a string containing download info.
        This function returns a list of the returned strings for each obsid.
    """

    # Use the multiprocessing library's Pool function to handle concurrency
    with mp.Pool(num_streams) as pool:
        # Using tqdm to display a progress bar. To use `tqdm` and `mp.Pool`
        # we use the `pool.imap` method. This differs from the `.map` in that
        # it is a 'lazy' version of map. We use the wrap it all in the `list`
        # type casting as it forces python to wait for `imap` to finish. `tqdm`
        # only works with `imap` not `map` so this work around is necessary.
        #
        # https://stackoverflow.com/a/26521507 for map vs imap differences
        
        results = list(tqdm(pool.imap(download_obsid, obsids),
                            total=len(obsids),
                            unit="observation",
                            desc="Downloading observations"
                            ))
    
    return results


def get_ciao_version(full: bool=False) -> str:
    """
    This function queries the running shell environment for the `CALDB` 
    variable. It then checks the `VERSION` file in the root of the environment
    folder. 
    
    Parameters
    ----------
    full : bool
        Return the full version string (including date) or just version number.
    
    Returns
    -------
    str
        A string containing the CIAO version loaded in the current shell env.
    """

    try:
        CALDB_location = Path(str(os.getenv('CALDB')))
    except:
        raise
    
    environment_directory = CALDB_location.parent

    version_file = f"{environment_directory}/VERSION"
    with open(version_file, 'r') as f:
        full_version = f.read()

    if full:
        return full_version

    # The version string looks something like this:
    # CIAO 4.13 Wednesday, December  2, 2020
    # As full==False if we are here, we only want the part before the date

    split_version = full_version.split(" ")
    version = " ".join(split_version[:2])

    return version


def run_command(command: CIAOToolParFile, **kwargs: dict) -> str:
    """
    This is the main function that should be called to run any CIAO command that
    is handled by the imported runtool. This is essentially every command EXCEPT
    the `download_chandra_obsids` command which is handled seperately in this 
    file. This function wraps the command and first clears any temporary files 
    from a previous run of the same command. Next, it runs the command with the
    supplied keyword arguments (**kwargs).
    
    Parameters
    ----------
    command : ciao_contrib.runtool.CIAOToolParFile
        This is the actual CIAO command you want to run. E.g. dmcopy (rt.dmcopy)
    kwargs : dictionary
        The keyword arguments to be passed to the command. can be passed as a
        double asterix'ed dictionary or typed out keyword=value, keyword=value 
        in the function call. I.e. both are valid:

    keyword_args = {"infile": "/some/file", "outfile": "/other/file.ext"}
    ciao_run_command(rt.dmcopy, **keyword_args)

    or

    ciao_run_command(rt.dmcopy, infile="/some/file", outfile="/other/file.ext")

    And they will run the same. 

    Note: If you don't specify full paths you need to change the directory 
          before running this function

    Returns
    -------
    str
        A string containing the results CIAO passes back to the console. 
    """
    logging.info(f"Running {command._toolname} with {kwargs}")
    try:
        command.punlearn()
        result = command(**kwargs)
        return str(result)
    except:
        logging.error(f"{command._toolname} failed with {kwargs} as arguments.")
        raise