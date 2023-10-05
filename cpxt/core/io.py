"""
File: cpxt/core/io.py
Author: Brian Alden
Date: 31 Mar 2022
Description: This file encodes the core input/output operations of ClusterPyXT.
             Input/output to both screen and storage are handled below.
"""

# Internal imports
from cpxt.core.stages import Stage, MessageType

# External imports
from os.path import exists
from pathlib import Path
from glob import glob
import shutil as sh
import logging
import json
import os
import re

logger = logging.getLogger(__name__)

class Colors:
    BLACK = "\u001b[30m"
    RED = "\u001b[31m"
    GREEN = "\u001b[32m"
    YELLOW = "\u001b[33m"
    BLUE = "\u001b[34m"
    MAGENTA = "\u001b[35m"
    CYAN = "\u001b[36m"
    WHITE = "\u001b[37m"
    RESET = "\033[0;0m"
    BOLD = "\033[;1m"

def color_string(string, color):
    return f"{color}{string}{Colors.RESET}"

def print_red(string):
    print(red_text(string))

def red_text(string):
    return color_string(string, Colors.RED)

def green_text(string):
    return color_string(string, Colors.GREEN)

def get_filenames_matching(pattern:str) -> list:
    """
    This function uses the `glob` library to find filenames matching the pattern
    provided. Simple wild cards such as an asterix, `*`, can be used to match 
    file names you only know pieces of. Types of patterns that can be used can 
    be found here https://docs.python.org/3/library/glob.html .

    Parameters
    ----------
    pattern : str
        The file pattern we are trying to match. For example:
            ..ching('/path/to/cluster/obs/analysis/acis_ccd[0-3].fits')
        would return a list of the 4 'acis_ccd#.fits' files for that particular
        observation. 
    
    Returns
    -------
    list
        A list of `Path` objects pointing to the matching files. 
    """
    try:
        return [Path(file) for file in glob(pattern)]
    except:
        logger.error(f"No files matching {pattern} was found.")
        raise    # To be implemented further as exceptions arise.


def get_filename_matching(pattern:str) -> Path:
    """
    This function calls the `get_filenames_matching` function and returns the 
    first result. Useful when only one file is expected and you want to deal
    with that `Path` object not figure out if you got a `list` returned and then
    having to figure out which index to use. This alleviates that concern. 
    
    Parameters
    ----------
    pattern : strr
        The pattern string matching the filename we are looking for.
    Returns
    -------
    Path
        A `pathlib.Path` object pointing to the first matching file from `glob`.
    """
    try:
        return get_filenames_matching(pattern)[0]
    except:
        logger.error(f"No file matching {pattern} was found.")
        raise    # To be implemented further as exceptions arise.

def get_temp_filename(filename: str | Path) -> Path:
    """
    This function takes a filename and returns a `Path` object pointing to a 
    temporary file in the same directory as the provided filename. 
    
    Parameters
    ----------
    filename : str
        The filename we are trying to get a temporary file for.
    
    Returns
    -------
    Path
        A `pathlib.Path` object pointing to a temporary file in the same 
        directory as the provided filename.
    """
    try:
        return Path(filename).with_suffix('.tmp')
    except:
        logger.error(f"Unable to get temporary filename for {filename}.")
        raise    # To be implemented further as exceptions arise.

def file_size(filename: str | Path) -> int:
    return os.path.getsize(filename)

def write_contents_to_file(contents, filename: str | Path, binary: bool=False) -> None:
    """
    This function writes `contents` to `filename`. The binary flag encodes 
    whether or not we are writing a binary data file or a text file. We use a
    context meanager to open the file for writing so that python handles closing
    the file when we are done with it. 
    
    Parameters
    ----------
    contents : Basically anything 
        The data/text we want to write to a file. Can be a variety of types. 
    filename : str
        The filename you want to write the contents to.
    binary : bool
        Flags if `contents` should be written in binary format or text format. 

    Returns
    -------
    None
    """    
    try:
        file_attributes = 'wb' if binary else 'w'    # binary file or text file 
        with open(filename, file_attributes) as f:
            f.write(contents)
    except FileExistsError:
        logger.error(f"{filename} already exists and I can't over write!")
        raise
    except FileNotFoundError:
        logger.error(f"{filename} not found!")
        raise
    except:
        raise # As more exceptions arise, implement them here. 


def read_contents_of_file(filename: str) -> str:
    """
    This function reads the contents of the file located at `filename`. If the
    file is not found, log it and raise the exception.
    
    Parameters
    ----------
    filename : str
        The filename you are looking to get the contents of. Also accepts a 
        `pathlib.Path` object. 
    
    Returns
    -------
    str
        A string containing the contents of the file located at `filename`.
    """
    try:
        with open(filename) as f:
            data = f.read()
        return str(data)
    except FileNotFoundError:
        logger.error(f"Unable to find {filename}. Did you complete the "
                     "previous stage of cpxt?"
            )
        raise
    except:
        logger.error("ClusterPyXT encountered an error. Please see the "
            "exception reported below. Feel free to try again and/or report "
            "the issue on GitHub."
            )
        raise # To be implemented as more exceptions arise.


def copy(source: str | Path, destination: str | Path, replace: bool=True):
    if not replace and exists(destination):
        logger.error(f"{destination} exists and replace set to False.")
        return
    try:
        sh.copy(source, destination)
    except:
        raise # To be implemented as exceptions arise.

def dates_and_versions_match(acis_filename: str, 
                             background_filename: str) -> bool:
    """
    Function description goes here.
    
    Parameters
    ----------
    <var> : <type>
        <description>
    
    Returns
    -------
    <type>
        <description>
    """                
    acis = {
        'date': get_date_from_filename(acis_filename),
        'version': get_version_from_filename(acis_filename)
        }

    background = {
        'date': get_version_from_filename(background_filename),
        'version': get_version_from_filename(background_filename)
        }

    return (acis['date'] == background['date']) and \
           (acis['version'] == background['version']) 


def get_date_from_filename(filename):
    """
    This function searches the given filename for a date based on how Chandra
    X-ray Observatory names their files. 
    
    Example filename: acisD2000-01-29gain_ctiN0006.fits
    
    Parameters
    ----------
    filename : str
        The given filename we are trying to parse a date from. 
    
    Returns
    -------
    str
        date
    """

    regex = r"(19|20)\d\d[- /.](0[1-9]|1[012])[- /.](0[1-9]|[12][0-9]|3[01])"

    return get_regex_result_from_string(regex, filename)


def get_regex_result_from_string(regex, string):
    """
    Function description goes here.
    
    Parameters
    ----------
    <var> : <type>
        <description>
    
    Returns
    -------
    <type>
        <description>
    """
    match = re.search(regex, string)
    if match:
        return match.group()                   

def get_version_from_filename(filename):
    """
    Function description goes here.
    
    Parameters
    ----------
    <var> : <type>
        <description>
    
    Returns
    -------
    <type>
        <description>
    """
    # acisD2000-01-29gain_ctiN0006.fits

    regex = r"N\d{0,}."

    return get_regex_result_from_string(regex, filename)

def str_to_bool(string: str) -> bool:
    """
    Takes an input string value such as "true" or "false" and returns the 
    according boolean value. Implemented here instead of using 
    `distutils.util.strtobool` as distutils is deprected in Python 3.11 and up.
    
    Parameters
    ----------
    string : str
        The string we are evaluating
    
    Returns
    -------
    bool
        The appropriate `bool` representation of the string.
    """
    if string.lower() in ['true', 't', 'yes', 'y', '1']:
        return True
    return False

def load_message(stage: Stage, message_type: MessageType, **kwargs) \
                                                                        -> str:
    """
    Loads the appropriate message from the pipeline_messages.json file and
    formats it with the provided kwargs.
    
    Parameters
    ----------
    stage : Stage
        The stage we are currently on.
    message_type : MessageType
        The type of message we are trying to load. (preparation or completion)
    
    **kwargs : dict
        The key word arguments we are using to format the message string.
    
    Returns
    -------
    str : The message string for the appropriate stage and message type formated
          with the provided kwargs.
    """

    # Get the current file's directory
    module_dir = os.path.dirname(os.path.realpath(__file__))

    # Construct the path to the assets directory
    assets_dir = os.path.abspath(os.path.join(module_dir, '..', 'assets'))

    # Construct the path to myfile.txt
    message_file = os.path.join(assets_dir, 'stage_messages.json')

    with open(message_file, 'r') as file:
        data = json.load(file)
    
    stage_key = f"stage_{stage.value}"
    if stage_key not in data['stages']:
        raise ValueError(f"Invalid stage: {stage}")

    if message_type.value not in data['stages'][stage_key]:
        raise ValueError(f"Invalid message type: {message_type.value}")
    
    # Format the template string with the actual values provided in kwargs
    return data['stages'][stage_key][message_type.value].format(**kwargs)


def file_exists(filename: str | Path) -> bool:
    """
    Checks if the given file exists. Returns True if it does, False if it does
    not. 
    
    Parameters
    ----------
    filename : str | Path
        The filename we are checking for existance. 
    
    Returns
    -------
    bool
        True if the file exists, False if it does not.
    """
    return Path(filename).exists()

