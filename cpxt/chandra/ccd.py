"""
File: cpxt/chandra/ccd.py 
Author: Brian Alden
Date: 1 Apr 2022
Description: This file contains the CCD dataclass. It is mainly used to keep 
             track of the CCD and observation/cluster with only a single 
             object.
"""

# Internal imports
from cpxt.chandra import observation as obs

# External imports
from pathlib import Path

class CCD:
    """
    
    Attributes
    ----------
    <var> : <type>
        Description
    
    Methods
    -------
    <method> : <return type>
    """
    def __init__(self, id : int, observation : obs.Observation) -> None:
        self.id = id
        self.observation = observation

################################################################################
################################################################################
#    
#                     ███████╗██╗██╗     ███████╗███████╗
#                     ██╔════╝██║██║     ██╔════╝██╔════╝
#                     █████╗  ██║██║     █████╗  ███████╗
#                     ██╔══╝  ██║██║     ██╔══╝  ╚════██║
#                     ██║     ██║███████╗███████╗███████║
#                     ╚═╝     ╚═╝╚══════╝╚══════╝╚══════╝
#                                                                       
################################################################################
################################################################################

    @property
    def reprocessed_evt2_file(self) -> Path:
        """
        This function modifies the observations level 2 event file with a string
        indicating to CIAO we only care about that specific CCD. When this is 
        passed as an infile, it filters out all other CCD's from processing.
        
        Parameters
        ----------
        ccd : int
            The specific ACIS CCD we care about. Should be single digit integer.
        
        Returns
        -------
        Path
            A `pathlib.Path` object modifying the level 2 event file for the 
            specific CCD we care about. 
        """
        return \
            Path(f"{self.observation.reprocessed_evt2_file}[ccd={self.id:d}]")
    
    @property
    def acis_ccd_evt2_file(self) -> Path:
        """
        This function generates a `Path` objection pointed at a fits file 
        filtered to this specific CCDs data. There is no need for any string
        modifier on this file to filter out any data, only the data from this
        CCD is contained within this file. This differs to how the 
        `reprocessed_evt2_file_for` method above generates its `Path` object. 
        
        Parameters
        ----------
        ccd : int
            The specific CCD we want the filtered 
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the reprocessed filtered level
            2 events file (evt2) for this specific CCD and observation.
        """
        return Path(f"{self.observation.analysis_dir}/acis_ccd{self.id}.fits")

    @property
    def background(self) -> Path:
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
        return Path(f"{self.observation.analysis_dir}/back_ccd{self.id}.fits")
    
################################################################################
################################################################################
##                                                                            ##
##                   ███╗   ███╗ █████╗  ██████╗ ██╗ ██████╗                  ##
##                   ████╗ ████║██╔══██╗██╔════╝ ██║██╔════╝                  ##
##                   ██╔████╔██║███████║██║  ███╗██║██║                       ##
##                   ██║╚██╔╝██║██╔══██║██║   ██║██║██║                       ##
##   ███████╗███████╗██║ ╚═╝ ██║██║  ██║╚██████╔╝██║╚██████╗███████╗███████╗  ##
##   ╚══════╝╚══════╝╚═╝     ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚═╝ ╚═════╝╚══════╝╚══════╝  ##
##                                                                            ##
################################################################################
################################################################################

    def __repr__(self):
        """
        This function returns a printable version of the `CCD` object. It 
        overrides pythons `magic` method (signified by double underscore).
        We really just want to dispaly the obsid and ccd as whatever is 
        displaying this likely is already indicating the cluster as well. If not
        the cluster can be identified using the observation id where only the
        CCD id number would not convey such information. 
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        return f"{self.observation.id}-CCD: {self.id}"