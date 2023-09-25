"""
File: cpxt/chandra/observation.py
Author: Brian Alden
Date: 31 Mar 2022
Description: This contains the `Observation` class and the ACIS enumeration. 
             The `Observation` class makes a python object pointing to all 
             aspects of the observation we care about in `ClusterPyXT`. 

             The `ACIS` enumeration codifies which type of observation we have,
             an ACIS-I or an ACIS-S. The CCD layout for each is shown below. 
             More information about the ACIS sensors can be found at:
      
          https://cxc.cfa.harvard.edu/proposer/POG/html/chap6.html#tth_chAp6

ACIS: Advanced CCD Imaging Spectrometer - Sensor layout
                               
                             ------- -------
                            |       |       |
                            | CCD:0 | CCD:1 |
                            |       |       |
                             ------- -------          <---  ACIS-I
                            |       |       |
                            | CCD:2 | CCD:3 |           
                            |       |       |           
                             ------- -------
                                                       
             ------- ------- ------- ------- ------- -------
            |       |       |       |       |       |       |  
            | CCD:4 | CCD:5 | CCD:6 | CCD:7 | CCD:8 | CCD:9 |  <--- ACIS-S
            |       |       |       |       |       |       |
             ------- ------- ------- ------- ------- -------
"""

# Internal imports
from cpxt.chandra import cluster as cl
from cpxt.chandra import ciao
from cpxt.core import io

# External imports
from ciao_contrib import runtool as rt
from pathlib import Path
from enum import IntEnum
import numpy as np
import logging

logger = logging.getLogger(__name__)

class ACIS(IntEnum):
    I = 0
    S = 1


class Observation:
    """
    The `Observation` class codifies various properties of a Chandra observation
    we need to access throughout `ClusterPyXT`. Things such as specific files
    and directories can be referenced in the observation class without having to
    type out the string each time. Further, specifics of the observation, such 
    as which CCDs were on during the observation, can be accessed through this
    observation class. 
    
    Attributes
    ----------
    obsid : int
        The Chandra obsid for this observation. 
    cluster : Cluster
        The `cpxt` `Cluster` object the observation belongs to.
    
    Methods
    -------
    directory : Path
        Returns a `pathlib.Path` object for the observations main directory.
    analysis_dir : Path
        Returns a `pathlib.Path` object for the observations analysis directory.
    secondar_dir : Path
        Returns a `pathlib.Path` object for the observations secondary directory
    reprocessed_dir : Path
        Returns a `pathlib.Path` object for directory containing the reprocessed
        files. 
    
    """

    def __init__(self, obsid: int, cluster: cl.Cluster):
        self.id = obsid
        self.cluster = cluster

################################################################################
################################################################################
#
#           ███████╗ ██████╗ ██╗     ██████╗ ███████╗██████╗ ███████╗
#           ██╔════╝██╔═══██╗██║     ██╔══██╗██╔════╝██╔══██╗██╔════╝
#           █████╗  ██║   ██║██║     ██║  ██║█████╗  ██████╔╝███████╗
#           ██╔══╝  ██║   ██║██║     ██║  ██║██╔══╝  ██╔══██╗╚════██║
#           ██║     ╚██████╔╝███████╗██████╔╝███████╗██║  ██║███████║
#           ╚═╝      ╚═════╝ ╚══════╝╚═════╝ ╚══════╝╚═╝  ╚═╝╚══════╝
#                                                                       
################################################################################
################################################################################
    
    @property
    def directory(self) -> Path:
        """
        This function encodes the observations main data directory. It should 
        be:

            /path/to/clusterpyxt-data/<clustername>/<obsid>/
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the observations main directory.
        """
        return Path(f"{self.cluster.data_directory}/{self.id}")

    @property
    def analysis_dir(self) -> Path:
        """
        This property returns a `pathlib.Path` object pointing to the 
        observations analysis directory. This directory is used and created by
        `ClusterPyXT`. It contains the processed products we use to do the 
        spectral fitting for each observation.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the observation's analysis dir.
        """
        return Path(f"{self.directory}/analysis/")

    @property
    def primary_dir(self) -> Path:
        """
        The clusters primary directory. This directory is created by CIAO and 
        
        Parameters
        ----------
        <var> : <type>
            <description>
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the observation's
        """
        return Path(f"{self.directory}/primary/")

    @property
    def secondary_dir(self) -> Path:
        """
        Function description goes here.
        
        Parameters
        ----------
        <var> : <type>
            <description>
        
        Returns
        -------
         Path
            A `pathlib.Path` object pointing to the observation's
        """
        return Path(f"{self.directory}/secondary/")

    @property
    def reprocessed_dir(self) -> Path:
        """
        Function description goes here.
        
        Parameters
        ----------
        <var> : <type>
            <description>
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the observation's 
        """
        return Path(f"{self.directory}/repro")
                                                                       
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
    def level_1_event_file(self) -> Path:
        """
        Returns a `Path` object pointed at the level 1 event file for the obsid.
        There is only one file that starts with `acis` and ends with `evt1.fits`
        in the `secondary` directory of the observation. This function finds 
        that file, and returns a `pathlib.Path` object pointed at it. 
        """
        return io.get_filename_matching(f"{self.secondary_dir}/acis*evt1.fits")
    
    @property
    def reprocessed_evt2_file(self) -> Path:
        """
        Returns a `Path` object pointed at the reprocessed level 2 event file
        for this particular observation. It is the only file in the reprocessing
        directory that starts with 'acis' and ends with 'evt2.fits'. 
        """
        return io.get_filename_matching(
            f"{self.reprocessed_dir}/acis*repro_evt2.fits"
        )

    @property
    def reprocessed_evt2_file_ccd_filter(self) -> str:
        """
        Returns a string pointing to the reprocessed level 2 event file
        for this particular observation filtered for ccds [0:3]. It is the only 
        file in the reprocessing directory that starts with 'acis' and ends with 
        'evt2.fits'. This would need to change for an ACIS-S observation
        """
        return f"{str(self.reprocessed_evt2_file)}[ccd_id=0:3]"

    @property
    def acis_ccd_fits_files(self) -> list[Path]:
        """
        This function returns a list of this observations reprocessed level 2
        event files per CCD. First we check which type of observation it is, 
        then we get a list of files matching the CCDs we expect (ACIS-I or S).
        
        Parameters
        ----------
        None
        
        Returns
        -------
        list
            A list of `Path` objects pointing to each of the `acis_ccd#.fits` 
            files. 
        """
        if self.acis_type == ACIS.I:
            ccds = "0-3"
            # The ACIS-I chip uses CCDs 0,1,2, and 3. We only want those.
        elif self.acis_type == ACIS.S:
            ccds = "4-9"
            # Similarly, the ACIS-S chip uses CCDs 4,5,6,7,8,9. 
        else:
            # No CCDs found. This is an error so let us handle it as such.
            logger.error(f"Unable to determine CCDs used for {self.id}")
            raise ValueError(f"Unable to determine CCDs used for {self.id}")
        return io.get_filenames_matching(
            f"{self.analysis_dir}/acis_ccd[{ccds}].fits"
        )

    @property
    def high_energy_data_file(self) -> Path:
        """
        This function returns a `Path` object pointing to the high energy data
        file for this observation. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the high energy data file for
            this observation.
        """
        return Path(f"{self.analysis_dir}/acisI_hiE.fits")
    
    @property
    def binned_high_energy_data_file(self) -> Path:
        """
        This function returns a `Path` object pointing to the binned high energy 
        data file for this observation. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the binned high energy data file 
            for this observation.
        """
        return Path(f"{self.analysis_dir}/img_acisI_hiE.fits")
    
    @property
    def binned_full_energy_nosrc_data_file(self) -> Path:
        """
        This function returns a `Path` object pointing to the binned full energy 
        data file for this observation. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the binned full energy data file 
            for this observation.
        """
        return Path(f"{self.analysis_dir}/img_acisI_nosrc_fullE.fits")
        
    @property
    def acis_nosrc_filename(self) -> Path:
        """
        This function returns a `Path` object pointing to the acis observation 
        with point sources removed. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the ACIS image with point 
            sources removed.
        """
        return Path(f"{self.analysis_dir}/acis_nosrc_{self.id}.fits")
    
    @property
    def background_nosrc_filename(self) -> Path:
        """
        This function returns a `Path` object pointing to the high energy data
        file for this observation. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the high energy data file for
            this observation.
        """
        return Path(f"{self.analysis_dir}/background_nosrc_{self.id}.fits")

    @property
    def high_energy_light_curve_file(self) -> Path:
        """
        This function returns a `Path` object pointing to the high energy light
        curve for this observation. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the high energy light curve for
            this observation.
        """
        return Path(f"{self.analysis_dir}/acisI_lcurve_hiE.lc")
    
    @property
    def good_time_interval_light_curve_file(self) -> Path:
        """
        This function returns a `Path` object pointing to the good time interval
        light curve for this observation. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the good time interval light 
            curve for this observation.
        """
        return Path(f"{self.analysis_dir}/acisI_gti_hiE.gti")

    @property
    def acis_nosrc_high_energy_filename(self) -> Path:
        """
        This function returns a `Path` object pointing to the acis observation 
        with point sources removed filtered to high energies. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the ACIS image with point 
            sources removed filtered to high energies.
        """
        return Path(f"{self.analysis_dir}/acisI_nosrc_hiEfilter.fits")



################################################################################
################################################################################
#
#   ███████╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗███████╗
#   ██╔════╝██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║██╔════╝
#   █████╗  ██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║███████╗
#   ██╔══╝  ██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║╚════██║
#   ██║     ╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║███████║
#   ╚═╝      ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝
#
################################################################################
################################################################################

    @property
    def detnam(self) -> str:
        """
        The detnam function defined here queries and returns the DETNAM 
        parameter for this observation. That paremeter is stored in the files
        CIAO downloads for the observation. We use this parameter to tell us 
        about what CCDs were used on the chip for this specific observation.

        To read more: https://cxc.harvard.edu/ciao/threads/ciao_intro/
        ...
        DETNAM: the chip number(s) of the observation, if relevant. The ACIS 
        chip numbering scheme is provided in Table 1. 
        ...
        
        After this passage on the above site, there is a table (Table 1) laying 
        out how each of the returned DETNAM numbers corresponds to a particular
        ACIS CCD. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        str
            A string containing the DETNAM parameter for the observation. A 
            typical ACIS-I result may look like `'ACIS-01235'`.
        """
        if not hasattr(self, '_detnam'):
            keyword_args = {
                "infile": self.level_1_event_file,
                "keyword": "DETNAM",
                "echo": True
            }
            try:
                self._detnam = \
                    ciao.run_command(rt.dmkeypar, **keyword_args)   #type:ignore
            except:
                # To be implemented as as exceptions arise. Raise all for now.
                raise
            logger.debug(f"DETNAM for {self.id} is {self._detnam}")
        return self._detnam

    @property
    def ccds(self) -> list:
        """
        Grabs uses the DETNAM parameter (pythonized above) and parses out the 
        CCDs used for this observation from that result. A typical ACIS-I 
        `detnam` looks like `'ACIS-01235'` (a str). This means that CCDs 0, 1, 
        2, 3, and 6 have events collected. For ACIS-I we will filter this to 
        CCDs -1,1,2, and 3. For now, we record which ccds are used in a list of 
        integers. 
        
        Parameters
        ----------
        None

        Returns
        -------
        list
            A list of integers codifying which CCDs on the ACIS chip were used.
        """
        from cpxt.chandra.ccd import CCD
        if not hasattr(self, '_ccds'):
            chip_ids = self.detnam[5:]
            if len(chip_ids) >= 1:
                self._ccds = [CCD(int(ccd), self) for ccd in chip_ids]
            else:
                logger.error(f"ClusterPyXT could not find the DETNAM parameter"
                    f" in {self.level_1_event_file}. Have you finished stage 1?"
                )
                raise ValueError(f"Error loading DETNAM for {self.id}")
        
        return self._ccds

    @property
    def acis_type(self)  -> ACIS:
        """
        This property returns an `ACIS` `enum` indicating which type of ACIS 
        chip was used during the observation, ACIS-I or ACIS-S. Often times the
        observations will have a CCD or two turned on from the ACIS chip that is
        not the primary one for the observation. To figure out which chip was 
        the main, we add up how many CCD's from each array were used. Whichever
        has more used, we return that one as the main sensor. 

                             ------- -------
                            |       |       |
                            | CCD:0 | CCD:1 |
                            |       |       |
                             ------- -------          <---  ACIS-I
                            |       |       |
                            | CCD:2 | CCD:3 |           
                            |       |       |           
                             ------- -------
                                                       
             ------- ------- ------- ------- ------- -------
            |       |       |       |       |       |       |  
            | CCD:4 | CCD:5 | CCD:6 | CCD:7 | CCD:8 | CCD:9 |  <--- ACIS-S
            |       |       |       |       |       |       |
             ------- ------- ------- ------- ------- -------
        
        Parameters
        ----------
        None
        
        Returns
        -------
        ACIS
            Returns an `ACIS(IntEnum)` (define at top of file) coding which type
            of observation this is, ACIS-I or ACIS-S.
        """
        if not hasattr(self, '_acis_type'):
            acis_i_ccds = [0, 1, 2, 3]
            acis_s_ccds = [4, 5, 6, 7, 8, 9]
            # back_illuminated_ids = [5, 7] # not currently implemented in cpxt
            
            i_ccds = np.array([1 for ccd in self.ccds if ccd.id in acis_i_ccds])
            s_ccds = np.array([1 for ccd in self.ccds if ccd.id in acis_s_ccds])

            # We tallied the total amount, now we need to add them all up to get
            # a total and see which type of observation it was. 
            i_ccd_count = np.sum(i_ccds)
            s_ccd_count = np.sum(s_ccds)

            # As ACIS-I and S observations both may contain some CCD's that 
            # belong to the other. We justed tallied up each sensor's used CCDs.
            # Now we assume whichever sensor has more CCDs used is the primary 
            # sensor for the observation. 
            if i_ccd_count > s_ccd_count:
                self._acis_type = ACIS.I
            elif s_ccd_count > i_ccd_count:
                self._acis_type = ACIS.S
            else:
                logger.error(
                    f"Unable to determine type of observation for {self.id}"
                )
                raise ValueError(f"Unable to determine type of observation.")
        return self._acis_type

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
        This function returns a printable version of the `Observation` object. 
        It overrides pythons `magic` method (signified by double underscore).
        All we need to convey to the user/developer is that it is an observation
        with the specific observation id. All other needed information should
        be able to be garnered from that.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        return f"Observation({self.id})"