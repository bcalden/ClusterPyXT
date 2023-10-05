"""
File: cpxt/chandra/cluster.py
Author: Brian Alden
Date: 29 Mar 2022
Description: This file contains the `ClusterPyXT` `Cluster` class definition. It 
             is the main way we interact with the CIAO cluster files and folders
             as to abstract some of that away.
"""

# Internal imports
from cpxt.chandra import observation as obs
from cpxt.core.config import CPXTConfig
from cpxt.core.stages import Stage

# External imports
from configparser import ConfigParser
from pathlib import Path
import logging
import sys


logger = logging.getLogger(__name__)


class Cluster:
    """
    ClusterPyXT main data class `Cluster`. This class encodes the properties
    of Chandra X-ray galaxy cluster observations needed for the project.
    
    Attributes
    ----------
    name : str
        The name of the cluster. This variable is used to name directories and 
        is the prefix for many filenames. There should be no spaces in the name.
        E.g. Abell 115 -> A115 or Abell_115
    observation_ids : str[]
        A list of chandra observation ids for download. These should all be for 
        the same galaxy cluster. E.g. 3233, 15175, 15144
    hydrogen_column_density : float
        The hydrogen column density for the cluster in /cm^2 . E.g. 5e19
    redshift : float
        The redshift of the cluster. E.g. 0.197
    abundance : float
        The solar abundance of the galaxy cluster. E.g. 0.2
    last_step_completed : cpxt.core.stages.Stage
        This is a state variable indicating the last step of the pypeline that 
        was completed. It indicates where the pypeline will begin if run without 
        any arguments.


    Methods
    -------
    xray_sb_map : Image (@property)

    temperature_map : Image

    pressure_map : Image 

    data_directory : Path
    
    """

    def __init__(self, name: str="", obsids:list=[], 
                 hydrogen_column_density: float=0.0, redshift : float=0.0, 
                 abundance: float=0.0, signal_to_noise: float=50,
                 last_step_completed: Stage=Stage.none) -> None:

        self.name = name
        self.obsids = obsids
        self.hydrogen_column_density = hydrogen_column_density
        self.abundance = abundance
        self.redshift = redshift
        self.signal_to_noise_threshold = signal_to_noise
        self.last_step_completed = last_step_completed

################################################################################
################################################################################
##                                                                            ##
##          ███████╗ ██████╗ ██╗     ██████╗ ███████╗██████╗ ███████╗         ##
##          ██╔════╝██╔═══██╗██║     ██╔══██╗██╔════╝██╔══██╗██╔════╝         ##
##          █████╗  ██║   ██║██║     ██║  ██║█████╗  ██████╔╝███████╗         ##
##          ██╔══╝  ██║   ██║██║     ██║  ██║██╔══╝  ██╔══██╗╚════██║         ##
##          ██║     ╚██████╔╝███████╗██████╔╝███████╗██║  ██║███████║         ##
##          ╚═╝      ╚═════╝ ╚══════╝╚═════╝ ╚══════╝╚═╝  ╚═╝╚══════╝         ##
##                                                                            ##
################################################################################
################################################################################

    @property
    def data_directory(self) -> Path:
        """
        This property returns a `pathlib.Path` object pointing to the cluster's
        data directory. This property checks first to see if the variable 
        `self._data_directory` is intialized yet. If not, create it, then 
        return it. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            Returns a `pathlib.Path` object pointing to the clusters directory.
        """
        if not hasattr(self, "_data_directory"):
            self._data_directory = \
                Path(f"{CPXTConfig().data_directory}/{self.name}")

        return self._data_directory

    @property
    def merged_directory(self) -> Path:
        """
        This property returns a `pathlib.Path` object pointing to the clusters's
        merged event files directory. That directory contains a .lis file with 
        the event files that will be merged.

        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            Returns a `pathlib.Path` object pointing to the clusters directory
            for the merged event file list.
        """

        return Path(f"{self.data_directory}/merged_obs_evt2")
    
    @property
    def main_output_dir(self) -> Path:
        """
        This property returns a `pathlib.Path` object pointing to the clusters's
        main output directory. This is where the merged images are stored. 

        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            Returns a `pathlib.Path` object pointing to the clusters main output
            directory.
        """
        return Path(f"{self.data_directory}/main_output/")

################################################################################
################################################################################
##                                                                            ##
##                     ███████╗██╗██╗     ███████╗███████╗                    ##
##                     ██╔════╝██║██║     ██╔════╝██╔════╝                    ##
##                     █████╗  ██║██║     █████╗  ███████╗                    ##
##                     ██╔══╝  ██║██║     ██╔══╝  ╚════██║                    ##
##                     ██║     ██║███████╗███████╗███████║                    ##
##                     ╚═╝     ╚═╝╚══════╝╚══════╝╚══════╝                    ##
##                                                                            ##
################################################################################
################################################################################

    @property
    def config_file(self) -> Path:
        """
        This return a `pathlib.Path` object poinging to the `ClusterPyXT` 
        coinfiguration file for this particular cluster. It contains information
        such as the observations ids used, hydrogen column density, redshift
        and the desired signal to noize ratio threshold. This information is 
        saved using the `configparser` class in similar format to a `.ini` file.
        The file is saved at the folloing location (<name> = cluster name):

        /path/to/clusterpyxt-data/<name>/<name>_pypeline_config.ini
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `pathlib.Path` object pointing to the clusters config file. 
        """
        if not hasattr(self, "_config_file"):
            self._config_file = get_cluster_config_file(self.name)
        return self._config_file

    @property
    def acis_ccd_evt2_file_list(self) -> list:
        """
        This function returns a list of the level 2 event files for each CCD,
        from each observation, for the cluster. The individual objects of the 
        list are all strings pointing to the respective files. This is usefull 
        when you want to perform an operation on every level 2 event file for 
        the particular cluster. 
        
        Parameters
        ----------
        None

        Returns
        -------
        list
            A list of `pathlib.Path` objects pointing to every `acis_ccd#.fits`
            from each of the clusters observations used.
        """
        all_obs_files = [obs.acis_ccd_fits_files for obs in self.observations]

        # Right now, `all_obs_files` is a list of lists. We want to flatten it. 
        # We will use a list comprehension to do it quickly.
        # flat_list = [item for sublist in list_of_lists for item in sublist]
        files = [str(ccd_file) for obs_files in all_obs_files \
                                                 for ccd_file in obs_files]
        
        
        return files

    @property
    def acis_ccd_evt2_filenames(self) -> Path:
        """
        This property returns a `Path` object pointing to the file containing a
        list of the clusters level 2 event files per observation per CCD. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            A `Path` object pointing to the text file containing each CCD level 
            2 event file for each observation for the cluster. This file is used
            to pass all of the level 2 event files we care about to CIAO 
            commands with a much simpler call. 
        """
        return Path(f"{self.data_directory}/acis_ccd_evt2.lis")

    @property
    def merged_acis_events(self) -> Path:
        """
        The is property returns a `Path` object that points to the merged events
        file. This file is created by merging the files listed in 
        `acis_ccd_evt_filenames`
        
        Parameters
        ----------
        <var> : <type>
            <description>
        
        Returns
        -------
        <type>
            <description>
        """
        # This used to be acisI.fits in older versions of the pipeline.
        return Path(f"{self.data_directory}/merged_acis_events.fits")
    

    @property
    def sources_file(self) -> Path:
        """
        The is property returns a `Path` object that points to the sources file.
        This file should contain the regions identified as point sources 
        outside of the cluster that should be removed from the analysis. Failure
        to remove these sources will have an impact on the spectral fitting.
                
        Parameters
        ----------
        <var> : <type>
            <description>
        
        Returns
        -------
        <type>
            <description>
        """
        return Path(f"{self.data_directory}/sources.reg")
    
    @property
    def exclude_file(self) -> Path:
        """
        The is property returns a `Path` object that points to the exclude file.
        This file should contain the region around the peak of cluster emission
        to be excluded from the deflaring process.
                
        Parameters
        ----------
        <var> : <type>
            <description>
        
        Returns
        -------
        <type>
            <description>
        """
        return Path(f"{self.data_directory}/exclude.reg")
    
    @property
    def xray_sb_map_filename(self) -> Path:
        """
        The is property returns a `Path` object that points to the merged X-ray 
        surface brightness map for the cluster. Use this file to identify the 
        point sources to be removed from the analysis.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            Path object pointing to the merged X-ray surface brightness map.
        """
        return Path(f"{self.data_directory}/{self.name}_broad_flux.img")

    @property
    def xray_sb_nosrc_map_filename(self) -> Path:
        """
        The is property returns a `Path` object that points to the merged X-ray 
        surface brightness map for the cluster with the point sources removed.
                
        Parameters
        ----------
        None
        
        Returns
        -------
        Path
            Path object pointing to the merged X-ray surface brightness map with
            point sources removed.
        """
        return Path(f"{self.main_output_dir}/" \
                    f"{self.name}_xray_surface_brightness_nosrc.fits")

################################################################################
################################################################################
##                                                                            ##
## ███████╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗███████╗ ##
## ██╔════╝██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║██╔════╝ ##
## █████╗  ██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║███████╗ ##
## ██╔══╝  ██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║╚════██║ ##
## ██║     ╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║███████║ ##
## ╚═╝      ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝ ##
##                                                                            ##
################################################################################
################################################################################

    @property
    def observations(self) -> list:
        """
        This function returns a list where each element is an `Observation` 
        class encoding each fo the particular clusters observation ids. 
        
        Parameters
        ----------
        <var> : <type>
            <description>
        
        Returns
        -------
        <type>
            <description>
        """
        if not hasattr(self, "_observations"):
            self._observations = \
                [obs.Observation(obsid, cluster=self) for obsid in self.obsids]
        
        return self._observations
    
    def save(self):
        """
        Uses the ConfigParser to save the cluster configuration paramters. 
        These parameters are generated by the `__iter__` method defined above.
        They are then turned into a dictionary and passed to the `ConfigParser`
        object to write them to `self.config_file` (defined as a property below)
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        cluster_config = ConfigParser()
        cluster_dict = dict(self)
        
        # Convert the list of integers to a string to save/load from disk 
        cluster_dict['obsids'] = ','.join([f"{obs}" for obs in self.obsids])
        cluster_dict['last_step_completed'] = \
            f"{self.last_step_completed.value}"
        cluster_config['cluster'] = cluster_dict

        self.config_file.parent.mkdir(parents=True, exist_ok=True)
        try:
            with open(str(self.config_file), 'w') as f:
                cluster_config.write(f)
        except:
            raise
        logger.info(f'Wrote config file: {str(self.config_file)}')


    def __iter__(self):
        """
        Implementing the __iter__ 'magic function' allows us to return an 
        python iterable if necessary. This implementation is necesssary for
        the `dict(self)` call in the `save` function.
        
        """
        yield 'name', self.name
        yield 'obsids', str(self.obsids)
        yield 'hydrogen_column_density', str(self.hydrogen_column_density)
        yield 'redshift', str(self.redshift)
        yield 'abundance', str(self.abundance)
        yield 'signal_to_noise', str(self.signal_to_noise_threshold)
        yield 'last_step_completed', str(self.last_step_completed)

        return

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
        This function returns a printable version of the `Cluster` object. 
        It overrides pythons `magic` method (signified by double underscore).
        All we need to convey to the user/developer is that it is an cluster
        with the specific number of observations. Listing the specific obsids
        associated with this project may seem find when len(obsids) < 5 to 10. 
        Some projects contain > 30 observations and listing out each specific
        observation id with the `__repr__` function would be cumbersome.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        return f"Cluster({self.name}, {len(self.obsids)} Obs.)"

################################################################################
################################################################################
##                                                                            ##
##    ██████╗██╗      █████╗ ███████╗███████╗                                 ##
##   ██╔════╝██║     ██╔══██╗██╔════╝██╔════╝                                 ##
##   ██║     ██║     ███████║███████╗███████╗                                 ##
##   ██║     ██║     ██╔══██║╚════██║╚════██║                                 ##
##   ╚██████╗███████╗██║  ██║███████║███████║                                 ##
##    ╚═════╝╚══════╝╚═╝  ╚═╝╚══════╝╚══════╝                                 ##
##                                                                            ##
##           ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗    ##
##           ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝    ##
##           ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗    ##
##           ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║    ##
##           ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║    ##
##           ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝    ##
##                                                                            ##
################################################################################
################################################################################       

def load_cluster(cluster_name: str) -> Cluster:
    """
    This function takes a `cluster_name` string, attempts to find the 
    configuration file, and then make a `Cluster` object out of it. 
    
    Parameters
    ----------
    cluster_name : str
        The name of the cluster we are trying to load. This is also the 
        directory name. 
    
    Returns
    -------
    Cluster
        A `cpxt` Cluster object pythonizing the desired cluster
    """
    # First, we need to find the cluster config filename
    cluster_config_file = get_cluster_config_file(cluster_name)
    
    # Next we use a `ConfigParser` to read the config file. 
    cluster_config = ConfigParser()
    cluster_config.read(cluster_config_file)

    # Now we try to read that config file. This may fail if the file doesn't
    # exist or does not contain that particular key (`cluster`). Either way
    # the ConfigParser class will throw a KeyError. If so, log it and `raise`.  
    try:
        cluster_dict = dict(cluster_config['cluster'])
    except KeyError:
        logger.error("Could not load {cluster_config_file}.")
        raise
    

    # Now that the configuration file is loaded, we grab the values from it
    name = cluster_dict['name']
    
    try:                                         # Check the configuration file
        raw_obsids = cluster_dict['obsids']          # for both possible ways the
    except KeyError:                             # observation ids may be 
        raw_obsids = cluster_dict['observation_ids'] # stored. 
    
    # Convert a long string of observation ids seperated by a comma into a list
    obsids = [int(obsid) for obsid in raw_obsids.split(',')]
    
    hydrogen_column_density = float(cluster_dict['hydrogen_column_density'])
    redshift = float(cluster_dict['redshift'])
    abundance = float(cluster_dict['abundance'])
    last_step_completed = Stage(int(cluster_dict['last_step_completed']))
    
    try:
        sig_to_noise_thresh = float(cluster_dict['signal_to_noise'])
    except KeyError:             # Early `cpxt` versions did not capture this
        sig_to_noise_thresh = 50.0 # value and just used `50`. 
    
    # Actually make the `Cluster` object then return it.
    cluster = Cluster(name=name,
                      obsids=obsids,
                      hydrogen_column_density=hydrogen_column_density,
                      redshift=redshift,
                      abundance=abundance,
                      signal_to_noise=sig_to_noise_thresh,
                      last_step_completed=last_step_completed
                      )
    return cluster
    
    
def get_cluster_config_file(cluster_name: str) -> Path:
    """
    This function tries to find the `cpxt` cluster configuration file for the
    passed `cluster_name`. First, it just checks the data
    
    Parameters
    ----------
    cluster_name : str
        The name of the cluster we are finding the configuration file for.
    
    Returns
    -------
    Path
        A Path object encoding the cluster configuration filename. 
    """
    data_dir = CPXTConfig().data_directory
    filename = f"{data_dir}/{cluster_name}/{cluster_name}_pypeline_config.ini"
    cluster_config_filename = Path(filename)
    
    return cluster_config_filename

