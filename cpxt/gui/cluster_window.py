"""
File: cpxt/gui/cluster_window.py
Author: Brian Alden
Date: 28 Mar 2022
Description: The main window of an individual cluster. Also used to create and 
             add the cluster to `cpxt` for processing. 
"""

# Internal imports
from cpxt.stages import stage_one, stage_two, stage_three,\
    stage_four, stage_five, stage_spectral_fitting
from cpxt.chandra.cluster import Cluster, load_cluster
from cpxt.gui.ui.cluster_ui import Ui_ClusterUI as ClusterUI
from cpxt.core.config import CPXTConfig
from cpxt.core.stages import Stage

# External imports
from PyQt5 import QtWidgets as qt
from PyQt5.QtCore import QObject
import logging
from astroquery.ned import Ned
from lxml import etree
from requests import get

logger = logging.getLogger(__name__)

class ClusterWindow(qt.QMainWindow, ClusterUI):
    """
    The `ClusterWindow` class displays a gui wrapping all of the processing
    capabilities of `ClusterPyXT` in to one main processing window. Initial 
    cluster creation is done in this window. After creation, all stages are 
    started from this window some launching commandline processes and some
    spawinging their own window for further user input. 
    
    Attributes
    ----------
    name : str
        The name of the cluster being passed. If no value is passed, 
        setup the GUI for cluster initialization.
    parent : QMainWindow
        The main `cpxt` window, a cpxt.gui.main_window.MainWindow object, that 
        instatiated this window. 
    
    Methods
    -------
    cluster_name_entered : None
        When the user inputs a cluster name, query NED to get properties of
        the cluster automatically. 
    run_stage_1 : bool
        Downloads the observation data, merges the observations, finds and 
        merges the backgrounds. 
    run_stage_2 : bool
        Removes the point sources from each observation. Filters high energy
        flares. 
    run_stage_3 : bool
        Extracts response files, RMF & ARF files are generated.
        RMF = https://cxc.cfa.harvard.edu/ciao/dictionary/rmf.html
        ARF = https://cxc.cfa.harvard.edu/ciao/dictionary/arf.html
    run_stage_5 : bool
        Calculates the adaptive circular bins, generates a bin map, and 
        calculates the exposure corrections necessary.
    run_stage_spectral : bool
        Launches the sherpa spectral fitting process to generate the temperature
        estimates. 
    make_final_products : None
        Uses the spectral fitting results to make temperature and pressure maps
        along with their resulting error maps. 

    """


    def __init__(self, name: str="", parent: QObject=QObject()):
    
        super().__init__()
        self.setupUi(self)
        
        self.setWindowTitle("New Cluster")
        self.parent = parent
        
        # When the user leaves the cluster name field, and astroquery is set
        # to be used in the settings, search for nH and redshift automatically.
        self.name_text.editingFinished.connect(self.cluster_name_entered)

        # Now that the buttons are created, we need to connect their `clicked` 
        # functionality to a function in our code. 
        self.save_update_button.clicked.connect(self.save_cluster)
        self.stage_1_button.clicked.connect(self.run_stage_1)
        self.stage_2_button.clicked.connect(self.run_stage_2)
        self.stage_3_button.clicked.connect(self.run_stage_3)
        self.stage_4_button.clicked.connect(self.run_stage_4)
        self.stage_5_button.clicked.connect(self.run_stage_5)
        self.spectral_fitting_button.clicked.connect(self.run_stage_spectral)
        self.products_button.clicked.connect(self.make_final_products)
        self.name_text.setFocus()

        # Do we need to create a new cluster or continue an already created one
        if name:
            try:
                self.cluster = load_cluster(name)
            except:
                logger.error(f"Tried to load cluster for {name} and failed.")
                raise
            logger.debug(f"Successfully loaded {self.cluster}")
            # We loaded the cluster successfully so now we can populate the GUI

            self.setWindowTitle(f'{self.cluster.name}')
            self.name_text.setText(self.cluster.name)
            
            # The obsids are stored in `cluster` as a list. We want a string 
            # with each obsid seperated by a " " (space)
            self.obsid_text.setPlainText(
                " ".join([f"{obsid}" for obsid in self.cluster.obsids])
                )
            
            # Set the rest of the attributes
            self.nH_text.setText(f"{self.cluster.hydrogen_column_density:0.3e}")
            self.redshift_text.setText(f"{self.cluster.redshift}")
            self.abundance_text.setText(f"{self.cluster.abundance}")
            self.sig_to_noise_text.setText(
                f"{self.cluster.signal_to_noise_threshold}"
                )
            self.stage = Stage(self.cluster.last_step_completed)
            self.update_buttons()
        else:
            stage_one.print_stage_1_prep()
            


    def cluster_name_entered(self) -> None:
        """
        Checks the `CPXTConfig` to see if automatic query is enabled. If so, 
        check for nH and redshift. If found, populate automatically. It uses
        astroquery and the NASA Extragalactic Database (NED). 

        TODO: Update this function so it runs asynchronously and thus allows
        the user to continue using the interface while this function runs. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        # The cluster name was just entered. Lets set it to the main windows
        # name property for easier access

        self.name = self.name_text.text()

        if not CPXTConfig().astro_query:
            return

        # This next `try` block uses astroquery to check NED for redshift
        # and RA/Dec information. 
        # 
        # Next, it uses a `request` object to query NASA HEASARC's nH calclator
        # website with the appropriate RA and Dec to estimate the clusters 
        # hydrogen column density (nH). 
        # 
        # The queried values are then used to populate the GUI for the user. 
        #
        # Note: The comments `#type: ignore` are for the automated type checker
        try: 
            ned_cluster = Ned.query_object(self.name)
            redshift = ned_cluster['Redshift'][0]  # type: ignore 
            ra = ned_cluster['RA'][0]              # type: ignore
            dec = ned_cluster['DEC'][0]            # type: ignore

            # This web page and query string are not guaranteed to remain valid 
            # as they are completely outside our control. 
            html = get(make_url(ra,dec)).content
            dom = etree.HTML(html, parser=None)                 
            raw_text = dom.xpath('//b[contains(text(), "Average nH")]')[0].text
            nh = float(raw_text.split(' ')[-1])

            # The above redshift query can gracefully fail and return a 0. 
            # If so, it is likely the rest failed and we shouldn't use the data.
            # If it succeeded, lets use the data and populate the GUI with it.
            if redshift != 0:
                self.redshift_text.setText(f"{redshift}")
                self.nH_text.setText(f"{nh:0.2e}")
        except:
            pass
    
    def save_cluster(self):
        """
        Makes a new `Cluster` object out of the cluster attributes the user 
        entered. After the object is created we save the cluster config data by
        calling the `Cluster.save` function.
        """
        name = self.name_text.text()
        obsids = \
            [obsid for obsid in self.obsid_text.toPlainText().split(' ')]
        hydrogen_column_density = float(self.nH_text.text())
        redshift = float(self.redshift_text.text())
        abundance = float(self.abundance_text.text())
        sig_to_noise_thresh = float(self.sig_to_noise_text.text())
        last_step_completed = Stage.init
        
        cluster = Cluster(
            name=name,
            obsids=obsids,
            hydrogen_column_density=hydrogen_column_density,
            redshift=redshift,
            abundance=abundance,
            signal_to_noise=sig_to_noise_thresh,
            last_step_completed=last_step_completed
        )
        self.cluster = cluster
        self.cluster.save()
        self.parent.load_clusters()
        self.setWindowTitle(f'{name}')
        self.update_buttons()


    def run_stage_1(self):
        """
        Runs stage one on the cluster shown in this window. Downloads the 
        observation data for each OBSID. Merge the observations to a combined
        surface brightness map. Last, finds and merges the backgrounds 
        in the CALDB.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        stage_one.run_on(cluster=self.cluster)
        self.update_stage(Stage.one)


    def run_stage_2(self):
        stage_two.run_on(cluster=self.cluster)
        self.update_stage(Stage.two)
    
    def run_stage_3(self):
        pass

    def run_stage_4(self):
        pass

    def run_stage_5(self):
        pass

    def run_stage_spectral(self):
        pass
    
    def make_final_products(self):
        pass

    def update_stage(self, stage: Stage):
        """
        """
        self.cluster.last_step_completed = stage
        self.cluster.save()

        # update the buttons to reflect the new stage available
        self.update_buttons()


    def update_buttons(self):
        """
        This function is called to enable/disable the save/stage buttons after
        every stage runs. 
        """
        gui_buttons = [
            self.save_update_button, 
            self.stage_1_button, 
            self.stage_2_button, 
            self.stage_3_button, 
            self.stage_4_button, 
            self.stage_5_button, 
            self.spectral_fitting_button, 
            self.products_button
        ]

        # Disable all the buttons, then enable them based on the last stage done    
        _ = [button.setEnabled(False) for button in gui_buttons]

        # This dictionary alleviates the need for a long set of if conditionals.
        # Simply indexing based on stage returns the list of buttons that should
        # be enabled. 
        buttons_to_enable = {
            Stage.init: gui_buttons[1:2],
            Stage.one: gui_buttons[1:3],
            Stage.two: gui_buttons[1:4],
            Stage.three: gui_buttons[1:5],
            Stage.four: gui_buttons[1:6],
            Stage.five: gui_buttons[1:7],
            Stage.spec: gui_buttons
        }
        
        # enable the buttons based on the last completed stage
        stage = self.cluster.last_step_completed
        _ = [button.setEnabled(True) \
                for button in buttons_to_enable[stage]]



def make_url(ra, dec):
    if dec > 0:
        return f"https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl?"\
        f"Entry={ra}%2C+{dec}&NR=GRB%2FSIMBAD%2BSesame%2FNED&CoordSys=" \
        "Equatorial&equinox=2000&radius=0.1&usemap=0"
    return "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl?Entry="\
        f"{ra}%2C{dec}&NR=GRB%2FSIMBAD%2BSesame%2FNED&CoordSys="\
            "Equatorial&equinox=2000&radius=0.1&usemap=0"


