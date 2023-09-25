"""
File: cpxt/gui/main_window.py
Author: Brian Alden
Date: 28 Mar 2022
Description: The main ClusterPyXT window. Displays all clusters saved within
             the `cpxt` data directory.
"""

# Internal imports
from cpxt.core.config import CPXTConfig
from cpxt.gui.cluster_window import ClusterWindow
from cpxt.gui.config_window import ConfigWindow
from cpxt.gui.ui.main_ui import Ui_MainUI as MainUI

# External imports
from PyQt5 import QtWidgets as qt
from PyQt5 import QtWidgets as qtw
from glob import glob
import sys


class MainWindow(qt.QMainWindow, MainUI):
    """
    The main ClusterPyXT window. Contains a list of clusters at various stages
    of completion. Allows the user to select a cluster project to continue an
    already created cluster, or create a new cluster project. Can also change
    system settings from this window. 

    Needs to implement: Add a menu bar to access an about box containing 
    citation information and an ability to launch a `cpxt` configuration 
    window (also needs to be implemented).

    Attributes
    ----------
    cluster_list_box : qt.QListWidget
        Main list box for the window to be populated with cluster names.
    continue_button : qt.QPushButton
        Button allowing the user to continue the selected cluster
    new_cluster_button : qt.QPushButton
        Button allowing the user to add a new cluster into `ClusterPyXT`.
    windows : list
        A list of open windows. Needed to close all open windows on exit.

    Methods
    -------
    load_clusters() -> None
        Scans the data directory and generates a list of `cpxt` clusters
    continue_button_clicked() -> None
        Opens a `ClusterWindow` for the selected cluster when clicked.
    new_button_clicked() -> None
        Opens a blank `ClusterWindow` to intialize a new cluster.

    Notes
    -----
    This class instantiates the MainUI class from `main_ui.py`. The `main_ui.py`
    file is generated automatically from `main_ui.ui`, a form made using 
    `Qt Designer`. Edits should be done there, then call (from the ui folder):

        `pyuic6 -o main_ui.py main_ui.ui`

    The only edit made to the `main_ui.py` file that is made from then, is for 
    the type checker to ignore it. 

    """

    def __init__(self) -> None:
        """
        The class initializer for `MainWindow`. This function creates the main
        GUI containing a list of all clusters `ClusterPyXT` is tracking. The 
        instantiation of this window is also the main application being tracked
        by the OS. When this window is closed, the whole program should close. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        super().__init__()
        self.setupUi(self)

        # Create a list of windows so we can close them all when the program
        # exits. It starts empty. Each cluster window opened will be added 
        self.windows = []
        
        self.cluster_list_box.itemDoubleClicked.connect( # Double clicking the
            self.continue_button_clicked                 # list box simulates
        )
        
        self.continue_button.clicked.connect(self.continue_button_clicked)
        self.new_cluster_button.clicked.connect(self.new_button_clicked)

        # Setup the menu bar
        self.menubar.triggered[qtw.QAction].connect(self.process_menu)

        # Populate the main list box with cluster names
        self.load_clusters()

        self.show()


    def load_clusters(self) -> None:
        """
        This function checks the `cpxt` data directory (stored in the config 
        file) for valid `cpxt` clusters. It does so using glob to pattern
        match the .ini file for each folder in the directory. The names of
        each cluster are then populated into the main list box of the window.
        
        Parameters
        ----------
        N/A
        
        Returns
        -------
        N/A
        """
        self.cluster_list_box.clear()
        try:
            cpxt_data_dir = CPXTConfig().data_directory
        except KeyError:
            config_win = ConfigWindow(self)
            config_win.show()
            config_win.exec()
            try:
                cpxt_data_dir = CPXTConfig().data_directory
            except:
                raise

        ini_files = glob(f"{cpxt_data_dir}/*/*.ini")
        
        # convert the long filename to cluster name for each file in ini_files
        # i.e. /path/data/A115/A115_pypeline_config.ini -> A115
        cluster_names = \
            [filename.split('/')[-1].split('_')[0] for filename in ini_files]
    
        # Populate the main list box in the window with each clusters name
        _ = [self.cluster_list_box.addItem(name) for name in cluster_names]

        # Called to redraw the newly populated list box. Without this call
        self.cluster_list_box.update()
    

    def continue_button_clicked(self) -> None:
        """
        When the continue button is clicked, or a user double clicks a cluster,
        this function is called. It figures out which cluster was selected and
        launches a `ClusterWindow` window for that cluster. 
        
        Parameters
        ----------
        N/A
        
        Returns
        -------
        N/A
        """
        
        selected_cluster = self.cluster_list_box.selectedItems()[0].text()
        cluster_window = ClusterWindow(name=selected_cluster, parent=self)
        self.windows.append(cluster_window)
        cluster_window.show()


    def new_button_clicked(self) -> None:
        """
        User clicked the `New Cluster` button in the GUI. This function 
        launches a `ClusterWindow` window without passing a cluster name.
        When a `ClusterWindow` is instantiated without the `name` parameter 
        passed, it assumes you are initializing a new cluster. 

        Parameters
        ----------
        N/A
        
        Returns
        -------
        N/A
        """

        new_cluster_window = ClusterWindow(parent=self)
        self.windows.append(new_cluster_window)
        new_cluster_window.show()

    def process_menu(self, menu_item):
        actions = {
            '&Settings': self.show_config_window,
            'E&xit': sys.exit
        }
        if menu_item.text() in actions:
            actions[menu_item.text()]()
    
    def show_config_window(self):
        win = ConfigWindow(parent=self)
        win.show()
        self.windows.append(win)