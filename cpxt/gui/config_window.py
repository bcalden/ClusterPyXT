"""
File: cpxt/gui/config_window.py
Author: Brian Alden
Date: 4 Apr 2022
Description: This file handles the GUI for the `cpxt` configuration settings. 
             It is run the first time the user runs `ClusterPyXT` and if the 
             end user wants to change settings from the main window. 
"""

# Internal imports
from cpxt.gui.ui.config_ui import Ui_ConfigUI as ConfigUI
from cpxt.core.io import str_to_bool
from cpxt.core.config import CFGKey
from cpxt.core import config as cfg

# External imports
from PyQt5.QtWidgets import QFileDialog, QDialog
from PyQt5 import QtWidgets as qt
from PyQt5.QtCore import QObject
from multiprocessing import cpu_count
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class ConfigWindow(qt.QDialog, ConfigUI):
    """
    The `ConfigWindow` class handles the GUI aspects of capturing configuration
    settings neccessary for the operation of `ClusterPyXT`. These settings are
    how many observations we can download at once, how many cpu cores to use
    during parallel processing tasks, and most importantly, where to save all of
    the data produced while running the software. 
    
    The configuration window has 5 rows of settings. Each row of settings
    contains usually contains multiple widgets such as a label, texbox, or
    button. These widgets are then used to populate a horizontal layout 
    for each row. At the end, we will add all of the horizontal layouts
    on to one big vertial layout.

    Attributes
    ----------
    parent : QObject
        The parent window opening this dialog. As implemented, this is the 
        `MainWindow` instantiation.
    
    Methods
    -------
    find_directory : None
        This function displays a file dialog box for the user to select a 
        directory.
    save : None
        Saves the the currently displayed settings to the `CPXTConfig.
    load : bool
        Returns whether or not it was able to successfully load the saved 
        configuration file. It will fail the first time run, triggering the 
        first creation of the file. All times after that should return True. 
    """

    def __init__(self, parent: QObject):
        super().__init__(parent)
        self.setupUi(self) # uses the ui/config_ui.py file to setup the object

        self.find_directory_button.clicked.connect(self.find_directory)

        # Setting up the number of cores to use widgets and layout
        cores_label = f"Number of CPU cores to use for parallel " \
                          f"processing? (Max: {cpu_count()}): "
        self.num_cores_label.setText(cores_label)
        
        self.save_button.clicked.connect(self.save)


        # We've made all of the GUI widgets at this point. Now its time to 
        # load the saved configuration data. If that call fails, we will set
        # up everything with default values. 
        if not self.load():
            # We were unable to load the `cpxt` config file therefore we need
            # to fill in some default values.
            
            self.directory_text.setText(f"{Path.home()}/cpxtdata")
            self.astro_query_check.setChecked(True)
            half_cpu_count = int(cpu_count())//2                  # type: ignore
            self.num_cores_text.setText(f"{half_cpu_count}")
            self.log_level_combo.setCurrentIndex(1)
        
    
    def find_directory(self) -> None:
        """
        Displays a message dialog and file dialog prompting the user to select
        a location for data storage. Updates the text in the GUI with the new
        location. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        dialog = QFileDialog()
        dialog.setFileMode(QFileDialog.DirectoryOnly)
        
        if dialog.exec_() == QDialog.Accepted: # Display the file open dialog
            # If the user selected a directory, set the textbox to it
            self.directory_text.setText(dialog.selectedFiles()[0])
        


    def save(self) -> None:
        """
        Saves all of the ClusterPyXT settings configured in this window. After
        saving the file, it closes the this window and continues on with the
        program.
        
        Files Generated
        ---------------
        /usershomedir/.config/cpxt/cpxt.cfg
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        data_directory = self.directory_text.text()
        astro_query = self.astro_query_check.isChecked()
        num_downloads = int(self.obs_download_combo.currentText())
        num_cores = int(self.num_cores_text.text())
        log_level = int(self.log_level_combo.currentText())

        cfg.write_cpxt_config(data_directory,
                              astro_query,
                              num_downloads,
                              num_cores,
                              log_level
                              )
        self.close()


    
    def load(self) -> bool:
        """
        This function attempts to load the `ClusterPyXT` settings from the saved
        configuration file on disk. If succesfull, it sets the values on the GUI
        to the saved parameters and returns True. If not, it returns False.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        bool
            Returns a bool characterizing whether or not it was able to load the
            configuration settings from disk.
        """

        try:
            config = cfg.CPXTConfig().load_cpxt_config()
            self.directory_text.setText(config[CFGKey.data_directory.value])
            self.astro_query_check.setChecked(
                str_to_bool(config[CFGKey.astro_query.value])
                )
            self.num_cores_text.setText(config[CFGKey.num_cores.value])
            self.obs_download_combo.setCurrentText(
                config[CFGKey.num_downloads.value]
                )
            self.log_level_combo.setCurrentText(config[CFGKey.log_level.value])
            return True
        except KeyError:
            return False
        except:
            raise # to be implemented as uncaught exceptions arise