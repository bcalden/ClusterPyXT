"""
File: cpxt/__main__.py
Author: Brian Alden
Date: 28 Mar 2022
Description: ClusterPyXT (`cpxt`) main file. First, create/loads the `cpxt`
             configuration files. Next, launch the main window of the app.
             Alternatively, process the command line arguments and execute
             accordingly.
"""

# Internal imports
from cpxt.cli.arguments import process_cli_args
from cpxt.chandra.ciao import get_ciao_version
from cpxt.gui.main_window import MainWindow
from cpxt.core import config as cfg

# External imports
from PyQt5.QtWidgets import QApplication, qApp
import logging
import sys
import os


def main():
    """
    Parameters
    ----------
    N/A

    Returns
    -------
    N/A

    This method runs the main program. First it checks to see if the user
    passed arguments at the command line. If none were passed, it is
    assumed the GUI application is desired and is launched. If arguments
    are passed at the command line, deal with them and execute accordingly.
    """

    logging.basicConfig(filename=cfg.CPXTConfig.log_filename, 
                        level=logging.DEBUG
        )
    logger = logging.getLogger(__name__) 
    logger.debug("Starting ClusterPyXT.")
    
    # Make sure we are in the CIAO conda environment before continuing. 
    try:
        ciao_version = get_ciao_version()
        logger.debug(f"Found CIAO Version info: {ciao_version}")
    except:
        logger.critical(
            "Must be in a conda environment for CIAO in order to use. Start " 
            "the CIAO conda environment and retry. See the following url for "
            "CIAO installation instructions and environment creation. Then "
            "run ClusterPyXT again:\n"
            "https://cxc.harvard.edu/ciao/threads/ciao_install_conda/")
        sys.exit(-1)
    
    allowed_versions = ['CIAO 4.13', 'CIAO 4.14']
    if ciao_version not in allowed_versions:
        logger.warning(f"The version of CIAO loaded, {ciao_version}, has not "
            "been tested with this version of ClusterPyXT. Unexpected errors "
            "may occur.")

    if len(sys.argv) == 1: 
        # No arguments other than the program name passed at command line. 
        # Start the main GUI. 
        app = QApplication(sys.argv)
        qApp.setApplicationName("CPXT")
        qApp.setApplicationDisplayName("ClusterPyXT")
        win = MainWindow()
        win.show()
        app.exec_()
    else: 
        # Arguments were passed at the command line and need to be processed
        # Assume command line interface (CLI).
        process_cli_args()

    logger.debug("Exiting ClusterPyXT")


if __name__ == "__main__":
    sys.exit(main())
    # This call to main is wrapped with sys.exit() allowing for scripting.
    # See the following link for further reason for the sys.exit wrapper: 
    #
    # https://docs.python.org/3/library/__main__.html#packaging-considerations
    