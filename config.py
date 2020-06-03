import os
import cluster
import pypeline_io as io
import configparser
import glob
import pypeline_io as io
from PyQt5.QtWidgets import QFileDialog, QDialog
from PyQt5 import QtWidgets


CONFIG_FILENAME = io.get_path("{clusterpy_dir}/pypeline_config.ini".format(
    clusterpy_dir=os.path.dirname(__file__)
    ))
print(CONFIG_FILENAME)
class ClusterPyConfig():
    def __init__(self, filename=CONFIG_FILENAME, data_directory="", ciao_directory=""):
        self.filename = filename
        self._data_directory = data_directory
        self._ciao_directory = ciao_directory
        self._current_cluster = ""
        if not io.file_exists(filename):
            io.make_initial_data_dir(self.data_directory)
            self.write_system_configuration()
        self.read_system_configuration()
        
    def get_data_dir(self, parent=None):
        alert = QtWidgets.QMessageBox()
        alert.setText("Please select the directory you would like to use for data storage.")
        alert.exec_()
        dialog = QFileDialog()
        dialog.setFileMode(QFileDialog.DirectoryOnly)
        if dialog.exec_() == QDialog.Accepted:
            self.data_directory = dialog.selectedFiles()[0]
    
    def get_ciao_dir(self):
        alert = QtWidgets.QMessageBox()
        alert.setText("Please select the CALDB installation directory. If installed using conda, likely in ../anaconda3/pkgs/caldb_main-<version>/")
        alert.exec_()
        dialog = QFileDialog()
        dialog.setFileMode(QFileDialog.DirectoryOnly)
        if dialog.exec_() == QDialog.Accepted:
            self.ciao_directory = dialog.selectedFiles()[0]
    
    @property
    def data_directory(self):
        return self._data_directory

    @data_directory.setter
    def data_directory(self, data_dir):
        self._data_directory = data_dir
        io.make_initial_data_dir(data_dir)
        self.write_system_configuration()


    @property
    def ciao_directory(self):
        return self._ciao_directory

    @ciao_directory.setter
    def ciao_directory(self, ciao_dir):
        self._ciao_directory = ciao_dir
        self.write_system_configuration()
    
    
    @property
    def current_cluster(self):
        return self._current_cluster

    @current_cluster.setter
    def current_cluster(self, cluster_name):
        self._current_cluster = cluster_name
        self.write_system_configuration()
    
    def read_system_configuration(self):
        config = configparser.ConfigParser()

        config.read(self.filename)
        try:
            config_dict = dict(config['pypeline_config'])
        except KeyError as err:
            print("Configuration file not setup correctly. Try deleting pypeline_config.ini and restarting."
             "If problem continues, please contact brian.alden@colorado.edu with the error message.")
            raise
        
        self._data_directory = config_dict['data_directory']
        self._ciao_directory = config_dict['ciao_directory'] 
        self.current_cluster = config_dict['current_cluster_name']

    def write_system_configuration(self):
        config = configparser.ConfigParser()

        config_dict = {
            'data_directory': os.path.normpath(self.data_directory),
            'ciao_directory': os.path.normpath(self.ciao_directory),
            'current_cluster_name': self.current_cluster
        }

        config['pypeline_config'] = config_dict

        with open(self.filename, 'w') as configfile:
            config.write(configfile)

sys_config = ClusterPyConfig()


def get_user_input(prompt):
    user_input_good = False
    user_input = None
    while not user_input_good:
        user_input = input(prompt)
        y_no_flag = input("You entered {}. Is this correct? (y/n):".format(user_input))

        if y_no_flag.lower() in ['y', 'yes']:
            user_input_good = True

    return user_input


def get_cluster_configs(data_dir=sys_config.data_directory):
    dirs = os.listdir(data_dir)
    configuration_files = []
    cluster_names = []
    for directory in dirs:
        #print("{}/{}".format(data_dir,directory))
        config_file = glob.glob("{data_dir}/{directory}/*_pypeline_config.ini".format(
            data_dir=data_dir,
            directory=directory
        ))
        if config_file:
            configuration_files.append(config_file[0])
            cluster_names.append(get_cluster_name_from_config_file(config_file[0]))
    return list(zip(cluster_names, configuration_files))

def get_cluster_name_from_config_file(config_file):
    return os.path.basename(os.path.split(config_file)[0])

