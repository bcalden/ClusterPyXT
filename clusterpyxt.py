# version 0.9.0a
# Jan 14, 2020

import sys
import argparse
import pypeline_io as io
import cluster
from errors import ClusterPyError
import config
from PyQt5 import QtCore, QtGui, QtWidgets
from enum import IntEnum
import acb
from astropy.io import fits
import data_operations as do
from matplotlib.colors import LogNorm
import matplotlib as mpl

from matplotlib.backends.backend_qt5agg import FigureCanvasQT as FigureCanvas, NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure
try:
    #from ciao_contrib.cda.data import download_chandra_obsids
    import ciao_contrib
except ImportError or ModuleNotFoundError:
    print("Failed to import CIAO python scripts. \n CIAO must be running prior to starting this script!")
    sys.exit(ClusterPyError.ciao_not_running)

import cluster
import ciao

class Stage(IntEnum):
    zero = 0
    one = 1
    two = 2
    three = 3
    four = 4
    five = 5
    tmap = 6


def process_commandline_arguments(cluster_obj):
    print("Processing commandline arguments")
    args = get_arguments()
    cluster_obj = None
    if None not in [args.name, args.abundance, args.nH, args.redshift, args.obsids]:
        ciao.initialize_cluster(name=args.name, obsids=args.obsids, abundance=args.abundance,
                                redshift=args.redshift, nH=args.nH)
    if args.init_cluster:
        cluster_obj = cluster.ClusterObj()
        cluster_obj.initialize_cluster()
    if args.config_file != "":
        if io.file_exists(args.config_file):
            cluster_obj = cluster.read_cluster_data(args.config_file)
        else:
            config_file = cluster.get_cluster_config(args.config_file)
            if config_file:
                cluster_obj = cluster.read_cluster_data(config_file)
            else:
                print("Error finding cluster configuration file. Try passing the full path to the file.")
                return
        if args.find_sources:
            ciao.find_sources(cluster_obj, ecf=args.ecf, energy=args.energy)

    if args.cont:
        if cluster_obj is None:
            cluster_obj = config.sys_config.current_cluster
            if cluster_obj is None:
                print("Cannot find a current working cluster.")
                exit(-1)

        if args.parallel and cluster_obj.last_step_completed == '1':
            ciao.run_stage_2_parallel(cluster_obj, args)
        else:
            ciao.start_from_last(cluster_obj, args)

    return cluster_obj


def get_arguments():
    from multiprocessing import cpu_count as max_cpu
    help_str = """
    pypeline help string
    """

    prog = 'clusterpyxt'

    # logger.debug("Getting commandline arguments.")
    parser = argparse.ArgumentParser(description=help_str, prog=prog)
    parser.add_argument("--initialize_cluster", dest="init_cluster", action="store_true",
                        help="First step in the pypeline. Creates the directory and initial configuration")

    parser.add_argument("--config_file", "-c", dest="config_file", action="store", default="",
                        help="The path and filename of the configuration file.")

    # parser.add_argument("--find_chandra_data", "-f", dest="find_data", action='store_true',
    #                    help="Find chandra observations to download")

    parser.add_argument("--continue", dest='cont', action='store_true',
                        help='Continue where you left off')

    parser.add_argument("--nH", dest='nH', action='store', default=None)
    parser.add_argument("--redshift", '-z', dest='redshift', action='store', default=None)
    parser.add_argument('--cluster_name', '-n', dest='name', action='store', default=None)
    parser.add_argument('--obsids', dest='obsids', action='store', nargs='+', default=None)
    parser.add_argument('--abundance', dest='abundance', action='store', default=None)
    parser.add_argument('--find_sources', dest='find_sources', action='store_true', default=False)
    parser.add_argument('--ecf', dest='ecf', action='store', default=0.3)
    parser.add_argument('--energy', dest='energy', action='store', default=0.3)
    parser.add_argument('--parallel', dest='parallel', action='store_true', default=False)
    parser.add_argument('--num_cpus', dest='num_cpus', action='store',type=int, default=max_cpu())
    parser.add_argument('--resolution', dest='resolution', action='store', default=2)

    args = parser.parse_args()

    # logger.debug("Finished getting commandline arguments")

    return args


class Stage3Window(QtWidgets.QMainWindow):

    def __init__(self, parent=None, cluster_obj=None):
        super(Stage3Window, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout()
        self.cluster = cluster_obj
        self.observations = self.cluster.observations


        
        
        #if False in acis_region_files_found:

        
        # self.obsid_label = QtWidgets.QLabel(self.observations[0].id)
        # self.canvas = FigureCanvas(Figure(figsize=(5,5)))
        
        # self.addToolBar(NavigationToolbar(self.canvas, self))

        # self.load_observation(self.observations[0])

        # layout.addWidget(self.obsid_label)
        # layout.addWidget(self.canvas)

        region_file_label = QtWidgets.QLabel(f"A region file (e.g. {self.observations[0].acisI_region_0_filename}) containing\n"
        "a small circular region covering each of the observation CCD chips is necessary to properly characterize the CCD.\n"
        "While any size upto the full CCD may be used, a region larger than 60 arc seconds is generally not necessary.")

        obs_string = self.get_obs_string()
        self.obsid_label = QtWidgets.QLabel(f"{obs_string}")
        
        layout.addWidget(region_file_label)
        layout.addWidget(self.obsid_label)

        self.run_stage_3_button = QtWidgets.QPushButton("Run Stage 3", self)
        self.run_stage_3_button.clicked.connect(self.run_stage_3)
        layout.addWidget(self.run_stage_3_button)
        
        widget = QtWidgets.QWidget(self)
        widget.setLayout(layout)

        self.setCentralWidget(widget)  
        

    def run_stage_3(self):
        if self.region_files_found():
            ciao.run_stage_3(self.cluster)
            self.cluster.last_step_completed = Stage.three.value
            ciao.finish_stage_3(self.cluster)

        else:
            print("Missing region files. Please see the documentation for instructions on how to create the region files.")

    def region_files_found(self):
        for obs in self.observations:
            if not io.file_exists(obs.acisI_region_0_filename):
                return False
        return True      

        

    def get_obs_string(self):
        obs_string_list = []
        for observation in self.observations:
            region_file = io.file_exists(observation.acisI_region_0_filename)
            obs_string_list.append(f'{observation.id}: {region_file}')
        
        return ", ".join(obs_string_list)

    def load_observation(self, observation):
        obs = observation
        print(f"Observation: {obs.clean}")
        self.ax = self.canvas.figure.subplots()
        data = obs.broad_flux #io.get_pixel_values(obs.broad_flux)
        

        cmap = mpl.cm.inferno
        cmap.set_bad(color='k', alpha=None)

        self.ax.imshow(data, norm=LogNorm(), cmap=cmap, origin='lower')
        self.ax.figure.canvas.draw()
 

class ProductMakingWindow(QtWidgets.QMainWindow):

    def __init__(self, parent=None, cluster_obj=None):
        super(ProductMakingWindow, self).__init__(parent)
        layout = QtWidgets.QVBoxLayout()

        self.cluster = cluster_obj

        self.temperature_map_button = QtWidgets.QPushButton("Make Temperature Map", self)
        self.temperature_map_button.clicked.connect(self.temperature_map_button_clicked)

        self.smoothed_xray_sb_button = QtWidgets.QPushButton("Make Smoothed X-ray SB Map", self)
        self.smoothed_xray_sb_button.clicked.connect(self.smoothed_xray_button_clicked)

        self.pressure_map_button = QtWidgets.QPushButton("Make Pressure Map", self)
        self.pressure_map_button.clicked.connect(self.pressure_button_clicked)
        self.pressure_map_button.setEnabled(False)

        self.entropy_map_button = QtWidgets.QPushButton("Make Entropy Map", self)
        self.entropy_map_button.clicked.connect(self.entropy_map_button_clicked)
        self.entropy_map_button.setEnabled(False)

        widgets = [self.temperature_map_button, self.smoothed_xray_sb_button, self.pressure_map_button, self.entropy_map_button]
        for widget in widgets:
            layout.addWidget(widget)

        widget = QtWidgets.QWidget(self)
        widget.setLayout(layout)

        self.setCentralWidget(widget)
        self.set_enabled()


    def smoothed_xray_button_clicked(self):
        self.disable_buttons()
        print("Making smoothed X-ray map.")
        acb.make_smoothed_xray_map(self.cluster)
        self.set_enabled()
        print("Done")
    
    def pressure_button_clicked(self):
        self.disable_buttons()
        print("Making pressure map")
        acb.make_pressure_map(self.cluster)
        self.set_enabled()
        print("Done")

    def temperature_map_button_clicked(self):
        self.disable_buttons()
        print("Making temperature map")
        acb.make_temperature_map(self.cluster, resolution=2)
        self.set_enabled()
        print("Done")

    def entropy_map_button_clicked(self):
        self.disable_buttons()
        print("Making entropy map")
        acb.make_entropy_map(self.cluster)    
        self.set_enabled()
        print("Done")

    def disable_buttons(self):
        buttons = [self.temperature_map_button, self.smoothed_xray_sb_button, self.pressure_map_button, self.entropy_map_button]
        for button in buttons:
            button.setEnabled(False)

    def set_enabled(self):
        self.temperature_map_button.setEnabled(True)
        self.smoothed_xray_sb_button.setEnabled(True)
        if io.file_exists(self.cluster.smoothed_xray_sb_cropped_nosrc_filename):
            self.pressure_map_button.setEnabled(True)
            self.entropy_map_button.setEnabled(True)
    
def alert_box(message_string, parent=None):
    alert = QtWidgets.QMessageBox(parent)
    alert.setIcon(QtWidgets.QMessageBox.Information)
    alert.setText(message_string)
    alert.setStandardButtons(QtWidgets.QMessageBox.Ok)
    return alert.exec_()

class ClusterWindow(QtWidgets.QMainWindow):

    def __init__(self, name=None, parent=None):
        super(ClusterWindow, self).__init__(parent)
        self.setWindowTitle(name)
        self.windows = list()
        layout = QtWidgets.QVBoxLayout()

        name_label = QtWidgets.QLabel('Cluster Name', self)
        obsid_label = QtWidgets.QLabel('Observation IDs', self)
        nH_label = QtWidgets.QLabel('Hydrogen Column Density', self)
        redshift_label = QtWidgets.QLabel('Redshift', self)
        abundance_label = QtWidgets.QLabel('Metallicity', self)
        signal_to_noise_label = QtWidgets.QLabel('Signal to Noise Ratio Desired', self)
        
        self.initialized = False

        if name: # cluster already initialized
            self.initialized = True
            self._cluster_obj = cluster.load_cluster(name)
            self.cluster_name = self._cluster_obj.name
            self.obsids = " ".join(self._cluster_obj.observation_list)
            
            if float(self._cluster_obj.hydrogen_column_density) < 1:
                self._cluster_obj.hydrogen_column_density = \
                    str(float(self._cluster_obj.hydrogen_column_density) * 1e22)
            
            self.hydrogen_column_density = self._cluster_obj.hydrogen_column_density
            self.abundance = self._cluster_obj.abundance
            self.redshift = self._cluster_obj.redshift
            self.signal_to_noise = self._cluster_obj.signal_to_noise
            self.last_step_completed = int(self._cluster_obj.last_step_completed)

        else: # cluster needs to be initialized
            self.cluster_name = ""
            self.obsids = ""
            self.hydrogen_column_density = ""
            self.abundance = "0.3"
            self.redshift = ""
            self.signal_to_noise = "50"
            save_update_text = 'Save Cluster Attributes'
            self.last_step_completed = 0

            obsid_label.setText("Observation IDs (seperated by a space)")
            nH_label.setText("Hydrogen Column Density (e.g. 5.24e22)")

        self.name_text = QtWidgets.QLineEdit(self.cluster_name, self)
        self.obsid_text = QtWidgets.QPlainTextEdit(self.obsids, self)       
        self.nH_text = QtWidgets.QLineEdit(self.hydrogen_column_density, self)
        self.redshift_text = QtWidgets.QLineEdit(self.redshift, self)
        self.abundance_text = QtWidgets.QLineEdit(self.abundance, self)
        self.signal_to_noise_text = QtWidgets.QLineEdit(self.signal_to_noise, self)
        
        if not self.initialized:
            save_update_button = QtWidgets.QPushButton(save_update_text, self)
            save_update_button.clicked.connect(self.save_update_button_clicked)
        else:
            cluster_attributes = [self.name_text, self.obsid_text, self.nH_text, 
                                self.redshift_text, self.abundance_text, self.signal_to_noise_text]
            for attr in cluster_attributes:
                attr.setEnabled(False)
        
        self.stage_1_button = QtWidgets.QPushButton('Run Stage 1', self)
        self.stage_1_button.clicked.connect(self.run_stage_1)
        
        
        self.stage_2_files = QtWidgets.QLabel('Sources file:\nExclude file:')
        self.stage_2_button = QtWidgets.QPushButton('Run Stage 2', self)
        self.stage_2_button.clicked.connect(self.run_stage_2)
        self.update_stage_2_text()
        
        self.stage_3_button = QtWidgets.QPushButton('Run Stage 3', self)
        self.stage_3_button.clicked.connect(self.run_stage_3)
        
        self.stage_4_button = QtWidgets.QPushButton('Run Stage 4', self)
        self.stage_4_button.clicked.connect(self.run_stage_4)
        
        self.stage_5_button = QtWidgets.QPushButton('Run Stage 5', self)
        self.stage_5_button.clicked.connect(self.run_stage_5)
        
        self.spectral_fitting_button = QtWidgets.QPushButton('Spectral Fitting', self)
        self.spectral_fitting_button.clicked.connect(self.run_spectral_fits)
        
        self.products_button = QtWidgets.QPushButton('Make Final Products', self)
        self.products_button.clicked.connect(self.make_products_clicked)


        self.buttons = [self.stage_1_button, self.stage_2_button, self.stage_3_button, self.stage_4_button, 
                        self.stage_5_button, self.spectral_fitting_button, self.products_button]

        stages = [
                self.stage_1_button, 
                self.stage_2_files, self.stage_2_button, 
                self.stage_3_button, 
                self.stage_4_button, 
                self.stage_5_button, 
                self.spectral_fitting_button,
                self.products_button
                ]

        for button in self.buttons:
            button.setEnabled(False)


        if not self.initialized:
            stages.insert(0, save_update_button)

        widgets = [name_label, self.name_text, 
                    obsid_label, self.obsid_text, 
                    nH_label, self.nH_text, 
                    redshift_label, self.redshift_text, 
                    abundance_label, self.abundance_text, 
                    signal_to_noise_label, self.signal_to_noise_text] + stages

        for widget in widgets:
            layout.addWidget(widget)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)

        self.setCentralWidget(widget)
        self.update_stages()


    # def closeEvent(self, event):
    #     self.parent.load_cluster_list()

    def save_update_button_clicked(self):
        cluster_name = self.name_text.text()
        obsids = self.obsid_text.toPlainText().split(' ')
        hydrogen_column_density = self.nH_text.text()
        redshift = self.redshift_text.text()
        abundance = self.abundance_text.text()
        signal_to_noise = self.signal_to_noise_text.text()
        
        cluster_obj = cluster.ClusterObj(name=cluster_name,
                                    observation_ids=obsids,
                                    data_directory=config.data_directory(),
                                    hydrogen_column_density=hydrogen_column_density,
                                    abundance=abundance,
                                    redshift=redshift,
                                    last_step_completed=0,
                                    signal_to_noise=signal_to_noise)
        io.make_directory(cluster_obj.directory)
        cluster_obj.write_cluster_data()
        io.make_initial_directories(cluster_obj)

        alert = QtWidgets.QMessageBox(self)
        alert.setIcon(QtWidgets.QMessageBox.Information)
        alert.setText(f"{cluster_name} initialized. Please restart ClusterPyXT to continue.")
        alert.setStandardButtons(QtWidgets.QMessageBox.Ok)
        return alert.exec_()

    def update_stages(self):
        self.update_stage_2_text()
        self.update_buttons()
        

    def update_buttons(self):
        self.last_step_completed = int(self._cluster_obj.last_step_completed)
        if self.initialized:
            for i in range(self.last_step_completed+1):
                self.buttons[i].setEnabled(True)
        
            if io.file_exists(self._cluster_obj.spec_fits_file):
                if io.num_lines_in(self._cluster_obj.spec_fits_file) > 10:
                    self.products_button.setEnabled(True)


    def run_stage_1(self):
        ciao.run_stage_1(self._cluster_obj)
        self._cluster_obj.last_step_completed = Stage.one.value
        ciao.finish_stage_1(self._cluster_obj)
        #self.display_stage_2_message()
        self.update_stage_2_text()

    def run_stage_2(self):
        if self.stage_2_files_exist:
            ciao.run_stage_2_parallel(self._cluster_obj, get_arguments())
            self._cluster_obj.last_step_completed = Stage.two.value
            ciao.finish_stage_2(self._cluster_obj)
        else:
            self.display_stage_2_message()
        return

    def run_stage_3(self):
        args = get_arguments()
        # ciao.run_stage_3(self._cluster_obj, args.num_cpus)
        # self._cluster_obj.last_step_completed = Stage.three.value
        # ciao.finish_stage_3(self._cluster_obj)
        win = Stage3Window(self, self._cluster_obj)
        self.windows.append(win)
        win.show()
        

    def run_stage_4(self):
        ciao.run_stage_4(self._cluster_obj)
        self._cluster_obj.last_step_completed = Stage.four.value
        ciao.finish_stage_4(self._cluster_obj)
        return

    def run_stage_5(self):
        args = get_arguments()
        ciao.run_stage_5(self._cluster_obj, args)
        self._cluster_obj.last_step_completed = Stage.five.value
        ciao.finish_stage_5(self._cluster_obj)
        return

    def run_spectral_fits(self):
        ciao.print_stage_tmap_prep(self._cluster_obj)

    def make_products_clicked(self):
        win = ProductMakingWindow(parent=self, cluster_obj=self._cluster_obj)
        self.windows.append(win)
        win.show()

    @property
    def stage_2_files_exist(self):
        sources_file_exists = False
        exclude_file_exists = False
        if hasattr(self, '_cluster_obj'):
            sources_file_exists = io.file_exists(self._cluster_obj.sources_file)
            exclude_file_exists = io.file_exists(self._cluster_obj.exclude_file)
        return sources_file_exists and exclude_file_exists

    def update_stage_2_text(self):
        if hasattr(self, '_cluster_obj'):
            sources_file_exists = io.file_exists(self._cluster_obj.sources_file)
            exclude_file_exists = io.file_exists(self._cluster_obj.exclude_file)
            
            text = f"Sources file: {sources_file_exists}\n Exclude file: {exclude_file_exists}"
            self.stage_2_files.setText(text)


    def update_stage_3_button_text(self):
        pass

    def display_stage_2_message(self):
        sb_map_filename = self._cluster_obj.xray_surface_brightness_filename
        sources_file = self._cluster_obj.sources_file
        exclude_file = self._cluster_obj.exclude_file
        cluster_name = self._cluster_obj.name
        message_string = f"""Data downloaded and the observations are merged into a surface brightness map {sb_map_filename}. 
    Now it is time to filter out point sources and high energy flares. To do so, first open the surface brightness
    map and create regions around sources you want excluded from the data analysis. These are typically foreground
    point sources one does not want to consider when analyzing the cluster. Save these regions as a DS9 region file
    named {sources_file}. 

    Additionally, you need to create a region file containing any regions you wanted excluded from the deflaring process.
    This would include areas such as the peak of cluster emission as these regions may contain high energy events we want
    to consider in this analysis. Save this region file as {exclude_file}. 

    After both files are saved, you can continue ClusterPyXT on {cluster_name}"""
        alert_box(message_string)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle('ClusterPyXT')
        self.setMinimumSize(500,400)
        self._selected_cluster_name = ""
        self._cluster_configs = config.get_cluster_configs()

        layout = QtWidgets.QVBoxLayout()
        
        label = QtWidgets.QLabel('Cluster List', self)
        
        self.cluster_list = QtWidgets.QListWidget(self)

        self.cluster_list.itemClicked.connect(self.list_selection_changed)
        self.cluster_list.itemDoubleClicked.connect(self.button_clicked)
        self.load_cluster_list()

        self.continue_cluster_button = QtWidgets.QPushButton('Continue cluster')
        self.continue_cluster_button.clicked.connect(self.button_clicked)
        self.continue_cluster_button.setEnabled(False)

        init_cluster_button = QtWidgets.QPushButton('New Cluster')
        init_cluster_button.clicked.connect(self.new_cluster_clicked)
        
        widgets = [label, self.cluster_list, self.continue_cluster_button, init_cluster_button]

        for widget in widgets:
            layout.addWidget(widget)
        
        widget = QtWidgets.QWidget()
        widget.setLayout(layout)

        self.setCentralWidget(widget)

        self.windows = list()

        if len(self._cluster_configs) == 0:
            window = ClusterWindow(parent=self)
            self.windows.append(window)
            window.show()

    def button_clicked(self):
        cluster_name = self._selected_cluster_name 
        window = ClusterWindow(name=cluster_name, parent=self)
        self.windows.append(window)
        window.show()

    def new_cluster_clicked(self):
        window = ClusterWindow(name=None, parent=self)
        self.windows.append(window)
        window.show()

    def list_selection_changed(self, item):
        try:
            selection_index = self.cluster_list.selectedIndexes()[0].row()
            cluster_name = self._cluster_configs[selection_index][0] # index 0 = name, index 1 = confi filename
            self._selected_cluster_name = cluster_name
            self.continue_cluster_button.setEnabled(True)
        except IndexError:
            self._selected_cluster_name = None
            self.continue_cluster_button.setEnabled(False)

    def load_cluster_list(self):
        self._cluster_configs = config.get_cluster_configs()
        for cluster_config in self._cluster_configs:
            self.cluster_list.addItem(cluster_config[0])

if __name__ == "__main__":
    if 1 == len(sys.argv):
        app = QtWidgets.QApplication([])
        if config.sys_config.data_directory in ["", "."]:
            config.sys_config.get_data_dir()
        if config.sys_config.ciao_directory in ["", "."]:
            config.sys_config.get_ciao_dir()
        win = MainWindow()
        win.show()
        app.exec_()

    else:
        cluster_obj = cluster.ClusterObj()
        args = process_commandline_arguments(cluster_obj)
