# version 0.9.0a
# Jan 14, 2020

import sys
import argparse
import pypeline_io as io
import cluster
from errors import ClusterPyError
import config
from PyQt5 import QtCore, QtGui, QtWidgets

import acb

config.initialize_pypeline()

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
            cluster_obj = config.current_cluster()
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

class ASmoothWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None, cluster_obj=None):
        super(ASmoothWindow, self).__init__(parent)

        self.cluster = cluster_obj

        self.min_label = QLabel('Minimum Smoothing Radius (pixels):', self)
        self.min_input = QLineEdit("", self)
        self.max_label = QLabel('Maximum Smoothing Radius (pixels):', self)
        self.max_input = QLineEdit("", self)
        self.num_iter_label = QLabel('Number of Iterations:',self)
        self.num_iter_input = QLineEdit("", self)
        self.counts_label = QLabel('Minimum Counts', self)
        self.counts_input = QLineEdit("", self)


        self.smooth_image_button = QPushButton("Adaptively Smooth Image")
        self.smooth_image_button.clicked.connect(self.smooth_image)
        layout = QVBoxLayout()

        layout.addWidget(self.canvas)
        layout.addWidget(self.smooth_file_button)
        layout.addWidget(self.min_label)
        layout.addWidget(self.min_input)
        layout.addWidget(self.max_label)
        layout.addWidget(self.max_input)
        layout.addWidget(self.num_iter_label)
        layout.addWidget(self.num_iter_input)
        layout.addWidget(self.counts_label)
        layout.addWidget(self.counts_input)
        layout.addWidget(self.smooth_image_button)

        self.setLayout(layout)

    def smooth_image(self):
        from ciao_contrib import runtool as rt
        infile = self.cluster.smoothed_xray_sb_cropped_nosrc_filename
        min_rad = float(self.min_input.text())
        max_rad = float(self.max_input.text())
        iterations = int(self.num_iter_input.text())
        counts = int(self.counts_input.text())
        directory = self.cluster.output_dir
        cluster = self.cluster.name

        outfile = "{dir}/{cluster}-smooth-min_{min:0.2f}-max_{max:0.2f}-iter_{iter}-counts_{counts}.fits".format(
            dir=directory,
            cluster=cluster,
            min=min_rad,
            max=max_rad,
            iter=iterations,
            counts=counts)

        print("Generating {}".format(outfile))
        rt.dmimgadapt(infile=infile, outfile=outfile, function='gaussian', minrad=min_rad, maxrad=max_rad, num=iterations,
                      counts=counts, radscale='log')


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
        
        initialized = False

        if name: # cluster already initialized
            initialized = True
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
        
        if not initialized:
            save_update_button = QtWidgets.QPushButton(save_update_text, self)
            save_update_button.clicked.connect(self.save_update_button_clicked)
        else:
            cluster_attributes = [self.name_text, self.obsid_text, self.nH_text, 
                                self.redshift_text, self.abundance_text, self.signal_to_noise_text]
            for attr in cluster_attributes:
                attr.setEnabled(False)
        
        stage_1_button = QtWidgets.QPushButton('Run Stage 1', self)
        stage_1_button.clicked.connect(self.run_stage_1)
        
        stage_2_button = QtWidgets.QPushButton('Run Stage 2', self)
        stage_2_button.clicked.connect(self.run_stage_2)
        
        stage_3_button = QtWidgets.QPushButton('Run Stage 3', self)
        stage_3_button.clicked.connect(self.run_stage_3)
        
        stage_4_button = QtWidgets.QPushButton('Run Stage 4', self)
        stage_4_button.clicked.connect(self.run_stage_4)
        
        stage_5_button = QtWidgets.QPushButton('Run Stage 5', self)
        stage_5_button.clicked.connect(self.run_stage_5)
        
        spectral_fitting_button = QtWidgets.QPushButton('Spectral Fitting', self)
        spectral_fitting_button.clicked.connect(self.run_spectral_fits)
        
        products_button = QtWidgets.QPushButton('Make Final Products', self)
        products_button.clicked.connect(self.make_products_clicked)

        stage_buttons = [stage_1_button, stage_2_button, stage_3_button, 
                        stage_4_button, stage_5_button, spectral_fitting_button, 
                        products_button]

        for button in stage_buttons:
            button.setEnabled(False)

        if initialized:
            for i in range(self.last_step_completed+1):
                stage_buttons[i].setEnabled(True)
        
            if io.file_exists(self._cluster_obj.spec_fits_file):
                if io.num_lines_in(self._cluster_obj.spec_fits_file) > 10:
                    products_button.setEnabled(True)

        if not initialized:
            stage_buttons.insert(0, save_update_button)

        widgets = [name_label, self.name_text, 
                    obsid_label, self.obsid_text, 
                    nH_label, self.nH_text, 
                    redshift_label, self.redshift_text, 
                    abundance_label, self.abundance_text, 
                    signal_to_noise_label, self.signal_to_noise_text] + stage_buttons

        for widget in widgets:
            layout.addWidget(widget)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)

        self.setCentralWidget(widget)

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
        retval = alert.exec_()

    def run_stage_1(self):
        ciao.run_stage_1(self._cluster_obj)
        self._cluster_obj.last_step_completed = Stage.one.value
        ciao.finish_stage_1(self._cluster_obj)

    def run_stage_2_parallel(self):
        ciao.run_stage_2_parallel(self._cluster_obj, get_arguments())
        self._cluster_obj.last_step_completed = Stage.two.value
        ciao.finish_stage_2(cluster)
        return

    def run_stage_3(self):
        args = get_arguments()
        ciao.run_stage_3(self._cluster_obj, args.num_cpus)
        self._cluster_obj.last_step_completed = Stage.three.value
        ciao.finish_stage_3(self._cluster_obj)
        return

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
        for cluster_config in self._cluster_configs:
            self.cluster_list.addItem(cluster_config[0])
        self.cluster_list.itemClicked.connect(self.list_selection_changed)
        self.cluster_list.itemDoubleClicked.connect(self.buttonClicked)

        continue_cluster_button = QtWidgets.QPushButton('Continue cluster')
        continue_cluster_button.clicked.connect(self.buttonClicked)

        init_cluster_button = QtWidgets.QPushButton('New Cluster')
        init_cluster_button.clicked.connect(self.buttonClicked)
        
        widgets = [label, self.cluster_list, continue_cluster_button, init_cluster_button]

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

    def buttonClicked(self):
        cluster_name = self._selected_cluster_name 
        window = ClusterWindow(name=cluster_name, parent=self)
        self.windows.append(window)
        window.show()

    def list_selection_changed(self, item):
        selection_index = self.cluster_list.selectedIndexes()[0].row()
        cluster_name = self._cluster_configs[selection_index][0] # index 0 = name, index 1 = confi filename
        self._selected_cluster_name = cluster_name 


if __name__ == "__main__":
    if 1 == len(sys.argv):
    
        cluster_configs = config.get_cluster_configs()
        app = QtWidgets.QApplication([])
        win = MainWindow()
        win.show()
        app.exec_()

    else:
        cluster_obj = cluster.ClusterObj()
        args = process_commandline_arguments(cluster_obj)
