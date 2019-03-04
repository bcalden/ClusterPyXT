# ClusterPyXT
## Introduction
ClusterPyXT is a software pipeline to automate the creation of x-ray temeprature maps, pressure maps, surface brightness maps, and density maps. It is open source and under active development. Please feel free to contribute!

## Requirements
This version of `ClusterPyXT` requires `CIAO-4.9` or later (4.11 recommended). The full calibration database (CALDB) is a requirement as well and can be installed with CIAO. 

### CIAO Installation
These instructions are for `CIAO 4.11`. Follow the installation instructions at the [Chandra X-ray Center (CXC)](http://cxc.harvard.edu/ciao/download/). Note, the custom installation option should be used as it allows for the full `CALDB` installation. Make sure to select all `CALDB` options before downloading the installation script.

Another requirement for `ClusterPyXT` is the `astropy` python library within the `CIAO` environment. `CIAO 4.11` allows for the easy installation of this library. 
After installation, start the `CIAO` environment and run `pip3 install astropy`. 

### Download ClusterPyXT
To download ClusterPyXT, simply run `git clone https://github.com/bcalden/ClusterPyXT.git`.

## Running ClusterPyXT
### System Configuration
After following the instructions above, go the `ClusterPyXT` directory and run `python clusterpyxt.py` to initialize the system configuration. 

### Cluster Initialization
Next, you must initialize a cluster. At a minimum, you need the name for the cluster and the Chandra Observation IDs you plan on using. Names can be whatever you want (e.g. A85, A115, Bullet, Toothbrush) just recognize the name needs to be valid in directory and file names. That means don't use slashes or other characters disallowed by your filesystem. Observation IDs can be found using the [Chandra Data Archive](https://cda.harvard.edu/chaser/). While you can start the pipeline with just this information, redshift, hydrogen column density, and the metallicity of the cluster are required to complete the spectral fitting. Redshift information can be found at the [NASA/IPAC Extragalactic Database (NED)](https://ned.ipac.caltech.edu). Hydrogen column density information can be found at [NASA HEASARC Tools](https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl). Solar abundance can usually be estimated at `0.3` although you can check the literature to see if a better value for your cluster should be used.

There are various ways to initialize a cluster. For a textual guided process, you can run `python clusterpyxt.py --initialize_cluster`. If you prefer to pass all information via command line arguments (can be useful to initialize multiple clusters within a script), the `--cluster_name`, `--nH`, `--redshift`, `--abundance`, and `--obsids` arguments can be used. 

### Processing a Cluster
After a cluster is initialized, you can simply start/continue the pipeline by running `python clusterpyxt.py --continue`. The pipeline will prompt for additional input (point source region files, exclusion areas, etc.). 

### Spectral Fitting
As there can be anywhere from 10<sup>3</sup> to 10<sup>5</sup> regions to fit, there are multiple ways to create a temperature map. One can either process every region in serial on their local computer, run it in parallel on a single, multicore computer, or even in parallel on a supercomputer. As scheduling on a supercomputer is highly specific to each environment, only a general description is provided. Feel free to contact us for help.

To do the spectral fitting, only a subset of the data is required. If processing on a remote machine, insure the remote machine has the required software and `ClusterPyXT` is configured (see above). Make a directory for your cluster in the remote cluster data directory (set in the system configuration on first run, also set in `ClusterPyXT\pypeline_config.ini`). The files required are the configuration file (`ClusterName_pypeline_config.ini`) and the `acb` folder within the cluster directory. Upload both of these to the remote machines cluster folder you just created. You are now ready for spectral fitting.

To run the spectral fitting portion of the pipeline in serial, run `python wrapper_serial.py --cluster_config_file 'path/to/cluster_config_file'`.

To run in parallel, run `python wrapper.py --parallel --num_cpus N --cluster_config_file 'path/to/cluster\_config\_file'`. Note, in the current version of the pipeline, there is a bug which can cause a memory leak when running in parallel on a single machine. If you need to restart for this reason, or any reason, simply add the `--continue` argument to the above `wrapper.py` command and `ClusterPyXT` will begin where it left off without having to re-fit any region.

To run on a supercomputer, you can make use of the command file generated (`commands_ClusterName.lis`) in the `acb` directory. This command file has a line for each region to be fit that directly calls the spectral fitting routine on that region. You can write a simple script to parse this command file and send it to each of the nodes used in the supercomputer. 

### Temperature Map Creation
After spectral fitting, the last thing to do is create the temperature map. If you did the spectral fitting on a remote machine, you need to download the three `.csv` files created within the remote `acb` directory. Next, simply run `python acb.py --temperature_map --cluster_config_file 'path/to/cluster_config_file'`.
Check the `clustername/main_output/` directory for the output.

### Pressure Map Creation
After generating the temperature map, there is now enough data to generate the pseudo pressure map. Simply run `python acb.py --make_pressure_map --cluster_config_file 'path/to/config/file'`. Check the `clustername/main_output` directory for the output.
