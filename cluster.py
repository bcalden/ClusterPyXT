import os
import sys
import configparser
import csv
import config
import pypeline_io as io
import numpy as np
import astropy.io.fits as fits # migrate pycrates commands to fits

acisI = [0, 1, 2, 3]
acisS = [4, 5, 6, 7, 8, 9]
back_illuminated_ids = [5, 7]


class CCD:
    def __init__(self, observation=None, id=None):
        self.observation = observation
        self.id = id
        if self.id in acisI:
            self.type = "I"
        elif self.id in acisS:
            self.type = "S"
            print("")
        else:
            print("Error with observation {} chip id {}. Cannot determine type (ACIS-I or ACIS-S).".format(
                observation.id, self.id))
            exit(-1)
        if self.id in back_illuminated_ids:
            self.back_illuminated = True
        else:
            self.back_illuminated = False

    def __repr__(self):
        return "OBSID {}:ACIS-{}:CCD {}: Back Illuminated={}".format(self.observation.id, self.type, self.id,
                                                                     self.back_illuminated)


class Observation:
    def __init__(self,
                 obsid=0,
                 cluster=None,
                 ):
        self.id = obsid
        self.cluster = cluster
        self.set_ccds()

    def set_ccds(self):
        try:
            self.ccds = self.get_ccds()
            self.acis_I_chips, self.acis_S_chips = self.get_acis_I_and_S_chips()

            if len(self.acis_S_chips) > 0:
                print("Observation {} has ACIS-S data that is not yet supported. "
                      "Feel free to implement and submit a pull request!".format(self.id))
                if len(self.acis_I_chips) > 0:
                    print("ACIS-I chips to be used:")
                    for ccd in self.acis_I_chips:
                        print(ccd)
        except FileNotFoundError:
            #print("Observation {} not yet downloaded".format(self.id))
            pass

    def get_ccds(self):
        chip_ids = self.oif_detnam

        ccds = []
        if chip_ids[:5].upper() == "ACIS-":
            ids = chip_ids[5:]

            for chip in ids:
                ccds.append(CCD(observation=self, id=int(chip)))
        else:
            print("Error determining ACIS CCDs used. Check observation {} and try again".format(self.id))
            exit(-1)

        return ccds

    def get_acis_I_and_S_chips(self):
        acis_I_ccds = []
        acis_S_ccds = []
        for ccd in self.ccds:
            if ccd.type == 'I':
                acis_I_ccds.append(ccd)
            else:
                acis_S_ccds.append(ccd)
        return acis_I_ccds, acis_S_ccds

    @property
    def directory(self):
        return io.get_path("{cluster_dir}/{obsid}/".format(cluster_dir=self.cluster.directory,
                                                           obsid=self.id))

    @property  # refactor so exclude points to the content of the file, not the file
    def exclude(self):
        return self.cluster.exclude_file

    @property
    def exclude_file(self):
        return self.cluster.exclude_file

    @property
    def oif_filename(self):
        return io.get_path("{obs_dir}/oif.fits".format(obs_dir=self.directory))

    @property
    def oif_fits(self):
        return fits.open(self.oif_filename)

    @property
    def oif_detnam(self):
        return self.oif_fits[1].header['detnam']

    @property
    def response_file_region_covering_ccds(self):
        return io.get_path("{analysis_dir}/acisI_region_0.reg".format(
            analysis_dir=self.analysis_directory
        ))

    @property
    def analysis_directory(self):
        return io.get_path("{obs_dir}/analysis/".format(obs_dir=self.directory))

    @property
    def reprocessing_directory(self):
        return io.get_path("{analysis_dir}/repro/".format(analysis_dir=self.analysis_directory))
        #return io.get_path("{obs_dir}/repro/".format(obs_dir=self.directory))

    @property
    def primary_directory(self):
        return io.get_path("{observation_dir}/primary/".format(observation_dir=self.directory))

    @property
    def secondary_directory(self):
        return io.get_path("{observation_dir}/secondary/".format(observation_dir=self.directory))

    @property
    def combined_directory(self):
        return self.cluster.combined_directory

    @property
    def global_response_directory(self):
        return io.get_path('{analysis_directory}/globalresponse/'.format(analysis_directory=self.analysis_directory))

    @property
    def data_filename(self):
        return io.get_path('{analysis_dir}/acisI.fits'.format(analysis_dir=self.analysis_directory))

    @property
    def back_filename(self):
        return io.get_path('{analysis_dir}/merged_back.fits'.format(analysis_dir=self.analysis_directory))

    @property
    def clean(self):
        return io.get_path("{analysis_dir}/acisI_clean.fits".format(analysis_dir=self.analysis_directory))

    @property
    def sc_clean(self):
        return io.get_path("{super_comp_dir}/acisI_clean_{obsid}.fits".format(
            super_comp_dir=self.cluster.super_comp_dir,
            obsid=self.id
        ))

    @property
    def back(self):
        return self.background_nosrc_filename
        #return io.get_path("{analysis_dir}/backI_clean.fits".format(analysis_dir=self.analysis_directory))

    @property
    def sc_back(self):
        return io.get_path("{super_comp_dir}/backI_clean_{obsid}.fits".format(
            super_comp_dir=self.cluster.super_comp_dir,
            obsid=self.id
        ))


    @property
    def fov_file(self):
        return io.get_filename_matching("{}/*{}*fov1.fits".format(
            self.analysis_directory,
            self.id
        ))[0]


    @property
    def acisI_comb_img(self):
        return io.get_path("{combined_dir}/acisI_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id))


    @property
    def backI_comb_img(self):
        return io.get_path("{combined_dir}/backI_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id)
        )

    @property
    def effbtime(self):
        return io.get_path("{acb_dir}/effbtime-{obsid}_circle.dat".format(
            acb_dir=self.cluster.acb_dir,
            obsid=self.id
        ))

    @property
    def effdtime(self):
        return io.get_path("{acb_dir}/effdtime-{obsid}_circle.dat".format(
            acb_dir=self.cluster.acb_dir,
            obsid=self.id
        ))

    @property
    def merged_back_lis(self):
        return io.get_path("{analysis_dir}/merged_back.lis".format(analysis_dir=self.analysis_directory))

    @property
    def acisI_combined_mask(self):
        return io.get_path("{combined_dir}/acisI_comb_mask-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def acisI_nosrc_combined_mask(self):
        return io.get_path("{combined_dir}/acisI_comb_mask_nosrc-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def acis_nosrc_filename(self):
        return io.get_path("{analysis_dir}/acis_nosrc_{obsid}.fits".format(
            analysis_dir=self.analysis_directory,
            obsid=self.id
        ))

    @property
    def background_nosrc_filename(self):
        return io.get_path("{analysis_dir}/back_nosrc_{obsid}.fits".format(
            analysis_dir=self.analysis_directory,
            obsid=self.id
        ))

    @property
    def temp_acis_comb_mask_filename(self):
        return io.get_path("{combined_dir}/acis_comb_mask_int-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def acisI_high_energy_combined_image(self):
        return io.get_path("{combined_dir}/acisI_hien_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def acisI_high_energy_temp_image(self):
        return io.get_path("{combined_dir}/acisI_high_energy_int_img.fits".format(
            combined_dir=self.combined_directory
        ))

    @property
    def backI_high_energy_combined_image(self):
        return io.get_path("{combined_dir}/backI_hien_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def backI_high_energy_temp_image(self):
        return io.get_path("{combined_dir}/backI_hien_int_img.fits".format(
            combined_dir=self.combined_directory
        ))

    @property
    def scale_map_region_list_filename(self):
        return "{acb_dir}/{cluster_name}_{obsid}_scale_map_region_list.reg".format(
            acb_dir=self.cluster.acb_dir,
            cluster_name=self.cluster.name,
            obsid=self.id
        )

    @property
    def scale_map_region_list(self):
        region_list = []
        with open(self.scale_map_region_list_filename, 'r') as f:
            reader = csv.reader(f, delimiter='#')
            region_list = list(reader)

        return region_list


    @scale_map_region_list.setter
    def scale_map_region_list(self, circle_list):
        with open(self.scale_map_region_list_filename, 'w') as f:
            writer = csv.writer(f, delimiter='#')
            writer.writerows(circle_list)


    @property
    def bad_pixel_file(self):
        return io.get_path("{obs_analysis_dir}/bpix1_new.fits".format(
            obs_analysis_dir=self.analysis_directory
        ))

    @property
    def aux_response_file(self):
        return io.get_path("{analysis_directory}/globalresponse/acisI.arf".format(
            analysis_directory=self.analysis_directory
        ))

    @property
    def original_reprocessed_evt2_filename(self):
        evt2_filename = io.get_filename_matching("{}/acis*repro_evt2.fits".format(self.reprocessing_directory))
        if isinstance(evt2_filename, list):
            if len(evt2_filename) >= 1:
                evt2_filename = evt2_filename[-1]
                return io.get_path(evt2_filename)
        return None

    @property
    def original_reprocessed_evt2_file_exists(self):
        if self.original_reprocessed_evt2_filename:
            return True
        else:
            return False

    @property
    def original_reprocessed_bad_pixel_filename(self):
        bpix1_filename = io.get_filename_matching("{}/*repro_bpix1.fits".format(self.reprocessing_directory))
        if isinstance(bpix1_filename, list):
            if len(bpix1_filename) >= 1:
                bpix1_filename = bpix1_filename[-1]
                return io.get_path(bpix1_filename)
        return None

    @property
    def original_reprocessed_bad_pixel_file_exists(self):
        if self.original_reprocessed_bad_pixel_filename:
            return True
        else:
            return False

    @property
    def level_1_event_filename(self):
        return io.get_filename_matching("{analysis_dir}/acis*evt1.fits".format(
            analysis_dir=self.analysis_directory
        ))[0]

    @property
    def reprocessed_evt2_filename(self):
        return self.original_reprocessed_evt2_filename
        #return io.get_path("{analysis_dir}/evt2.fits".format(analysis_dir=self.analysis_directory))

    @property
    def reprocessed_bad_pixel_filename(self):
        return self.original_reprocessed_bad_pixel_filename
        #return io.get_path("{analysis_dir}/bpix1_new.fits".format(analysis_dir=self.analysis_directory))

    @property
    def redistribution_matrix_file(self):
        return io.get_path("{analysis_directory}/globalresponse/acisI.rmf".format(
            analysis_directory=self.analysis_directory
        ))

    @property
    def rmf_sc(self):
        return io.get_path("{super_comp_dir}/{name}_{obsid}.rmf".format(
            super_comp_dir=self.cluster.super_comp_dir,
            name=self.cluster.name,
            obsid=self.id
        ))

    @property
    def arf_sc(self):
        return io.get_path("{super_comp_dir}/{name}_{obsid}.arf".format(
            super_comp_dir=self.cluster.super_comp_dir,
            name=self.cluster.name,
            obsid=self.id
        ))

    @property
    def exposure_time_file(self):
        return io.get_path("{super_comp_dir}/{name}_{obsid}_exptime.dat".format(
            super_comp_dir=self.cluster.super_comp_dir,
            name=self.cluster.name,
            obsid=self.id
        ))

    @property
    def exposure_time(self):
        return float(io.read_line_number(self.exposure_time_file, 1))

    @property
    def pcad_asol(self):
        return io.get_path("{analysis_dir}/pcad_asol1.lis".format(
            analysis_dir=self.analysis_directory
        ))

    @property
    def acis_mask(self):
        filename = io.get_filename_matching("{analysis_dir}/*_msk1.fits".format(
            analysis_dir=self.analysis_directory
        ))[0]
        return filename

    @property
    def acis_mask_sc(self):
        return io.get_path("{super_comp_dir}/{obsid}_msk1.fits".format(
            super_comp_dir=self.cluster.super_comp_dir,
            obsid=self.id
        ))

    @property
    def effective_data_time_file(self):
        return io.get_path('{acb_dir}/{name}_effective_data_time_{obs}.npy'.format(
            acb_dir=self.cluster.acb_dir,
            name=self.cluster.name,
            obs=self.id
        ))

    @property
    def effective_data_time(self):
        return np.load(self.effective_data_time_file)

    @effective_data_time.setter
    def effective_data_time(self, effective_time):
        np.save(self.effective_data_time_file, effective_time)

    @property
    def effective_background_time_file(self):
        return io.get_path('{acb_dir}/{name}_effective_background_time_{obs}.npy'.format(
            acb_dir=self.cluster.acb_dir,
            name=self.cluster.name,
            obs=self.id
        ))

    @property
    def effective_background_time(self):
        return np.load(self.effective_background_time_file)

    @effective_background_time.setter
    def effective_background_time(self, effective_time):
        np.save(self.effective_background_time_file, effective_time)

    @property
    def ccd_merge_list(self):
        return io.get_path('{obs_analysis_dir}/acisI.lis'.format(
            obs_analysis_dir=self.analysis_directory)
        )

    def coordinates_for_scale_map_region(self, region, scale_map_regions):
        return np.where(scale_map_regions == region)

    def coordinates_for_big_region_index(self, region, region_index_map):
        return np.where(region_index_map == region)

    def get_region_from_region_number(self, region_number):
        with open(self.scale_map_region_list_filename, 'r') as f:
            reader = csv.reader(f, delimiter='#')
            for row in reader:
                if int(row[1]) == int(region_number):
                    return row[0]
            print('Error: region {reg} not found. exiting'.format(reg=region_number))
            return -1

    def reprocessed_evt2_for_ccd(self, ccd_id):
        return io.get_path("{evt2_file}[ccd_id={ccd_id}]".format(evt2_file=self.reprocessed_evt2_filename,
                                                                   ccd_id=ccd_id))

    def acis_ccd(self, ccd_id):
        return io.get_path("{analysis_dir}/acis_ccd{id}.fits".format(analysis_dir=self.analysis_directory,
                                                                     id=ccd_id))


class ClusterObj:
    """Cluster objects are intended to be pythonic representations of
    a galaxy cluster for processing through the Xray pypeline.

    """
    def __init__(self,
                 name="",
                 observation_ids=[],
                 data_directory="",
                 hydrogen_column_density=0,
                 redshift=0,
                 abundance=0,
                 last_step_completed=0,
                 ):
        """
        Initialization method.

        Parameters
        ----------

        name : str
            The name of the cluster. This variable is used to name directories and is the prefix for many filenames.
            There should be no spaces in the name.
            E.g. Abell 115 -> A115 or Abell_115
        observation_ids : str[]
            A list of chandra observation ids for download. These should all be for the same galaxy cluster.
            E.g. 3233, 15175, 15144
        data_directory : str
            The main data repository directory. This is the directory where the script will create a cluster.name
            subdirectory.
            E.g. /home/user/username/cluster_data/
        hydrogen_column_density : float
            The hydrogen column density for the cluster.
            E.g. 0.2
        redshift : float
            The redshift of the cluster
            E.g. 0.197
        abundance : float
            The solar abundance of the galaxy cluster.
            E.g. 0.2
        last_step_completed : int
            This is a state variable indicating the last step of the pypeline that was completed. It indicates where the
            pypeline will begin if run without any arguments.
        """
        self.name = name
        self.data_directory = data_directory
        self.hydrogen_column_density = hydrogen_column_density
        self.redshift = redshift
        self.abundance = abundance
        self._last_step_completed = last_step_completed
        self.observation_ids = observation_ids
        self.observations = [Observation(obsid=x, cluster=self) for x in self.observation_ids]

        #self.write_cluster_data()


    def write_cluster_data(self):
        """

        :return:
        """
        if ("" == self.name) or ("" == self.data_directory):
            assert "Trying to write before any work done"

        cluster_config = configparser.ConfigParser()
        cluster_dict = dict(self)
        cluster_dict['observation_ids'] = ','.join(self.observation_ids)
        cluster_config['cluster'] = cluster_dict

        try:
            with open(self.configuration_filename, 'w') as configfile:
                cluster_config.write(configfile)
        except FileExistsError:
            print("File exists and I can't overwrite! File: {}".format(self.configuration_filename))
            sys.exit(1)
        except FileNotFoundError:
            print("Cannot write cluster config to {}!".format(self.configuration_filename))
            print("Try updating your configuration file to reflect its current path.")
            sys.exit(1)
        print("Cluster data written to {}".format(self.configuration_filename))

        return

    def initialize_cluster(self):
        print("Initializing cluster object")
        self.get_cluster_info_from_user()
        io.make_directory(self.directory)
        self.write_cluster_data()
        io.make_initial_directories(self)
        print("Initialization complete.")  # Next step is to run the following command: ")
        print("Please continue running ClusterPyXT on {name}".format(name=self.name))

    def get_cluster_info_from_user(self):
        self.name = io.get_user_input("Enter the cluster name: ", "cluster name")
        self.data_directory = os.path.normpath(config.data_directory())

        self.observation_ids = get_observation_ids()
        self.observations = [Observation(obsid=x, cluster=self) for x in self.observation_ids]
        print()
        get_fitting_values = \
            io.check_yes_no("Enter values for fitting (nH, z, abundance) now? [y/n]")
        if get_fitting_values:
            self.hydrogen_column_density = io.get_user_input(
                "Enter the hydrogen column density for {} (on order of 10^22, e.g. 0.052 for 5.2e20): ".format(self.name),
                "hydrogen column density")
            self.redshift = io.get_user_input("Enter the redshift of {}: ".format(self.name), "redshift")

            self.abundance = io.get_user_input("Enter the abundance: ", "abundance")
        else:
            print("Before completing the ACB portion of the pypeline, you need "
                  "to edit the configuration file ({config}) "
                  "and update the values for hydrogen column density, redshift, "
                  "and abundance.".format(config=self.configuration_filename))
            self.hydrogen_column_density = "Update me! (on order of 10^22 e.g. 0.052 for 5.2e20)"
            self.redshift = "Update me! (e.g. 0.192)"
            self.abundance = "Update me! (e.g. 0.2)"
        self._last_step_completed = 0

        return

    def obs_directory(self, obsid):
        return io.get_path("{}/{}".format(self.directory, obsid))

    def obs_analysis_directory(self, obsid):
        return io.get_path("{}/analysis/".format(self.obs_directory(obsid)))

    def __repr__(self):
        return "<ClusterObj name:{}>".format(self.name)

    def __str__(self):
        ret_str = """Cluster: {}
Observations IDs: {}
Data Directory: {}
Hydrogen Column Density: {}
Redshift (z): {}
Abundance: {}
Last Step Completed: {}""".format(self.name,
                                  self.observation_ids,
                                  self.data_directory,
                                  self.hydrogen_column_density,
                                  self.redshift,
                                  self.abundance,
                                  self._last_step_completed
                                  )
        return ret_str

    def __iter__(self):
        yield 'name', self.name
        yield 'observation_ids', str(self.observation_ids)
        yield 'data_directory', str(self.data_directory)
        yield 'hydrogen_column_density', str(self.hydrogen_column_density)
        yield 'redshift', str(self.redshift)
        yield 'abundance', str(self.abundance)
        yield 'last_step_completed', str(self._last_step_completed)

        return

    @property
    def directory(self):
        """The directory containing the cluster data"""
        return io.get_path("{data_dir}/{cluster_name}/".format(
            data_dir=self.data_directory,
            cluster_name=self.name
        ))

    @property
    def configuration_filename(self):
        """The configuration file name"""
        return io.get_path("{}/{}/{}_pypeline_config.ini".format(self.data_directory, self.name, self.name))

    @property
    def observation_list(self):
        if isinstance(self.observation_ids, list):
            return self.observation_ids
        else:
            return [self.observation_ids]

    @property
    def last_step_completed(self):
        return self._last_step_completed

    @property
    def merged_directory(self):
        return io.get_path('{}/merged_obs_evt2/'.format(self.directory))

    @property
    def combined_directory(self):
        return io.get_path("{}/combined/".format(self.directory))

    @property
    def wvt_directory(self):
        return io.get_path("{}/wvt/".format(self.directory))

    @property
    def counts_image(self):
        return io.get_path("{combined_dir}/{cluster_name}_comb_ctsimage.fits".format(
            combined_dir=self.combined_directory,
            cluster_name=self.name
        ))

    @property
    def combined_signal(self):
        return io.get_path("{combined_dir}/{cluster_name}_comb_signal.fits".format(
            combined_dir=self.combined_directory,
            cluster_name=self.name
        ))

    @property
    def back_rescale(self):
        return io.get_path("{combined_dir}/{cluster_name}_comb_backrescl.fits".format(
            combined_dir=self.combined_directory,
            cluster_name=self.name
        ))

    @property
    def combined_mask(self):
        return io.get_path("{combined_dir}/acisI_comb_mask.fits".format(
            combined_dir=self.combined_directory
        ))

    @property
    def master_crop_file(self):
        return io.get_path("{combined_dir}/master_crop-ciaowcs.reg".format(
            combined_dir=self.combined_directory))


    @property
    def temp_acisI_comb(self):
        return io.get_path("{combined_dir}/acisI_comb_img_int.fits".format(combined_dir=self.combined_directory))


    @property
    def temp_backI_comb(self):
        return io.get_path("{combined_dir}/backI_comb_img_int.fits".format(combined_dir=self.combined_directory))


    @property
    def wvt_label_map(self):
        return io.get_path("{combined_dir}/{cluster}_wvt_label.fits".format(combined_dir=self.combined_dir,
                                                                            cluster=self.cluster_name))

    @property
    def mach_map_filename(self):
        return io.get_path("{output_dir}/{name}_mach_map.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def angle_map_filename(self):
        return io.get_path("{output_dir}/{name}_angle_map.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def mach_histogram_filename(self):
        return io.get_path("{output_dir}/{name}_mach_histogram.png".format(
            output_dir=self.output_dir,
            name=self.name
        ))


    @property
    def voronoi_region_file(self):
        return io.get_path("{wvt_dir}/voronoi.reg".format(wvt_dir=self.wvt_directory))

    @property
    def acb_dir(self):
        return io.get_path("{dir}/acb/".format(dir=self.directory))

    @property
    def scale_map_file(self):
        return io.get_path("{acb_dir}/{cluster_name}_pype_scale_map.fits".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def scale_map(self):
        #return fits.open(self.scale_map_file)[0].data
        return io.get_pixel_values(self.scale_map_file)

    @property
    def scale_map_header(self):
        from astropy.io import fits
        return fits.open(self.scale_map_file)[0].header

    @property
    def scale_map_region_file(self):
        return io.get_path("{acb_dir}/{cluster_name}_scale_map_region_index.fits".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def broad_flux_data(self):
        return fits.open(self.broad_flux_filename)[0].data

    @property
    def broad_flux_filename(self):
        return io.get_path("{cluster_dir}/{name}_broad_flux.img".format(
            cluster_dir=self.directory,
            name=self.name
        ))

    @property
    def scale_map_region_index(self):
        #return fits.open(self.scale_map_region_file)[0].data
        return io.get_pixel_values(self.scale_map_region_file)

    @property
    def number_of_regions(self):
        return self.scale_map_mask[np.where(self.scale_map_mask == 1)].size


    @property
    def scale_map_mask(self):
        scale_map = self.scale_map
        scale_map[np.nonzero(scale_map)] = 1
        return scale_map


    @property
    def sn_map(self):
        return io.get_path("{acb_dir}/{cluster_name}_pype_SN_map.fits".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @combined_mask.setter
    def combined_mask(self, filename):
        self._combined_mask = filename
        self.write_cluster_data()

    @last_step_completed.setter
    def last_step_completed(self, step):
        self._last_step_completed = step
        self.write_cluster_data()

    @property
    def region_list(self):
        return io.get_path("{acb_dir}/{cluster_name}_bin3_regionlist.list".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def filtered_region_list(self):
        return io.get_path("{acb_dir}/{cluster_name}_filtered_regions.list".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def region_to_index(self):
        return io.get_path("{acb_dir}/{cluster_name}_region_to_index.fits".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def combined_mask_region_index_map(self):
        #return fits.open(self.region_to_index)[0].data
        return io.get_pixel_values(self.region_to_index)


    @property
    def sources_file(self):
        return io.get_path("{cluster_dir}/sources.reg".format(cluster_dir=self.directory))

    @property
    def exclude_file(self):
        return io.get_path("{cluster_dir}/exclude.reg".format(cluster_dir=self.directory))

    @property
    def target_sn(self):
        return 40

    def spec_lis(self, region_number):
        return io.get_path("{super_comp_dir}/spec_{region_number}.lis".format(
            super_comp_dir=self.super_comp_dir,
            region_number=region_number
        ))

    @property
    def pi_directory(self):
        return io.get_path("{acb_dir}/pi_files/".format(
            acb_dir=self.acb_dir
            ))

    @property
    def spec_log_directory(self):
        return io.get_path("{acb_dir}/log_files/".format(
            acb_dir=self.acb_dir
            ))

    @property
    def super_comp_dir(self):
        return io.get_path("{acb_dir}".format( #/super_computer/".format(
            acb_dir=self.acb_dir
            ))

    @property
    def super_comp_cluster_config(self):
        return io.get_path("{super_comp_dir}/{name}_pypeline_config.ini".format(
            super_comp_dir=self.super_comp_dir,
            name=self.name
        ))

    @property
    def command_lis(self):
        return io.get_path("{super_comp_dir}/commands_{name}.lis".format(
            super_comp_dir=self.super_comp_dir,
            name=self.name
        ))

    @property
    def spec_fits_file(self):
        return io.get_path("{super_comp_dir}/{cluster_name}_best_spectral_fits.csv".format(
            super_comp_dir=self.super_comp_dir,
            cluster_name=self.name
        ))

    @property
    def all_fits_file(self):
        return io.get_path("{super_comp_dir}/{cluster_name}_all_spectral_fits.csv".format(
            super_comp_dir=self.super_comp_dir,
            cluster_name=self.name
        ))

    @property
    def bad_fits_file(self):
        return io.get_path("{super_comp_dir}/{cluster_name}_worst_spectral_fits.csv".format(
            super_comp_dir=self.super_comp_dir,
            cluster_name=self.name
        ))

    def acisI_clean_obs(self, observation_id):
        return io.get_path("{super_comp_dir}/acisI_clean_{obsid}.fits".format(
            super_comp_dir=self.super_comp_dir,
            obsid=observation_id
        ))

    def backI_clean_obs(self, observation_id):
        return io.get_path("{super_comp_dir}/backI_clean_{obsid}.fits".format(
            super_comp_dir=self.super_comp_dir,
            obsid=observation_id
        ))

    def effbtime_file_obs(self, observation_id):
        return io.get_path("{super_comp_dir}/effbtime-{obsid}_circle.dat".format(
            super_comp_dir=self.super_comp_dir,
            obsid=observation_id
        ))

    def effdtime_file_obs(self, observation_id):
        return io.get_path("{super_comp_dir}/effdtime-{obsid}_circle.dat".format(
            super_comp_dir=self.super_comp_dir,
            obsid=observation_id
        ))

    def scalemap_regionlist_file_obs(self, observation_id):
        return io.get_path("{super_comp_dir}/{name}_scalemap_regionlist_phys_{obsid}.reg".format(
            super_comp_dir=self.super_comp_dir,
            name=self.name,
            obsid=observation_id
        ))

    def sherpa_save_region(self, region_num):
        return io.get_path("{sherpa}/{region_num}.p".format(
            sherpa=self.sherpa_save_dir,
            region_num=region_num
        ))

    @property
    def sherpa_save_dir(self):
        return io.get_path("{super_comp_dir}/sherpa/".format(
            super_comp_dir=self.super_comp_dir
        ))

    def initialize_best_fits_file(self):
        self.initialize_fits_file(self.spec_fits_file)

    def initialize_all_fits_file(self):
        self.initialize_fits_file(self.all_fits_file)

    def initialize_worst_fits_file(self):
        self.initialize_fits_file(self.bad_fits_file)

    def initialize_fits_file(self, filename):
        with open(filename, 'w') as f:
            fieldnames = ['region', 'T', 'T_err_+', 'T_err_-',
                          'Norm', 'Norm_err_+', 'Norm_err_-',
                          'reduced_x2', 'obs_id']
            writer = csv.writer(f)
            writer.writerow(fieldnames)

    def write_best_fits_to_file(self, region=0,
                           T=0.0, T_err_plus=0.0,
                           T_err_minus=0.0, norm=0.0,
                           norm_err_plus=0.0, norm_err_minus=0.0,
                           reduced_x2=0.0, observation_id=0):
        self.write_fits_to_file(self.spec_fits_file, region=region,
                                T=T, T_err_plus=T_err_plus, T_err_minus=T_err_minus,
                                norm=norm, norm_err_plus=norm_err_plus,
                                norm_err_minus=norm_err_minus, reduced_x2=reduced_x2,
                                observation_id=observation_id)

    def write_fits_to_file(self, filename, region=0,
                           T=0.0, T_err_plus=0.0,
                           T_err_minus=0.0, norm=0.0,
                           norm_err_plus=0.0, norm_err_minus=0.0,
                           reduced_x2=0.0, observation_id=0):
        with open(filename, 'a') as f:
            writer = csv.writer(f)

            writer.writerow([region, T, T_err_plus, T_err_minus,
                             norm, norm_err_plus, norm_err_minus,
                             reduced_x2, observation_id
                             ])

    def write_bad_fits_to_file(self, region=0,
                               T=0.0, T_err_plus=0.0,
                               T_err_minus=0.0, norm=0.0,
                               norm_err_plus=0.0, norm_err_minus=0.0,
                               reduced_x2=0.0, observation_id=0):
        self.write_fits_to_file(self.bad_fits_file, region=region,
                                T=T, T_err_plus=T_err_plus, T_err_minus=T_err_minus,
                                norm=norm, norm_err_plus=norm_err_plus,
                                norm_err_minus=norm_err_minus, reduced_x2=reduced_x2,
                                observation_id=observation_id)

    def write_all_fits_to_file(self, region_number, fit_results,
                               confidences, observations):
        for i, observation in enumerate(observations):
            print("Writing fits from observation {} to file.".format(observation.id))
            T = confidences[i].parvals[0]
            T_err_plus = confidences[i].parmaxes[0]
            T_err_minus = confidences[i].parmins[0]
            norm = confidences[i].parvals[1]
            norm_err_plus = confidences[i].parmaxes[1]
            norm_err_minus = confidences[i].parmins[1]
            reduced_x2 = fit_results[i].rstat
            observation_id = observation.id

            self.write_fits_to_file(self.all_fits_file, region=region_number,
                                    T=T, T_err_plus=T_err_plus, T_err_minus=T_err_minus,
                                    norm=norm, norm_err_plus=norm_err_plus,
                                    norm_err_minus=norm_err_minus, reduced_x2=reduced_x2,
                                    observation_id=observation_id)

    @property
    def temperature_fits(self):
        temps = {'region': [],
                 'temperature': [],
                 'temp_err_plus': [],
                 'temp_err_minus': []
                 }

        with open(self.spec_fits_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                temps['region'].append(int(row['region']))
                temps['temperature'].append(float(row['T']))
                temps['temp_err_plus'].append(float(row['T_err_+']))
                temps['temp_err_minus'].append(float(row['T_err_-']))
        return temps

    @property
    def average_temperature_fits(self):
        temps = {'region': [],
                 'temperature': [],
                 'temp_err_plus': [],
                 'temp_err_minus': [],
                 'reduced_x2': []
                 }

        with open('/Users/brian/local/data/A115/acb/A115_average_fits_file.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                temps['region'].append(int(row['region']))
                temps['temperature'].append(float(row['T']))
                temps['temp_err_plus'].append(float(row['T_err_+']))
                temps['temp_err_minus'].append(float(row['T_err_-']))
                temps['reduced_x2'].append(float(row['reduced_x2']))
        return temps

    @property
    def pressure_map_filename(self):
        return io.get_path('{output_dir}/{name}_pressure.fits'.format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def temperature_map_filename(self):
        return io.get_path('{output_dir}/{cluster_name}_temperature_map.fits'.format(
            output_dir=self.output_dir,
            cluster_name=self.name
        ))

    @property
    def temperature_map(self):
        return io.get_pixel_values(self.temperature_map_filename)

    @property
    def temperature_map_header(self):
        from astropy.io import fits
        return fits.open(self.temperature_map_filename)[0].header

    @property
    def temperature_error_map_filename(self):
        return io.get_path('{output_dir}/{cluster_name}_temperature_error_map.fits'.format(
            output_dir=self.output_dir,
            cluster_name=self.name
        ))

    @property
    def temperature_fractional_error_map_filename(self):
        return io.get_path('{output_dir}/{cluster_name}_temperature_fractional_error_map.fits'.format(
            output_dir=self.output_dir,
            cluster_name=self.name
        ))

    def coordinates_for_scale_map_region(self, region, scale_map_regions):
        return np.where(scale_map_regions == region)

    @property
    def xray_surface_brightness_filename(self):
        return io.get_path("{cluster_dir}/{name}_broad_flux.img".format(
            cluster_dir=self.directory,
            name=self.name
        ))

    @property
    def xray_surface_brightness_cropped_filename(self):
        return io.get_path("{output_dir}/{name}_xray_surface_brightness_cropped.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def cropped_nosrc_xray_surface_brightness(self):
        return io.get_pixel_values(self.xray_surface_brightness_nosrc_filename)

    @property
    def density_map_filename(self):
        return io.get_path("{output_dir}/{name}_density.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))


    @property
    def xray_surface_brightness_nosrc_filename(self):
        return io.get_path("{output_dir}/{name}_xray_surface_brightness_nosrc.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def xray_surface_brightness_nosrc_cropped_filename(self):
        return io.get_path("{output_dir}/{name}_xray_surface_brightness_nosrc_cropped.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def output_dir(self):
        return io.get_path("{cluster_dir}/main_output/".format(
            cluster_dir=self.directory
        ))

    @property
    def regions_to_fit(self):
        with open(self.filtered_region_list, 'r') as f:
            raw_data = f.read()
        string_list = raw_data.split('\n')
        region_list = [int(x) for x in string_list]
        return np.array(region_list)

    def scale_map_regions_to_fit(self, resolution):
        resolutions = [None, 5, 3, 1]
        step_size = resolutions[resolution]
        region_map = self.scale_map_region_index
        shape = region_map.shape
        nx = shape[0]
        ny = shape[1]
        regions_to_fit = []
        for x in range(nx):
            for y in range(ny):
                if x % step_size == y % step_size == 0 and not np.isnan(region_map[x, y]):
                    regions_to_fit.append(int(region_map[x, y]))
        return regions_to_fit

    def observation(self, id):
        for obs in self.observations:
            if obs.id == id:
                return obs
        raise KeyError("Observation not found for observation {obsid}".format(obsid=id))

    def unfinished_regions_to_fit(self, resolution):
        all_regions = self.scale_map_regions_to_fit(resolution)
        complete_good = []
        complete_bad = []
        with open(self.spec_fits_file, 'r') as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                complete_good.append(int(row[0]))
        with open(self.bad_fits_file, 'r') as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                complete_bad.append(int(row[0]))
        complete_good = np.array(complete_good)
        complete_bad = np.array(complete_bad)

        print("For {reg} total regions, {good} good fits found, {bad} bad fits found.".format(
            reg=len(all_regions),
            good=len(complete_good),
            bad=len(complete_bad)
        ))

        for i in complete_good:
            all_regions = np.delete(all_regions, np.where(all_regions == i))
        for i in complete_bad:
            all_regions = np.delete(all_regions, np.where(all_regions == i))
        print("Number of regions left to fit: {reg}.".format(reg=len(all_regions)))
        return all_regions

    def write_effective_times_to_fits(self):
        header = self.scale_map_header

        import astropy.io.fits as fits

        for obs in self.observations:
            fits_data_filename = "{}.fits".format(obs.effective_data_time_file)
            fits_background_filename = "{}.fits".format(obs.effective_background_time_file)
            fits.writeto(fits_data_filename, obs.effective_data_time)
            fits.writeto(fits_background_filename, obs.effective_background_time)


def get_observation_ids():
    user_input_good = False
    observations = None
    while not user_input_good:
        user_input = input("Enter the observation ids to use (with a space to separate): ")
        split_user_input = user_input.split(' ')
        if [user_input] != split_user_input:
            y_no_flag = input("Are {} the observations you want to use? (y/n): ".format(user_input))
            if y_no_flag.lower() in ['y', 'yes']:
                user_input_good = True
            observations = split_user_input
        else:
            y_no_flag = input("Is {} the observation you want to use? (y/n): ".format(user_input))
            if y_no_flag.lower() in ['y', 'yes']:
                user_input_good = True
            observations = [user_input]
    return observations


def read_cluster_data(filename):

    cluster_config = configparser.ConfigParser()

    cluster_config.read(filename)

    try:
        cluster_dict = dict(cluster_config['cluster'])
    except KeyError as keyerror:
        print("Problem loading the cluster configuration file. Is {} correct?".format(filename))
        sys.exit(1)

    cluster = ClusterObj()
    cluster.name = cluster_dict['name']
    cluster.observation_ids = cluster_dict['observation_ids'].split(',')
    cluster.observation_ids = [obsid for obsid in cluster.observation_ids if obsid is not ',']
    cluster.data_directory = cluster_dict['data_directory']
    cluster.hydrogen_column_density = cluster_dict['hydrogen_column_density']
    cluster.redshift = cluster_dict['redshift']
    cluster.abundance = cluster_dict['abundance']
    cluster.last_step_completed = cluster_dict['last_step_completed']

    cluster.observations = [Observation(obsid=x, cluster=cluster) for x in cluster.observation_ids]

    return cluster


def get_cluster_config(clstr_name):
    data_dir = config.data_directory()
    config_file = io.get_filename_matching('{0}{1}/{1}_pypeline_config.ini'.format(data_dir, clstr_name))

    if len(config_file) >= 1:
        return config_file[-1]
    else:
        return None


