import os
import sys
import configparser
import csv
import config
import pypeline_io as io
import numpy as np
import astropy.io.fits as fits # migrate pycrates commands to fits
import logging

try:
    from ciao_contrib import runtool as rt
    from sherpa.astro import ui as sherpa
except ImportError:
    print("Unable to load CIAO. Fitting calls will fail")
    rt = None
    sherpa = None

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

            # if len(self.acis_S_chips) > 0:
            #     print("Observation {} has ACIS-S data that is not yet supported. "
            #           "Feel free to implement and submit a pull request!".format(self.id))
            #     # if len(self.acis_I_chips) > 0:
                #     print("ACIS-I chips to be used:")
                #     for ccd in self.acis_I_chips:
                #         print(ccd)
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
        if len(acis_I_ccds) > len(acis_S_ccds):
            self.acis_type = 0 # ACIS-I
            self.ccd_filter = "0:3"
        else:
            self.acis_type = 1 # ACIS-S
            self.ccd_filter = "4:9"
        return acis_I_ccds, acis_S_ccds

    def effective_data_time_for_region(self, region_number):
        coordinates_for_region = self.coordinates_for_scale_map_region(region_number,
                                                                       self.cluster.scale_map_region_index)
        return self.effective_data_time[coordinates_for_region][0]

    def effective_background_time_for_region(self, region_number):
        coordinates_for_region = self.coordinates_for_scale_map_region(region_number,
                                                                       self.cluster.scale_map_region_index)
        return self.effective_background_time[coordinates_for_region][0]

    def extract_spec_for_region(self, region_number):
        data_time = self.effective_data_time_for_region(region_number)
        back_time = self.effective_background_time_for_region(region_number)

        self.write_temp_region(region_number)
        region_file = self.temp_region_filename(region_number)

        infile = "{clean}[sky=region({region_file})][bin pi]".format(
            clean=self.sc_clean,
            region_file=region_file
        )

        outfile = io.get_path(
            "{super_comp_dir}/{obsid}_{region_number}.pi".format(
                super_comp_dir=self.cluster.super_comp_dir,
                obsid=self.id,
                region_number=region_number
            ))

        rt.dmextract(infile=infile, outfile=outfile, clobber=True)

        infile = "{back}[sky=region({region_file})][bin pi]".format(
            back=self.sc_back,
            region_file=region_file
        )

        outfile = io.get_path(
            "{super_comp_dir}/{obsid}_back_{region_number}.pi".format(
                super_comp_dir=self.cluster.super_comp_dir,
                obsid=self.id,
                region_number=region_number
            ))

        rt.dmextract(infile=infile, outfile=outfile, clobber=True)

        data_pi = "{super_comp_dir}/{obsid}_{region_number}.pi".format(
            super_comp_dir=self.cluster.super_comp_dir,
            obsid=self.id,
            region_number=region_number
        )

        back_pi = "{super_comp_dir}/{obsid}_back_{region_number}.pi".format(
            super_comp_dir=self.cluster.super_comp_dir,
            obsid=self.id,
            region_number=region_number
        )

        warf = "'{super_comp_dir}/{name}_{obsid}.arf'".format(
            super_comp_dir=self.cluster.super_comp_dir,
            name=self.cluster.name,
            obsid=self.id
        )

        wrmf = "'{super_comp_dir}/{name}_{obsid}.rmf'".format(
            super_comp_dir=self.cluster.super_comp_dir,
            name=self.cluster.name,
            obsid=self.id
        )

        # Put this background file into the 'grouped' data file for the region

        # rt.dmhedit(infile=data_pi, filelist="", operation="add", key="BACKFILE", value=back_pi)

        rt.dmhedit(infile=data_pi, filelist="", operation="add", key="EXPOSURE", value=data_time)
        rt.dmhedit(infile=data_pi, filelist="", operation="add", key="RESPFILE", value=wrmf)
        rt.dmhedit(infile=data_pi, filelist="", operation="add", key="ANCRFILE", value=warf)
        rt.dmhedit(infile=data_pi, filelist="", operation="add", key="BACKFILE", value=back_pi)
        rt.dmhedit(infile=back_pi, filelist="", operation="add", key="EXPOSURE", value=back_time)

        io.append_to_file(self.cluster.spec_lis(region_number), "{}\n".format(data_pi))
        io.delete(self.temp_region_filename(region_number))

        return data_pi, back_pi

    @property
    def directory(self):
        return io.get_path("{cluster_dir}/{obsid}/".format(cluster_dir=self.cluster.directory,
                                                           obsid=self.id))

    @property
    def ccd_filtered_reprocessed_evt2_filename(self):
        return "{evt2}[ccd_id={ccd_filter}]".format(
            evt2=self.reprocessed_evt2_filename,
            ccd_filter=self.ccd_filter
        )

    @property  # refactor so exclude points to the content of the file, not the file
    def exclude(self):
        return self.cluster.exclude_file

    @property
    def exclude_file(self):
        return self.cluster.exclude_file

    @property
    def cropped_clean_infile_string(self):
        return "{}[sky=region({})]".format(self.clean,
                                           self.cluster.master_crop_file)

    @property
    def cropped_background_infile_string(self):
        return "{}[sky=region({})]".format(self.back,
                                           self.cluster.master_crop_file)

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
    def mask_file(self):
        return io.get_filename_matching("{}/*msk1.fits".format(self.reprocessing_directory))[0]

    @property
    def analysis_directory(self):
        return io.get_path("{obs_dir}/analysis/".format(obs_dir=self.directory))

    @property
    def reprocessing_directory(self):
        # return io.get_path("{analysis_dir}/repro/".format(analysis_dir=self.analysis_directory))
        return io.get_path("{obs_dir}/repro/".format(obs_dir=self.directory))

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
    def broad_flux_filename(self):
        return io.get_path("{cluster_dir}/{cluster_name}_{obsid}_broad_flux.img".format(
            cluster_dir=self.cluster.directory,
            cluster_name=self.cluster.name,
            obsid=self.id
        ))

    @property
    def broad_flux(self):
        return io.get_pixel_values(self.broad_flux_filename)

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
            self.reprocessing_directory,
            self.id
        ))[0]

    @property
    def acisI_comb_img(self):
        return io.get_path("{combined_dir}/acisI_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id))

    @property
    def acisI_combined_image_filename(self):
        return self.acisI_comb_img

    @property
    def acisI_combined_image(self):
        if not hasattr(self, "_acisI_comb_image"):
            self._acisI_combined_image = io.get_pixel_values(self.acisI_comb_img)

        return self._acisI_combined_image

    @property
    def acisI_combined_image_header(self):
        if not hasattr(self, "_acisI_combined_image_header"):
            self._acisI_combined_image_header = fits.open(self.acisI_comb_img)[0].header

        return self._acisI_combined_image_header

    @property
    def backI_comb_img(self):
        return io.get_path("{combined_dir}/backI_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id)
        )

    @property
    def backI_comb_temp_img(self):
        return io.get_path("{combined_dir}/backI_comb_temp_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id)
        )

    @property
    def backI_combined_image(self):
        if not hasattr(self, "_backI_combined_image"):
            self._backI_combined_image = io.get_pixel_values(self.backI_comb_img)

        return self._backI_combined_image

    @property
    def backI_combined_image_header(self):
        if not hasattr(self, "_backI_combined_image_header"):
            self._backI_combined_image_header = fits.open(self.backI_comb_img)[0].header

        return self._backI_combined_image_header

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

    def reproject_combined_mask(self, file_to_match):
        from ciao import reproject
        print("Reprojecting combined mask for {}".format(self.id))
        temp_file = io.temp_filename(self.acisI_combined_mask_file)
        reproject(infile=self.acisI_combined_mask_file,
                  matchfile=file_to_match,
                  outfile=temp_file,
                  overwrite=True)
        io.move(temp_file, self.acisI_combined_mask_file)
        self._acisI_combined_mask = io.get_pixel_values(self.acisI_combined_mask_file)

    @property
    def acisI_combined_mask_file(self):
        return io.get_path("{combined_dir}/acisI_comb_mask-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def acisI_combined_mask(self):
        if not hasattr(self, "_acisI_combined_mask"):
            self._acisI_combined_mask = io.get_pixel_values(self.acisI_combined_mask_file)
        return self._acisI_combined_mask

    def reproject_nosrc_combined_mask(self, file_to_match):
        from ciao import reproject
        print("Reprojecting source removed combined mask for {}".format(self.id))
        temp_file = io.temp_filename(self.acisI_nosrc_combined_mask_file)
        reproject(infile=self.acisI_nosrc_combined_mask_file,
                  matchfile=file_to_match,
                  outfile=temp_file,
                  overwrite=True)
        io.move(temp_file, self.acisI_nosrc_combined_mask_file)
        self._acisI_nosrc_combined_mask = io.get_pixel_values(self.acisI_nosrc_combined_mask_file)

    @property
    def acisI_nosrc_combined_mask_file(self):
        return io.get_path("{combined_dir}/acisI_comb_mask_nosrc-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def acisI_nosrc_combined_mask(self):
        if not hasattr(self, '_acisI_nosrc_combined_mask'):
            self._acisI_nosrc_combined_mask = io.get_pixel_values(self.acisI_nosrc_combined_mask_file)
        return self._acisI_nosrc_combined_mask

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
    def temp_acis_comb_filename(self):
        return io.get_path("{combined_dir}/acis_comb_temp-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def temp_back_comb_filename(self):
        return io.get_path("{combined_dir}/back_comb_temp-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def temporary_acis_combined_energy_filtered_infile_string(self):
        return "{}[bin sky=4][energy=700:8000]".format(self.temp_acis_comb_filename)

    @property
    def temporary_back_combined_energy_filtered_infile_string(self):
        return "{}[bin sky=4][energy=700:8000]".format(self.temp_back_comb_filename)

    @property
    def acisI_high_energy_combined_image_file(self):
        return io.get_path("{combined_dir}/acisI_hien_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def acisI_high_energy_combined_image(self):
        if not hasattr(self, "_acisI_high_energy_combined_image"):
            self._acisI_high_energy_combined_image = \
                io.get_pixel_values(self.acisI_high_energy_combined_image_file)

        return self._acisI_high_energy_combined_image

    @property
    def acisI_high_energy_combined_image_header(self):
        if not hasattr(self, "_acisI_high_energy_combined_image_header"):
            self._acisI_high_energy_combined_image_header = \
                fits.open(self.acisI_high_energy_combined_image_file)[0].header
        return self._acisI_high_energy_combined_image_header


    @property
    def acisI_high_energy_temp_image(self):
        return io.get_path("{combined_dir}/acisI_high_energy_int_img_{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def backI_high_energy_combined_image_file(self):
        return io.get_path("{combined_dir}/backI_hien_comb_img-{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
        ))

    @property
    def backI_high_energy_combined_image(self):
        if not hasattr(self, "_backI_high_energy_combined_image"):
            self._backI_high_energy_combined_image = \
                io.get_pixel_values(self.backI_high_energy_combined_image_file)
        return self._backI_high_energy_combined_image

    @property
    def backI_high_energy_temp_image(self):
        return io.get_path("{combined_dir}/backI_hien_int_img_{obsid}.fits".format(
            combined_dir=self.combined_directory,
            obsid=self.id
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
        if not hasattr(self, "_scale_map_region_list"):
            region_list = []
            with open(self.scale_map_region_list_filename, 'r') as f:
                reader = csv.reader(f, delimiter='#')
                region_list = list(reader)
            self._scale_map_region_list = region_list

        return self._scale_map_region_list

    @scale_map_region_list.setter
    def scale_map_region_list(self, circle_list):
        with open(self.scale_map_region_list_filename, 'w') as f:
            writer = csv.writer(f, delimiter='#')
            writer.writerows(circle_list)

    @property
    def point_spread_function_map_filename(self):
        return io.get_path("{cluster_dir}/{name}_{obsid}_psfmap.fits".format(
            cluster_dir=self.cluster.directory,
            name=self.cluster.name,
            obsid=self.id
        ))

    @property
    def broad_threshold_image_filename(self):
        return io.get_path("{cluster_dir}/{name}_{obsid}_broad_thresh.img".format(
            cluster_dir=self.cluster.directory,
            name=self.cluster.name,
            obsid=self.id
        ))

    @property
    def source_map_filename(self):
        return io.get_path("{obs_dir}/{obsid}_source_events.fits".format(
            obs_dir=self.directory,
            obsid=self.id
        ))

    @property
    def source_image_filename(self):
        return io.get_path("{obs_dir}/{obsid}_image.fits".format(
            obs_dir=self.directory,
            obsid=self.id
        ))

    @property
    def normalized_background_without_sources_filename(self):
        return io.get_path("{obs_dir}/{obsid}_nbgd.fits".format(
            obs_dir=self.directory,
            obsid=self.id
        ))

    @property
    def source_cell_map_filename(self):
        return io.get_path("{obs_dir}/{obsid}_source_cell.fits".format(
            obs_dir=self.directory,
            obsid=self.id
        ))

    @property
    def source_region_filename(self):
        # return io.get_path("{obs_dir}/{obsid}_sources.reg".format(
        #     obs_dir=self.directory,
        #     obsid=self.id
        # ))
        return self.cluster.sources_file

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
        return io.get_filename_matching("{secondary_dir}/acis*evt1.fits".format(
            secondary_dir=self.secondary_directory
        ))
        # return io.get_filename_matching("{analysis_dir}/acis*evt1.fits".format(
        #     analysis_dir=self.analysis_directory
        # ))[0]

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
        filename = io.get_filename_matching("{repro_dir}/*_msk1.fits".format(
            repro_dir=self.reprocessing_directory
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
        if not hasattr(self, "_effective_data_time"):
            self._effective_data_time = np.load(self.effective_data_time_file)
        return self._effective_data_time

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
        if not hasattr(self, "_effective_background_time"):
            self._effective_background_time = np.load(self.effective_background_time_file)

        return self._effective_background_time

    @effective_background_time.setter
    def effective_background_time(self, effective_time):
        np.save(self.effective_background_time_file, effective_time)

    @property
    def ccd_merge_list(self):
        return io.get_path('{obs_analysis_dir}/acisI.lis'.format(
            obs_analysis_dir=self.analysis_directory)
        )

    @property
    def acisI_region_0_filename(self):
        return io.get_path('{obs_analysis_dir}/acisI_region_0.reg'.format(obs_analysis_dir=self.analysis_directory))

    @property
    def acisI_region_0_size(self):
        if not hasattr(self, "_acisI_region_0_size"):
            try:
                acisI_region_0 = io.read_contents_of_file(self.acisI_region_0_filename)
                radius = acisI_region_0.split(',')[-1][:-2]
                self._acisI_region_0_size = radius
            except FileNotFoundError:
                pass

        return self._acisI_region_0_size


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

    def temp_region_filename(self, region_number):
        if region_number == -1:
            return
        temp_region_file_name = io.get_path("{acb_dir}/temp_{region_number}_{obsid}.reg".format(
            acb_dir=self.cluster.acb_dir,
            region_number=region_number,
            obsid=self.id
        ))
        return temp_region_file_name

    def write_temp_region(self, region_number):
        region = self.get_region_from_region_number(region_number)
        io.write_contents_to_file(region, self.temp_region_filename(region_number), False)

    def reprocessed_evt2_for_ccd(self, ccd_id):
        return io.get_path("{evt2_file}[ccd_id={ccd_id}]".format(evt2_file=self.reprocessed_evt2_filename,
                                                                   ccd_id=ccd_id))

    def acis_ccd(self, ccd_id):
        return io.get_path("{analysis_dir}/acis_ccd{id}.fits".format(analysis_dir=self.analysis_directory,
                                                                     id=ccd_id))

# ********************************************************************
#
#   ClusterObj
#
# ********************************************************************


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
                 signal_to_noise=50
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
        self.signal_to_noise = float(signal_to_noise)

    def write_cluster_data(self):
        """

        :return:
        """
        if ("" == self.name) or ("" == self.data_directory):
            print("Trying to write cluster data file before any work complete")
            sys.exit(1)

        cluster_config = configparser.ConfigParser()
        cluster_dict = dict(self)
        cluster_dict['observation_ids'] = ','.join(self.observation_ids)
        cluster_config['cluster'] = cluster_dict

        try:
            with open(self.configuration_filename, 'w') as configfile:
                cluster_config.write(configfile)
                print("Cluster data written to {}".format(self.configuration_filename))
        except FileExistsError:
            print("File exists and I can't overwrite! File: {}".format(self.configuration_filename))
            sys.exit(1)
        except FileNotFoundError:
            if self.data_directory != config.sys_config.data_directory:
                self.data_directory = config.sys_config.data_directory
                self.write_cluster_data()
            else:
                print("Cannot write cluster config to {}!".format(self.configuration_filename))
                print("Try updating your configuration file to reflect its current path.")
                sys.exit(1)
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
        self.data_directory = os.path.normpath(config.sys_config.data_directory)

        self.observation_ids = get_observation_ids()
        self.observations = [Observation(obsid=x, cluster=self) for x in self.observation_ids]
        print()
        get_fitting_values = \
            io.check_yes_no("Enter values for fitting (nH, z, abundance, S/N) now? [y/n]")
        if get_fitting_values:
            self.hydrogen_column_density = io.get_user_input(
                "Enter the hydrogen column density for {} (on order of 10^22, e.g. 0.052 for 5.2e20): ".format(self.name),
                "hydrogen column density")
            self.redshift = io.get_user_input("Enter the redshift of {}: ".format(self.name), "redshift")

            self.abundance = io.get_user_input("Enter the abundance: ", "abundance")
            self.signal_to_noise = io.get_user_input("Enter the desired signal to noise ratio: ",
                                                     "the signal to noise ratio")
        else:
            print("Before completing the ACB portion of the pypeline, you need "
                  "to edit the configuration file ({config}) "
                  "and update the values for hydrogen column density, redshift, "
                  "and abundance.".format(config=self.configuration_filename))
            self.hydrogen_column_density = "Update me! (on order of 10^22 e.g. 0.052 for 5.2e20)"
            self.redshift = "Update me! (e.g. 0.192)"
            self.abundance = "Update me! (e.g. 0.2)"
            self.signal_to_noise = "Update me! (e.g. 50)"
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
Signal to Noise Ratio: {}
Last Step Completed: {}""".format(self.name,
                                  self.observation_ids,
                                  self.data_directory,
                                  self.hydrogen_column_density,
                                  self.redshift,
                                  self.abundance,
                                  self.signal_to_noise,
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
        yield 'signal_to_noise', str(self.signal_to_noise)
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
    def counts_image_filename(self):
        return io.get_path("{combined_dir}/{cluster_name}_comb_ctsimage.fits".format(
            combined_dir=self.combined_directory,
            cluster_name=self.name
        ))

    @property
    def counts_image(self):
        if not hasattr(self, '_counts_image'):
            cts_image = np.zeros(self.combined_mask_data.shape)
            try:
                for obs in self.observations:
                    cts_image += obs.acisI_combined_image
                self._counts_image = cts_image
            except FileNotFoundError:
                print("No combined image file yet.")
                return None
        return self._counts_image
        

    @property
    def combined_signal(self):
        return io.get_path("{combined_dir}/{cluster_name}_comb_signal.fits".format(
            combined_dir=self.combined_directory,
            cluster_name=self.name
        ))

    @property
    def back_rescale_filename(self):
        return io.get_path("{combined_dir}/{cluster_name}_comb_backrescl.fits".format(
            combined_dir=self.combined_directory,
            cluster_name=self.name
        ))

    @property
    def back_rescale(self):
        if not hasattr(self, "_back_rescale"):
            back_rescale = np.zeros(self.combined_mask_data.shape)
            try:
                for obs in self.observations:
                    t_obs = obs.acisI_combined_image_header['EXPOSURE']
                    t_back = obs.backI_combined_image_header['EXPOSURE']
                    back_rescale += (t_obs / t_back) * obs.backI_combined_image
                self._back_rescale = back_rescale
            except FileNotFoundError:
                print("No background combined images created yet.")
                return None
        return self._back_rescale

    @property
    def combined_mask(self):
        return io.get_path("{combined_dir}/acisI_comb_mask.fits".format(
            combined_dir=self.combined_directory
        ))

    @property
    def combined_mask_header(self):
        if not hasattr(self, "_combined_mask_header"):
            self._combined_mask_header = fits.open(self.combined_mask)[0].header

        return self._combined_mask_header

    @property
    def combined_mask_data(self):
        if not hasattr(self, "_combined_mask_data"):
            self._combined_mask_data = io.get_pixel_values(self.combined_mask)

        return self._combined_mask_data

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
        return io.get_path("{combined_dir}/{cluster}_wvt_label.fits".format(combined_dir=self.combined_directory,
                                                                            cluster=self.name))

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
    def scale_map_csv_filename(self):
        return io.get_path("{acb_dir}/{cluster_name}_scale_map_values.csv".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def scale_map_file(self):
        return io.get_path("{acb_dir}/{cluster_name}_pype_scale_map.fits".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def scale_map(self):
        if not hasattr(self, "_scale_map"):
            self._scale_map = io.get_pixel_values(self.scale_map_file)
        return self._scale_map

    @property
    def scale_map_indices(self):
        if not hasattr(self, "_scale_map_indices"):
            mask = self.combined_mask_data
            self._scale_map_indices = np.vstack(np.where(mask==1)).T
        return self._scale_map_indices

    @property
    def scale_map_header(self):
        if not hasattr(self, "_scale_map_header"):
            self._scale_map_header = fits.open(self.scale_map_file)[0].header

        return self._scale_map_header

    @property
    def scale_map_region_file(self):
        return io.get_path("{acb_dir}/{cluster_name}_scale_map_region_index.fits".format(
            acb_dir=self.acb_dir,
            cluster_name=self.name
        ))

    @property
    def broad_flux_data(self):
        if not hasattr(self, "_broad_flux_data"):
            self._broad_flux_data = io.get_pixel_values(self.broad_flux_filename)

        return self._broad_flux_data

    @property
    def broad_flux_filename(self):
        return io.get_path("{cluster_dir}/{name}_broad_flux.img".format(
            cluster_dir=self.directory,
            name=self.name
        ))

    @property
    def scale_map_region_index(self):
        if not hasattr(self, "_scale_map_region_index"):
            self._scale_map_region_index = io.get_pixel_values(self.scale_map_region_file)

        return self._scale_map_region_index

    @property
    def number_of_regions(self):
        if not hasattr(self, "_number_of_regions"):
            self._number_of_regions = self.scale_map_mask[np.where(self.scale_map_mask == 1)].size
        return self._number_of_regions


    @property
    def scale_map_mask(self):
        if not hasattr(self, "_scale_map_mask"):
            scale_map = np.copy(self.scale_map)
            scale_map[np.nonzero(scale_map)] = 1
            self._scale_map_mask = scale_map

        return self._scale_map_mask


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
        if not hasattr(self, "_combined_mask_region_index_map"):
            self._combined_mask_region_index_map = io.get_pixel_values(self.region_to_index)
        #return fits.open(self.region_to_index)[0].data
        return self._combined_mask_region_index_map


    @property
    def sources_file(self):
        return io.get_path("{cluster_dir}/sources.reg".format(cluster_dir=self.directory))

    @property
    def exclude_file(self):
        return io.get_path("{cluster_dir}/exclude.reg".format(cluster_dir=self.directory))

    @property
    def target_sn(self):
        return float(self.signal_to_noise)

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

    def initialize_scale_map_csv(self):
        with open(self.scale_map_csv_filename, 'w') as f:
            fieldnames = ['x', 'y', 'radius', 'signal_to_noise']
            writer = csv.writer(f)
            writer.writerow(fieldnames)

    def write_scale_map_radius(self, x=0, y=0, radius=0, signal_to_noise=0):
        with open(self.scale_map_csv_filename, 'a') as f:
            writer = csv.writer(f)
            writer.writerow([x, y, radius, signal_to_noise])

    def write_best_fits_to_file(self, region=0,
                           T=0.0, T_err_plus=0.0,
                           T_err_minus=0.0, norm=0.0,
                           norm_err_plus=0.0, norm_err_minus=0.0,
                           reduced_x2=0.0, observation_ids=""):
        self.write_fits_to_file(self.spec_fits_file, region=region,
                                T=T, T_err_plus=T_err_plus, T_err_minus=T_err_minus,
                                norm=norm, norm_err_plus=norm_err_plus,
                                norm_err_minus=norm_err_minus, reduced_x2=reduced_x2,
                                observation_ids=observation_ids)

    def write_fits_to_file(self, filename, region=0,
                           T=0.0, T_err_plus=0.0,
                           T_err_minus=0.0, norm=0.0,
                           norm_err_plus=0.0, norm_err_minus=0.0,
                           reduced_x2=0.0, observation_ids=""):
        with open(filename, 'a') as f:
            writer = csv.writer(f)

            writer.writerow([region, T, T_err_plus, T_err_minus,
                             norm, norm_err_plus, norm_err_minus,
                             reduced_x2, observation_ids
                             ])

    def write_bad_fits_to_file(self, region=0,
                               T=0.0, T_err_plus=0.0,
                               T_err_minus=0.0, norm=0.0,
                               norm_err_plus=0.0, norm_err_minus=0.0,
                               reduced_x2=0.0, observation_ids=""):
        self.write_fits_to_file(self.bad_fits_file, region=region,
                                T=T, T_err_plus=T_err_plus, T_err_minus=T_err_minus,
                                norm=norm, norm_err_plus=norm_err_plus,
                                norm_err_minus=norm_err_minus, reduced_x2=reduced_x2,
                                observation_ids=observation_ids)

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
                                    observation_ids=observation_id)

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

    def get_fits_from_file_for(self, fit_type='Norm'):
        err_high = "{fit_type}_err_+".format(fit_type=fit_type)
        err_low = "{fit_type}_err_-".format(fit_type=fit_type)

        fits = {'region': [],
                fit_type: [],
                err_high: [],
                err_low: []
        }
        with open(self.spec_fits_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                fits['region'].append(int(row['region']))
                fits[fit_type].append(float(row[fit_type]))
                try:
                    fits[err_high].append(float(row[err_high]))
                except ValueError:
                    fits[err_high].append(np.nan)               
                try:
                    fits[err_low].append(float(row[err_low]))
                except ValueError:
                    fits[err_low].append(np.nan)
            return fits

    @property
    def norm_fits(self):
        return self.get_fits_from_file_for(fit_type='Norm')


    @property
    def scale_map_csv_values(self):
        scale_map_values = {'x': [],
                            'y': [],
                            'radius': [],
                            'signal_to_noise': []
                            }
        with open(self.scale_map_csv_filename, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                scale_map_values['x'].append(int(row['x']))
                scale_map_values['y'].append(int(row['y']))
                scale_map_values['radius'].append(float(row['radius']))
                scale_map_values['signal_to_noise'].append(float(row['signal_to_noise']))

        return scale_map_values


    def write_scale_map_csv_to_fits(self):
        scale_map_values = self.scale_map_csv_values

        scale_map = np.zeros(self.combined_mask_data.shape)
        s_to_n_map = np.zeros(self.combined_mask_data.shape)

        for i in range(len(scale_map_values['x'])):
            x = scale_map_values['x'][i]
            y = scale_map_values['y'][i]
            radius = scale_map_values['radius'][i]
            s_to_n = scale_map_values['signal_to_noise'][i]

            scale_map[x,y] = radius
            s_to_n_map[x,y] = s_to_n

        header = self.combined_mask_header

        io.write_numpy_array_to_fits(scale_map, self.scale_map_file, header)
        io.write_numpy_array_to_fits(s_to_n_map, self.sn_map, header)

    def fit_map_filename(self, fit_type):
        return io.get_path("{output_dir}/{name}_{fit_type}.fits".format(
            output_dir=self.output_dir,
            name=self.name,
            fit_type=fit_type
        ))
    
    def fit_error_map_filename(self, fit_type):
        return io.get_path("{output_dir}/{name}_{fit_type}_error_map.fits".format(
            output_dir=self.output_dir,
            name=self.name,
            fit_type=fit_type
        ))

    def fit_fractional_error_map_filename(self, fit_type):
        return io.get_path("{output_dir}/{name}_{fit_type}_fractional_error_map.fits".format(
            output_dir=self.output_dir,
            name=self.name,
            fit_type=fit_type
        ))

    def fit_err_map_high_filename(self, fit_type):
        return io.get_path("{output_dir}/{name}_{fit_type}_high_error.fits".format(
            output_dir=self.output_dir,
            name=self.name,
            fit_type=fit_type
        ))
    
    def fit_err_map_low_filename(self, fit_type):
        return io.get_path("{output_dir}/{name}_{fit_type}_low_error.fits".format(
            output_dir=self.output_dir,
            name=self.name,
            fit_type=fit_type
        ))

    
    @property
    def norm_map(self):
        pass

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
    def pressure_map_high_filename(self):
        return io.get_path("{output_dir}/{name}_pressure_high_error.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def pressure_map_low_filename(self):
        return io.get_path("{output_dir}/{name}_pressure_low_error.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def entropy_map_filename(self):
        return io.get_path('{output_dir}/{name}_entropy.fits'.format(
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
    def temperature_map_temp_filename(self):
        return io.get_path('{output_dir}/{cluster_name}_temp_temperature_map.fits'.format(
            output_dir=self.output_dir,
            cluster_name=self.name
        ))

    @property
    def temperature_map(self):
        if not hasattr(self, "_temperature_map"):
            self._temperature_map = io.get_pixel_values(self.temperature_map_filename)

        return self._temperature_map

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
    def temperature_max_error_filename(self):
        return io.get_path("{output_dir}/{name}_temperature_map_max_error.fits".format(
            output_dir=self.output_dir,
            cluster_name=self.name    
        ))

    @property
    def temperature_min_error_filename(self):
        return io.get_path("{output_dir}/{name}_temperature_map_min_error.fits".format(
            output_dir=self.output_dir,
            name=self.name    
        ))

    @property
    def temperature_error_map(self):
        if not hasattr(self, "_temperature_map"):
            self._temperature_error_map = io.get_pixel_values(self.temperature_error_map_filename)
        return self._temperature_error_map

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
    def smoothed_xray_sb_cropped_nosrc_filename(self):
        return io.get_path("{output_dir}/{name}_acb_xray_sb_nosrc_cropped.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def cropped_nosrc_xray_surface_brightness(self):
        if not hasattr(self, "_cropped_nosrc_xray_surface_brightness"):
            self._cropped_nosrc_xray_surface_brightness = \
                io.get_pixel_values(self.xray_surface_brightness_nosrc_filename)

        return self._cropped_nosrc_xray_surface_brightness

    @property
    def density_map(self):
        if not hasattr(self, "_density_map"):
            self._density_map = io.get_pixel_values(self.density_map_filename)

        return self._density_map

    @property
    def density_map_filename(self):
        return io.get_path("{output_dir}/{name}_density.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def density_map_temp_filename(self):
        return io.get_path("{output_dir}/{name}_temp_density.fits".format(
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
    def xray_surface_brightness_nosrc_header(self):
        return fits.open(self.xray_surface_brightness_nosrc_filename)[0].header

    def parallel_observation_lists(self, num_cpus=1):
        import multiprocessing as mp
        num_cpus = num_cpus
        num_observations = len(self.observations)
        observation_lists = np.array_split(self.observations, (num_observations // num_cpus)+1)
        return observation_lists

    @property
    def xray_surface_brightness_nosrc_cropped_filename(self):
        return io.get_path("{output_dir}/{name}_xray_surface_brightness_nosrc_cropped.fits".format(
            output_dir=self.output_dir,
            name=self.name
        ))

    @property
    def xray_surface_brightness_nosrc_cropped_header(self):
        return fits.open(self.xray_surface_brightness_nosrc_cropped_filename)[0].header

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
            if int(obs.id) == int(id):
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

    def write_effective_times_to_fits(self, pdf=False, png=False):
        for obs in self.observations:
            fits_data_filename = "{}.fits".format(obs.effective_data_time_file)
            fits_background_filename = "{}.fits".format(obs.effective_background_time_file)
            io.write_numpy_array_to_fits(obs.effective_data_time, fits_data_filename, self.scale_map_header)
            io.write_numpy_array_to_fits(obs.effective_background_time, fits_background_filename, self.scale_map_header)

            if pdf:
                pdf_data_filename = "{}.pdf".format(obs.effective_data_time_file)
                pdf_background_filename = "{}.pdf".format(obs.effective_background_time_file)
                io.write_numpy_array_to_image(obs.effective_data_time, pdf_data_filename)
                io.write_numpy_array_to_image(obs.effective_background_time, pdf_background_filename)

            if png:
                png_data_filename = "{}.png".format(obs.effective_data_time)
                png_background_filename = "{}.png".format(obs.effective_background_time)
                io.write_numpy_array_to_image(obs.effective_data_time, png_data_filename)
                io.write_numpy_array_to_image(obs.effective_background_time, png_background_filename)

    def fit_region_number(self, region_number, output_pdf=False):
        sherpa_log = logging.getLogger("Sherpa")
        sherpa_log.setLevel(logging.ERROR)

        print("Processing region number: {reg_num}".format(
            reg_num=region_number
        ))

        data_pi_files = []
        background_pi_files = []
        good_observations = []
        for observation in self.observations:
            print("{region}:\tWorking on observation id: {obsid}".format(
                obsid=observation.id,
                region=region_number
            ))

            effective_data_time_for_region = observation.effective_data_time_for_region(region_number)

            exposure_time = observation.exposure_time
            signal_to_noise = 10000 * (effective_data_time_for_region / exposure_time)
            if signal_to_noise >= self.signal_to_noise_threshold:
                print('{region}:\tMet S/N threshold for {obsid}'.format(
                    obsid=observation.id,
                    region=region_number
                ))
                # extract cleaned spec & background
                print("{region}:\tExtracting PI files for {obsid}".format(region=region_number,
                                                                          obsid=observation.id))
                data_pi, back_pi = observation.extract_spec_for_region(region_number)
                good_observations.append(observation.id)
                data_pi_files.append(data_pi)
                background_pi_files.append(back_pi)

        #cluster = observation.cluster
        number_of_observations = len(good_observations)

        print("{region}:\tLoading data pulse invariant files (PI files)".format(region=region_number))
        for i, data_pi in enumerate(data_pi_files):
            sherpa.load_pha(i, data_pi)
            # this should load the arf and rmf files automatically
            # they are set in the previous function call (ciao.extract_spec())

        print("{region}:\tLoading background PI files".format(region=region_number))
        for i, background_pi in enumerate(background_pi_files):
            sherpa.load_bkg(i, background_pi)


        #print("Subtracting the background")
        # subtract the background from the observation
        for i in range(number_of_observations):
            sherpa.subtract(i)

        sherpa.set_analysis('energy')

        # usually the data is rather noisey below 0.7 keV and above 8.0 keV.
        # For low signal to noise regions, the high energy cutoff may need to
        # be lowered down (to 5 keV for example, look at the data).
        low_energy_cutoff = 0.7  # [units: keV]
        high_energy_cutoff = 8.0  # [units: keV]
        sherpa.ignore(":{loE}, {hiE}:".format(
            loE=low_energy_cutoff,
            hiE=high_energy_cutoff
        ))

        #setting the model for each observation
        continue_fit = False
        for i in range(number_of_observations):
            continue_fit = True
            sherpa.set_source(i, sherpa.xsphabs.phabs*sherpa.xsapec.apec)

        if continue_fit:
            print("{region}:\tCreating the model and defining initial fit parameters".format(region=region_number))
            phabs.nH = float(self.hydrogen_column_density) / 1e22 
            apec.kT = 8.0
            apec.Abundanc = self.abundance
            apec.redshift = self.redshift
            apec.norm = 1.0

            print("{region}:\tFreezing and thawing parameters".format(region=region_number))
            sherpa.freeze(phabs.nH, apec.Abundanc, apec.redshift)
            sherpa.thaw(apec.kT, apec.norm)

            print("Fitting region: {}".format(region_number))
            sherpa.fit()
            sherpa.conf()

            fit_results = sherpa.get_fit_results()
            confidences = sherpa.get_conf_results()

            # parameter names given in kT, Norm order
            T = confidences.parvals[0]
            T_err_plus = confidences.parmaxes[0]
            T_err_minus = confidences.parmins[0]
            norm = confidences.parvals[1]
            norm_err_plus = confidences.parmaxes[1]
            norm_err_minus = confidences.parmins[1]
            reduced_x2 = fit_results.rstat
            observations = ','.join(good_observations)

            print("{region}:\tObservations used:\t{obs}\n"
                  "Reduced X2:\t{rx2}\n"
                  "Temperature:\t{T} keV\n".format(region=region_number,
                                                   obs=good_observations,
                                                   rx2=reduced_x2,
                                                   T=T))

            if T_err_plus == None:
                self.write_bad_fits_to_file(region=int(region_number),
                                            T=T,
                                            T_err_plus=T_err_plus,
                                            T_err_minus=T_err_minus,
                                            norm=norm,
                                            norm_err_plus=norm_err_plus,
                                            norm_err_minus=norm_err_minus,
                                            reduced_x2=reduced_x2,
                                            observation_ids=observations
                                            )
            else:
                self.write_best_fits_to_file(region=int(region_number),
                                             T=T,
                                             T_err_plus=T_err_plus,
                                             T_err_minus=T_err_minus,
                                             norm=norm,
                                             norm_err_plus=norm_err_plus,
                                             norm_err_minus=norm_err_minus,
                                             reduced_x2=reduced_x2,
                                             observation_ids=observations)

    # cluster.write_all_fits_to_file(int(region_number),
            #                                fit_results,
            #                                confidences,
            #                                cluster.observations)

            for data_pi in data_pi_files:
                io.delete(data_pi)
            for background_pi in background_pi_files:
                io.delete(background_pi)
            io.delete(self.spec_lis(region_number))

            print("{region}:\tFinished".format(region=region_number))


            ### TODO update the below code to exclude pychips
            if output_pdf:
                print("Update to exclude pychips")
            #     sherpa.plot_fit()
            #     pychips.set_preference('export.orientation', 'landscape')
            #     pychips.set_preference('export.clobber', 'True')
            #     pychips.set_preference('export.fittopage', 'True')
            #     pychips.print_window("{acb_dir}/{region}.pdf".format(
            #         acb_dir=self.acb_dir,
            #         region=region_number)
            #     )

            return "{region}: {temperature} keV".format(region=region_number,
                                                        temperature=T)
        else:
            io.print_red("No obs to fit for region {}".format(region_number))
   
    @property
    def signal_to_noise_threshold(self):
        return 900


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

    try:
        cluster_config.read(filename)
    except TypeError:
        print("Problem loading configuration file {}.".format(filename))

    try:
        cluster_dict = dict(cluster_config['cluster'])
    except KeyError as keyerror:
        print("Problem loading the cluster configuration file. Is {} correct?".format(filename))
        sys.exit(1)

    cluster = ClusterObj()
    cluster.name = cluster_dict['name']
    cluster.observation_ids = cluster_dict['observation_ids'].split(',')
    cluster.observation_ids = [obsid for obsid in cluster.observation_ids if obsid != ',']
    cluster.data_directory = cluster_dict['data_directory']
    cluster.hydrogen_column_density = cluster_dict['hydrogen_column_density']
    cluster.redshift = cluster_dict['redshift']
    cluster.abundance = cluster_dict['abundance']
    cluster.last_step_completed = cluster_dict['last_step_completed']
    try:
        cluster.signal_to_noise = cluster_dict['signal_to_noise']
    except KeyError:
        cluster.signal_to_noise = 50

    cluster.observations = [Observation(obsid=x, cluster=cluster) for x in cluster.observation_ids]

    # print("Loaded {}".format(cluster))

    return cluster


def load_cluster(cluster_name: str):
    if io.file_exists(cluster_name):
        return read_cluster_data(cluster_name)

    config_file = get_cluster_config(cluster_name)
    return read_cluster_data(config_file)


def get_cluster_config(clstr_name):
    filename = "{data_dir}/{name}/{name}_pypeline_config.ini".format(
            data_dir=config.sys_config.data_directory,
            name=clstr_name
        )
    config_file = io.get_filename_matching(filename)
    
    if len(config_file) >= 1:
        return config_file[-1]
    else:
        return None
