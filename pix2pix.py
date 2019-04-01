import ciao
import pypeline_io as io
import numpy as np
import cluster
import sys
from sherpa.astro import ui as sherpa


def pix_to_pix(cluster: cluster.ClusterObj, region_number, process='main'):
    print("Processing region number: {reg_num}".format(
        reg_num=region_number
    ))

    scale_map_region_index = cluster.scale_map_region_index

    data_pi_files = []
    background_pi_files = []
    good_observations = []
    for observation in cluster.observations:
        print("{proc}:\tWorking on obsid: {obsid}, region: {region}".format(
            proc=process,
            obsid=observation.id,
            region=region_number
        ))

        effective_data_time_for_region = observation.effective_data_time_for_region(region_number)

        exposure_time = observation.exposure_time
        signal_to_noise = 10000 * (effective_data_time_for_region / exposure_time)
        print('Effective signal to noise for observation: {obs_id} in region {region}'.format(
            obs_id=observation.id,
            region=region_number
        ))
        if signal_to_noise >= 900:
            #extract cleaned spec & background
            data_pi, back_pi = ciao.extract_spec(observation,
                                                 region_number)
            good_observations.append(observation.id)
            data_pi_files.append(data_pi)
            background_pi_files.append(back_pi)

    #cluster = observation.cluster
    number_of_observations = len(good_observations)

    print("Loading data pulse invariant files (pi files)")
    for i, data_pi in enumerate(data_pi_files):
        sherpa.load_pha(i, data_pi)
        # this should load the arf and rmf files automatically
        # they are set in the previous function call (ciao.extract_spec())

    print("Loading background files")
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
    for i in range(number_of_observations):
        sherpa.set_source(i, sherpa.xsphabs.phabs*sherpa.xsapec.apec)

    print("{proc}:\tCreating the model and defining initial fit parameters".format(proc=process))
    phabs.nH = cluster.hydrogen_column_density
    apec.kT = 8.0
    apec.Abundanc = cluster.abundance
    apec.redshift = cluster.redshift
    apec.norm = 1.0

    print("{proc}:\tFreezing and thawing parameters".format(proc=process))
    sherpa.freeze(phabs.nH, apec.Abundanc, apec.redshift)
    sherpa.thaw(apec.kT, apec.norm)

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

    print("Observations used:\t{obs}\n"
          "Reduced X2:\t{rx2}\n"
          "Temperature:\t{T} keV\n".format(obs=good_observations,
                                           rx2=reduced_x2,
                                           T=T))

    if T_err_plus == None:
        cluster.write_bad_fits_to_file(region=int(region_number),
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
        cluster.write_best_fits_to_file(region=int(region_number),
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
    io.delete(cluster.spec_lis(region_number))

    print("{proc}:\tFinished region number: {region}".format(proc=process, region=region_number))


def index_of_best_fit(fits, confs):
    best_fit = 0
    best_fit_index = 0
    for i, fit in enumerate(fits):
        if confs[i].parmaxes[0] is not None:  # T_error_plus is non-zero (i.e. the fit was ok)
            if np.abs(1-fit.rstat) < np.abs(1-best_fit):
                best_fit_index = i
                best_fit = fit.rstat

    return best_fit_index


if __name__ == "__main__":
    if len(sys.argv) == 3:
        cluster_file = sys.argv[1]
        region = int(sys.argv[2])
        cluster_obj = cluster.read_cluster_data(cluster_file)

        pix_to_pix(cluster_obj, region)
    else:
        print("Error when calling pix2pix.")
        print("Should be called python pix2pix.py /path/to/cluster_config.ini region_number")
