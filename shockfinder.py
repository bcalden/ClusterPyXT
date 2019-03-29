'''
This version of the shockfinder first identifies the angle at each pixel which gives the largest Mach number.
The X-ray jump is then tested at this angle to verify if it satisfies the criterion that both X-ray and
temperature should increase in the same direction. The other version (newshock_weak.py) tests the criterion at
every angle and saves the maximum jump which satisfies the angle. This makes it so no jump is detected only
if no direction satisfies the criterion. In this code, the criterion must be satisfied for the direction
which gives the maximum Mach number.
'''

import numpy as np
import astropy.io.fits as pf
try:
    from matplotlib import pylab as pl
except ImportError:
    print("Shockfinder cannot run with CIAO running. Please reset your environment and try again.")
    raise

from math import *
import cluster


def find_shock_in(clstr: cluster.ClusterObj):
    # Read in data
    # xray is smoothed x-ray map
    # temp is smoothed temperature map
    # res is the resolution map used for smoothing
    # res is rescaled by 1.25 because that was done for smoothing

    # xray=pf.getdata('obs_sim5_sm.fits')
    # temp=pf.getdata('Tmap_vor5_sm.fits')
    # res=pf.getdata('acisI_scale.fits')
    # res=1.25*res

    xray = pf.getdata(clstr.combined_signal)  # combined_signal
    # temp=pf.getdata('A3667_shocks_8p/A3667_Tmap_c8_xspec_mg_hybCF.fits')
    temp = pf.getdata(clstr.temperature_map_filename)  # temperature map
    thead = pf.getheader(clstr.temperature_map_filename)  # temperature map
    res = pf.getdata(clstr.scale_map_file)
    # res=pf.getdata('A3667_shocks_8p/A3667_c_kmap.fits')
    # data_dir='../shockfinder_new/Obs/observational_shockfinder/simulated/projected/spin_cluster/'


    # xray=pf.getdata('simulated/projections_XRayEmissivity_1_wvtbin_asm_rb.fits')
    # temp=pf.getdata('simulated/projections_Temperature_1_wvtbin_asm_rb.fits')
    # thead=pf.getheader('simulated/projections_Temperature_1_wvtbin_asm_rb.fits')
    # res=pf.getdata('simulated/sim_scale_rb.fits')

    xray[temp == 0] = 0;

    res = 1.25 * res

    # The sizes of the arrays should be the same
    nx = temp.shape[0]
    ny = temp.shape[1]

    mach = np.zeros((nx, ny), dtype='d')
    angle = np.zeros((nx, ny), dtype='d')
    i_plus = np.zeros((nx, ny), dtype='i')
    i_minus = np.zeros((nx, ny), dtype='i')
    j_plus = np.zeros((nx, ny), dtype='i')
    j_minus = np.zeros((nx, ny), dtype='i')
    Tjump = np.zeros((nx, ny), dtype='d')
    Xjump = np.zeros((nx, ny), dtype='d')
    max_Tjump = np.ones((nx, ny), dtype='d')
    max_theta = np.zeros((nx, ny), dtype='d')
    Xyes = np.zeros((nx, ny))
    mask = np.ones((nx, ny))
    XT = np.zeros((nx, ny), dtype='d')
    i = np.zeros((nx, ny), dtype='i')
    j = np.zeros((nx, ny), dtype='i')
    for k in range(ny):
        i[:, k] = np.arange(nx)
    for k in range(nx):
        j[k:, ] = np.arange(ny)

    theta = 0.0
    while theta < np.pi:
        print(theta)
        # Calculate the plus and minus indices given the angle and resolution map
        i_plus = np.rint(i + res * sin(theta))
        i_minus = i + i - i_plus
        j_plus = np.rint(j + res * cos(theta))
        j_minus = j + j - j_plus
        i_plus = i_plus.astype(int)
        i_minus = i_minus.astype(int)
        j_plus = j_plus.astype(int)
        j_minus = j_minus.astype(int)

        # Make sure all of the indices are contained in the image
        # If any are not, set the index to zero
        i_plus[i_plus > nx - 1] = 0
        i_plus[i_plus < 0] = 0
        i_minus[i_minus > nx - 1] = 0
        i_minus[i_minus < 0] = 0
        j_plus[j_plus > ny - 1] = 0
        j_plus[j_plus < 0] = 0
        j_minus[j_minus > ny - 1] = 0
        j_minus[j_minus < 0] = 0

        # Calculate the temperature jump and X-ray jump if both temperatures are non-zero
        pix1 = temp[i_plus, j_plus] > 0.0
        pix2 = temp[i_minus, j_minus] > 0.0
        pix = pix1 * pix2
        print("temp = ", temp[i_minus[pix], j_minus[pix]])
        print("xray = ", xray[i_minus[pix], j_minus[pix]])
        Tjump[pix] = temp[i_plus[pix], j_plus[pix]] / temp[i_minus[pix], j_minus[pix]]
        Xjump[pix] = xray[i_plus[pix], j_plus[pix]] / xray[i_minus[pix], j_minus[pix]]

        # Assign maximum values to the arrays
        plus_pix = Tjump > max_Tjump
        minus_pix = 1.0 / Tjump > max_Tjump
        max_Tjump[plus_pix] = Tjump[plus_pix]
        max_theta[plus_pix] = theta
        max_Tjump[minus_pix] = 1.0 / Tjump[minus_pix]
        max_theta[minus_pix] = theta + np.pi

        # Test if the temperature and X-ray jumps are in the same direction
        # XT equals 1 when they are in the same direction and 0 when they are not
        # XT is only calculated for pixels which achieved a max Mach number for the current theta
        XT[plus_pix] = np.sign((Tjump[plus_pix] - 1.0) * (Xjump[plus_pix] - 1.0)) * 0.5 + 0.5
        XT[minus_pix] = np.sign((Tjump[minus_pix] - 1.0) * (Xjump[minus_pix] - 1.0)) * 0.5 + 0.5

        theta = theta + np.pi / 32.0

    mach[XT == 1.0] = 0.2 * (-7.0 + 8.0 * max_Tjump[XT == 1.0] + 4.0 * np.sqrt(
        4.0 - 7.0 * max_Tjump[XT == 1.0] + 4.0 * max_Tjump[XT == 1.0] ** 2))
    angle[XT == 1.0] = max_theta[XT == 1.0] + 180. / np.pi;  # in degrees
    angle[XT != 1.0] = -1.0
    mach = np.sqrt(mach)
    angle[np.isnan(mach)] = -1.0
    mach[np.isnan(mach)] = 0.0
    mach[temp == 0.0] = np.nan
    angle[temp == 0.0] = np.nan

    for ii in range(0, nx):
        for jj in range(0, ny):
            if mach[ii, jj] > 5: mach[ii, jj] = 0

    # pf.writeto('A3667_hyb_strong_mach_new_rb.fits',mach,header=thead,clobber=True)
    # pf.writeto('A3667_hyb_strong_angle_new_rb.fits',angle,header=thead,clobber=True)

    pf.writeto(clstr.mach_map_filename, mach, header=thead, overwrite=True)
    pf.writeto(clstr.angle_map_filename, angle, header=thead, overwrite=True)

    # The code below was lifted from Sam's shockfinder
    # It is used to make a histogram of the Mach number

    # # Make histogram
    sa2d, mach2d = pl.histogram(mach, bins=50, range=[0.1, 4.0])
    # sa2d = 1.0*sa2d*(1.0*fov/pixels)**2

    #fig = pl.figure(figsize=(6.0, 6.0), dpi=200)
    #
    pl.clf()
    pl.semilogy(mach2d[1:], sa2d / (mach2d[1] - mach2d[0]), 'k', ls='steps-')
    #
    pl.xlabel('Mach')
    # pl.ylabel(r'd(Surface Area [$Mpc^2$])/d(Mach)')
    pl.ylabel(r'pixels')
    pl.xlim(0.1, 4.0)
    pl.ylim(1e2, 1.5e6)
    # pl.legend(['2D Shocks'],loc='lower left')
    pl.savefig(clstr.mach_histogram_filename, dpi=200, bbox_inches='tight')
