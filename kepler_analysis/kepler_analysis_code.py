from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import math
import astropy.io.fits as pyfits
import pandas as pd
from astropy.wcs import WCS
from scipy import signal

from astroquery.mast import Mast
from astroquery.mast import Observations


### DEFINITIONS ###

def open_file(filename, add_back = 'no'):
    
    h = pyfits.open(filename)
    data = h[1].data
    hdr = h[1].header
    #print(h.info())
    h.close() # close the file
    
    time = np.array(data['TIME']).astype(float)
    flux = np.array(data['FLUX']).astype(float)
    flux_err = np.array(data['FLUX_ERR']).astype(float)
    
    if add_back != 'no':
        print('Adding background')
        back_flux = np.array(data['FLUX_BKG']).astype(float)
        mean_back_flux_frame = np.nanmean(np.array(back_flux), axis = 0)
        flux = flux+mean_back_flux_frame
    
    return time, flux, flux_err, hdr
    
def filter_parameters(time, flux, flux_err):
    
    '''remove nan from time, flux, and flux_err as well
    as flux = zero data points'''
    
    #finding the index of the middle of the image
    half_row_len = int(len(flux[0])/2)
    half_col_len = int(len(flux[0][0])/2)
    
    #filtering nans from time array
    time_isnan = []
    for i in time:
        time_isnan.append(~np.isnan(i))
    #applying to arrays
    time = time[time_isnan]
    flux = flux[time_isnan]
    flux_err = flux_err[time_isnan]
    
    flux_isnan = []
    for i in flux:
        value = i[half_row_len][half_col_len]
        if value == 0:
            value = np.nan
        flux_isnan.append(~np.isnan(value))
    #applying to arrays
    time = time[flux_isnan]
    flux = flux[flux_isnan]
    flux_err = flux_err[flux_isnan]
    
    flux_err_isnan = []
    for i in flux_err:
        flux_err_isnan.append(~np.isnan(i[half_row_len][half_col_len]))
    #applying to arrays
    time = time[flux_err_isnan]
    flux = flux[flux_err_isnan]
    flux_err = flux_err[flux_err_isnan]
    
    return time, flux, flux_err
    
def get_wcs_coords(hdr):
    
    wcs_input_dict = {
        'CTYPE1': hdr['1CTYP4'],
        'CUNIT1': hdr['1CUNI4'],
        'CDELT1': hdr['1CDLT4'],
        'CRPIX1': hdr['1CRPX4'],
        'CRVAL1': hdr['1CRVL4'],
        'NAXIS1': hdr['NAXIS1'],
        'CTYPE2': hdr['2CTYP4'],
        'CUNIT2': hdr['2CUNI4'],
        'CDELT2': hdr['2CDLT4'],
        'CRPIX2': hdr['2CRPX4'],
        'CRVAL2': hdr['2CRVL4'],
        'NAXIS2': hdr['NAXIS2']
        }
    wcs_dict = WCS(wcs_input_dict)
    
    return wcs_dict
    
def sum_img_flux(time, flux, flux_err, plot = 'no', savefig_filename = 'no'):
    
    time_filt, flux_filt, flux_err_filt = filter_parameters(time, flux, flux_err)
    
    flux_filt_sum = []
    for i in flux_filt:
        flux_sum_row = []
        for j in i:
            flux_sum_row.append(sum(j))
        flux_sum_row = np.nan_to_num(flux_sum_row, nan = 0.0)
        flux_filt_sum.append(sum(flux_sum_row))
    
    if plot == 'yes':
        plt.figure(figsize = (9, 5))
        #plt.title('Light Curve Pixel(row, col) = (%s, %s)' %(row_ind, col_ind))
        plt.plot(time_filt, flux_filt_sum, color = 'purple')
        plt.xlabel('Time', fontsize = 15)
        plt.ylabel('Sum of Image Flux', fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        if savefig_filename != 'no':
            plt.savefig(savefig_filename)
        plt.show()
    
    #frame_ave_flux_sum = np.array(flux_sum)/(len(flux))
    
    return time_filt, flux_filt_sum
    

def calc_flux_over_area(time, flux, flux_err, ind, sum_area):
    
    time_filt, flux_filt, flux_err_filt = filter_parameters(time[ind], flux[ind], flux_err[ind])
    
    #finding middle coord of first image
    half_row_len = int(len(flux[0])/2)
    half_col_len = int(len(flux[0][0])/2)

    summed_pupil_fluxes = []
    for i in range(len(flux_filt)):
        img_flux = flux_filt[i]
        
        row_fluxes = img_flux[int(half_row_len-sum_area/2):int(half_row_len+sum_area/2)]
        pupil = row_fluxes.T[int(half_col_len-sum_area/2):int(half_col_len+sum_area/2)]
        summed_pupil_fluxes.append(np.sum(pupil))
        #print('percent done:', i/len(flux))
    
    return summed_pupil_fluxes
    
def calc_sample_frequency(time, flux, flux_err, ind):
    
    time_filt, flux_filt, flux_err_filt = filter_parameters(time[ind], flux[ind], flux_err[ind])
    
    time_steps = []
    for i in range(len(time_filt)-1):
        time_steps.append(time_filt[i+1]-time_filt[i])
    
    frequency = 1/np.median(time_steps)
    
    return frequency

def make_periodogram(time, flux, flux_err, filename_all, start_ind, file_num, sum_area, savefig = 'no'):

    f_all = []
    Pxx_den_all = []
    fnames = []

    ind = start_ind
    for i in range(start_ind, start_ind+file_num):
        #calculating light curves over pupil area
        lc = calc_flux_over_area(time, flux, flux_err, ind, sum_area)
        #calculating observation frequency
        freq = calc_sample_frequency(time, flux, flux_err, ind)
        #calculating periodogram signal
        f, Pxx_den = signal.periodogram(lc, freq)
        f_all.append(f)
        fnames.append(filenames_all[ind])
        Pxx_den_all.append(Pxx_den)
        ind += int(len(time)/3)
        
        
    
    fig, axs = plt.subplots(1, file_num, figsize = (13,3))
    
    axs_all = []
    for i in range(file_num):
        a = axs[i].plot(f_all[i], Pxx_den_all[i], color = 'black', linewidth = 1)
        axs[i].set_title(fnames[i][:-9])
        axs[i].set_xlabel('frequency [Hz]')
        axs_all.append(a)
        
    axs[0].set_ylabel('PSD [V**2/Hz]')
    plt.show()





######## Plotting functions ###########

def plot_image(flux, num_im, hdr):
    
    wcs_dict = get_wcs_coords(hdr)

    for i in range(num_im):
        #plt.xlabel('Column')
        #plt.ylabel('Row')
        #plt.imshow(flux[i], cmap=plt.cm.YlGnBu_r)
        fig = plt.figure(figsize=(10, 10))
        ax = plt.subplot(projection=wcs_dict)
        plt.imshow(flux[i], origin='lower', cmap=plt.cm.YlGnBu_r, aspect='equal')
        plt.xlabel(r'RA')
        plt.ylabel(r'Dec')


        overlay = ax.get_coords_overlay('icrs')
        overlay.grid(color='white', ls='dotted')
        plt.colorbar()
        #plt.clim(0,300)
        plt.title('Flux Image: %s' %(i))
        plt.show()
        
def plot_light_curve(time, flux, flux_err, row_ind, col_ind):

    time_filt, flux_filt, flux_err_filt = filter_parameters(time, flux, flux_err)
    
    star_pix = []
    star_pix_err = []
    star_time = []
    for i in range(len(flux_filt)):
        star_pix.append(flux_filt[i][row_ind][col_ind])
        star_pix_err.append(flux_err_filt[i][row_ind][col_ind])
        star_time.append(time_filt[i])
    
    plt.figure(figsize = (9, 5))
    plt.title('Light Curve Pixel(row, col) = (%s, %s)' %(row_ind, col_ind))
    plt.errorbar(star_time, star_pix, yerr = star_pix_err, marker = '.', markersize = '2', label = 'After NaN & 0 Cuts')
    plt.xlabel('Time', fontsize = 15)
    plt.ylabel('Flux', fontsize = 15)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    #plt.legend()
    plt.show()


def plot_summed_lc_multiple_surveys(time, flux, flux_err, start_ind, file_num, filenames, file_paths, plot = 'yes', savefig = 'no'):

    fsum_all = []
    time_all = []
    fnames = []

    ind = start_ind
    for i in range(start_ind, start_ind+file_num):
        time_filt, fsum_filt = sum_img_flux(time[ind], flux[ind], flux_err[ind])
        fsum_all.append(fsum_filt)
        time_all.append(time_filt)
        fnames.append(filenames[ind])
        ind += int(len(time)/3)

    colors = ['slateblue', 'darkorange', 'pink', 'darkgreen', 'cyan', 'blue']
    plt.figure(figsize = (11,8))
    for i in range(len(fsum_all)):
        plt.plot(time_all[i], fsum_all[i], color = colors[i], alpha = 0.7, label = fnames[i])
    
    plt.xlabel('Time (days)', fontsize = 21)
    plt.ylabel('Flux (e-/s)', fontsize = 21)
    plt.xticks(fontsize = 19)
    plt.yticks(fontsize = 19)
    plt.legend(fontsize = 15)
    if savefig != 'no':
        plt.savefig(os.path.join(file_path_all[0], fnames[0][:-13]+'_full_obs_summed_flux.png'))
    plt.show()

######### execute codes ##########


def execute_lc(filename, row_ind, col_ind, num_im = 1, add_back = 'no'):
    
    time, flux, flux_err, hdr = open_file(filename, add_back)
    plot_image(flux, num_im, hdr)
    plot_light_curve(time, flux, flux_err, row_ind, col_ind)
    
def execute_img_sum(filename, plot = 'no', savefig_filename = 'no', add_back = 'no'):
    
    time, flux, flux_err, hdr = open_file(filename, add_back)
    time_filt, fsum_filt = sum_img_flux(time, flux_err, plot, savefig_filename)
    
    return time_filt, fsum_filt

def execute_all(filename, row_ind, col_ind, num_im = 1, plot = 'no'):
    
    execute_lc(filename, row_ind, col_ind, num_im = 1)
    execute_img_sum(filename, plot)


