from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import math
import astropy.io.fits as pyfits
import pandas as pd
from astropy.wcs import WCS

from astroquery.mast import Mast
from astroquery.mast import Observations


### DEFINITIONS ###

def open_file(filename):
    
    h = pyfits.open(filename)
    data = h[1].data
    hdr = h[1].header
    #print(h.info())
    h.close() # close the file
    
    return data, hdr

def filter_nan(data):
    
    '''remove nan from time, flux, and flux_err as well
    as flux = zero data points'''
    
    time = np.array(data['TIME']).astype(float)
    flux = np.array(data['FLUX']).astype(float)
    flux_err = np.array(data['FLUX_ERR']).astype(float)
    back_flux = np.array(data['FLUX_BKG']).astype(float)
    
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
    back_flux = back_flux[time_isnan]
    
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
    back_flux = back_flux[flux_isnan]
    
    flux_err_isnan = []
    for i in flux_err:
        flux_err_isnan.append(~np.isnan(i[half_row_len][half_col_len]))
    #applying to arrays
    time = time[flux_err_isnan]
    flux = flux[flux_err_isnan]
    flux_err = flux_err[flux_err_isnan]
    back_flux = back_flux[flux_err_isnan]
    
    mean_back_flux_frame = np.nanmean(np.array(back_flux), axis = 0)
    
    return time, flux, flux_err, mean_back_flux_frame
    
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


def plot_image(binarytab, num_im, hdr):
    
    wcs_dict = get_wcs_coords(hdr)
    
    flux = binarytab['FLUX']

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
        
def plot_light_curve(binarytab, row_ind, col_ind):

    time, flux, flux_err, back_flux = filter_nan(binarytab)
    
    star_pix = []
    star_pix_err = []
    star_time = []
    for i in range(len(flux)):
        star_pix.append(flux[i][row_ind][col_ind])
        star_pix_err.append(flux_err[i][row_ind][col_ind])
        star_time.append(time[i])
    
    plt.figure(figsize = (9, 5))
    plt.title('Light Curve Pixel(row, col) = (%s, %s)' %(row_ind, col_ind))
    plt.errorbar(star_time, star_pix, yerr = star_pix_err, marker = '.', markersize = '2', label = 'After NaN & 0 Cuts')
    plt.xlabel('Time', fontsize = 15)
    plt.ylabel('Flux', fontsize = 15)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    #plt.legend()
    plt.show()
    

'''def sum_flux(img_set, plot = 'no', savefig_filename = 'no'):
    
    time, flux, flux_err, back_flux  = filter_nan(img_set)
    
    flux_sum = []
    for i in flux:
        flux_sum_row = []
        for j in i:
            flux_sum_row.append(sum(j))
        flux_sum_row = np.nan_to_num(flux_sum_row, nan = 0.0)
        flux_sum.append(sum(flux_sum_row))
    
    if plot == 'yes':
        plt.figure(figsize = (9, 5))
        #plt.title('Light Curve Pixel(row, col) = (%s, %s)' %(row_ind, col_ind))
        plt.plot(time, flux_sum, color = 'purple')
        plt.xlabel('Time', fontsize = 15)
        plt.ylabel('Sum of Image Flux', fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        if savefig_filename != 'no':
            plt.savefig(savefig_filename)
        plt.show()
        
    #frame_ave_flux_sum = np.array(flux_sum)/(len(flux))
    
    
    return time, flux_sum'''
    
def sum_flux(img_set, plot = 'no', savefig_filename = 'no', add_back = 'no'):
    
    time, flux, flux_err, back_flux = filter_nan(img_set)
    
    back_add_flux = flux+back_flux
    
    if add_back == 'yes':
        flux_sum = []
        for i in back_add_flux:
            flux_sum_row = []
            for j in i:
                flux_sum_row.append(sum(j))
            flux_sum_row = np.nan_to_num(flux_sum_row, nan = 0.0)
            flux_sum.append(sum(flux_sum_row))
    else:
        flux_sum = []
        for i in flux:
            flux_sum_row = []
            for j in i:
                flux_sum_row.append(sum(j))
            flux_sum_row = np.nan_to_num(flux_sum_row, nan = 0.0)
            flux_sum.append(sum(flux_sum_row))
    
    if plot == 'yes':
        plt.figure(figsize = (9, 5))
        #plt.title('Light Curve Pixel(row, col) = (%s, %s)' %(row_ind, col_ind))
        plt.plot(time, flux_sum, color = 'purple')
        plt.xlabel('Time', fontsize = 15)
        plt.ylabel('Sum of Image Flux', fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        if savefig_filename != 'no':
            plt.savefig(savefig_filename)
        plt.show()
    
    #frame_ave_flux_sum = np.array(flux_sum)/(len(flux))
    
    return time, flux_sum
    
def execute_lc(filename, row_ind, col_ind, num_im = 1):
    
    data, hdr = open_file(filename)
    plot_image(data, num_im, hdr)
    plot_light_curve(data, row_ind, col_ind)
    

'''def execute_img_sum(filename, plot = 'no', savefig_filename = 'no'):
    
    data, hdr = open_file(filename)
    time, fsum = sum_flux(data, plot, savefig_filename)
    
    return time, fsum'''
    
def execute_img_sum(filename, add_back = 'no', plot = 'no', savefig_filename = 'no'):
    
    data, hdr = open_file(filename)
    time, fsum = sum_flux(data, plot, savefig_filename, add_back)
    
    return time, fsum

def execute_all(filename, row_ind, col_ind, num_im = 1, plot = 'no'):
    
    execute_lc(filename, row_ind, col_ind, num_im = 1)
    execute_img_sum(filename, plot)
    



