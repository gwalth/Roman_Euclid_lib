

import math,glob,sys

from astropy.table import Table
import astropy.io.fits as pyfits
#import astropy.wcs as pywcs
from astropy.wcs import WCS

import numpy as np
from scipy import integrate, interpolate
from scipy import linalg
from scipy.optimize import leastsq, least_squares

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.colors import Normalize


#import grizli_aXe_glw
from grizli_aXe_glw import aXeConf_glw, EuclidData, EuclidFilters, zodi_subtract, insidePolygon, flam_to_fnu
from grizli_functions import euclid_wcs, remove_basic_wcs


import matplotlib as mpl    
mpl.rcParams['font.size'] = 8

# ## highSNR
# [top](#Table-of-Contents)

def MAD(x):
    #return np.median(np.abs(x-np.median(x)))
    return np.nanmedian(np.abs(x-np.nanmedian(x)))

def stats(X):
    x_min = np.nanmin(X)
    x_max = np.nanmax(X)
    x_mean = np.nanmean(X)
    x_median = np.nanmedian(X)
    x_std = np.nanstd(X)
    x_mad = MAD(X)
    x_q23 = np.quantile(X, 0.023)
    x_q159 = np.quantile(X, 0.159)
    x_q841 = np.quantile(X, 0.841)
    x_q977 = np.quantile(X, 0.977)        

    stat_dict = {"min": x_min,
                 "max": x_max,
                 "mean": x_mean,
                 "median": x_median,
                 "std": x_std,
                 "mad": x_mad,
                 "q(2.3)": x_q23,
                 "q(15.9)": x_q159,
                 "q(84.1)": x_q841,
                 "q(97.7)": x_q977}

    print("min =", x_min)
    print("max =", x_max)
    print("mean =", x_mean)
    print("median =",  x_median)
    print("std =", x_std)
    print("std_mad =", x_mad)
    print("q( 2.3) =", x_q23)
    print("q(15.9) =", x_q159)
    print("q(84.1) =", x_q841)
    print("q(97.7) =", x_q977)

    return stat_dict

def gauss_fn(param, x):
    # From Dan's gaussfit
    """Parameters: a0 = height of exp, a1 = center of exp, a2 = sigma
    (the width)"""
    a0, a1, a2, a3 = param
    return a0 * np.exp(-(((x - a1) / a2) ** 2.0) / 2.0) + a3


noise = 1.
err_fn = lambda p, x, z: (gauss_fn(p, x) - z) / noise




simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/Prep/"
raw_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/frame1/"

#gain = 3500.0
gain = 1.0
#aper = 8.
aper = 12.
#aper = 20.

gfit = False

cmap = cm.viridis
#cmap = cm.gist_ncar
#cmap = cm.jet



#sci_file = raw_path  + "NISPS_TIPS_TestPoints_highSNR_mod1_14324_2023_05_26_frame1.fits"
#zodi_file = raw_path  + "NISPS_TIPS_TestPoints_highSNR_JustZodi_14324_2023_05_26_frame1.fits"
#out_file = simdata_path  + "NISPS_TIPS_TestPoints_highSNR_subZodi_14324_2023_05_26_frame1.fits"


# In[14]:


#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=6.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=1.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=2.0)


# In[86]:


#prefix = "/Users/gwalth/data/aXeSIM/Roman/aXeSIM_Roman/"

conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv1/"

#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS000p0/"    # Frame1   
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS180p4/"    # Frame2
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS180p0/"    # Frame3
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS000m4/"    # Frame4


#roman_path = prefix + "roman_setup/"
simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/Prep/"

#conf = "config_axesim_det_11.conf"
#conf = "NISP_RGS180p4_12.conf"
conf = "test_NISP_RGS180p4_12.conf"
sens_file = 'SENS_A.fits'
#input_spec_list = "input_spectra.lis"

direct_file = simdata_path + "Euclid-NISP_H_ref.fits"
#slitless_file = simdata_path + "NISP_SLITLESS_FRAME1_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"    # RGS000p0
slitless_file = simdata_path + "NISP_SLITLESS_FRAME2_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"   # RGS180p4
#slitless_file = simdata_path + "NISP_SLITLESS_FRAME3_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"   # RGS180p0
#slitless_file = simdata_path + "NISP_SLITLESS_FRAME4_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"   # RGS000m4



#idx = 0
#xc = cat[idx]['X_IMAGE']
#yc = cat[idx]['Y_IMAGE']
#mag = cat[idx]['MAG_AUTO']
#spectemp = cat[idx]['SPECTEMP']

#print(xc, yc, mag)

##################################################

aXe = aXeConf_glw(conf_file=conf_path + conf)

print(aXe.conf)
#aXe.get_beam(xc, yc)
#aXe.load_files(direct_file, slitless_file, offset=0)
aXe.load_files(direct_file, slitless_file, offset=-1024.)

aXe.load_sensitivity(conf_path+sens_file)

#idx = 0
#xc = cat[idx]['X_IMAGE']
#yc = cat[idx]['Y_IMAGE']
#mag = cat[idx]['MAG_AUTO']
#spectemp = cat[idx]['SPECTEMP']

#print(xc, yc, mag)

##################################################

aXe = aXeConf_glw(conf_file=conf_path + conf)

print(aXe.conf)
#aXe.get_beam(xc, yc)
#aXe.load_files(direct_file, slitless_file, offset=0)
aXe.load_files(direct_file, slitless_file, offset=-1024., fix_mosa=1)
#aXe.load_files(direct_file, slitless_file, pre_process_path=conf_path, offset=-1024, fix_mosa=1, det='11', remove_dc=1, correct_gain=0, diag=1)

aXe.load_sensitivity(conf_path+sens_file)


# In[47]:


euc_dat = EuclidData(
    catalog1 = "FSpatch_mod3_16183_TESTarrayCALIB_V1.fits",
    catalog2 = "NisInputConfiguration_1_catalog_11.cat",
    model_spectra = "catalog_11_11.1.spc.fits",
    euclud_sims_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/",
    catalog1_path = "FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/",
    catalog2_path = "FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/frame1/Catalogs/",
    model_spectra_path = "FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/frame1/Input_Spectra/",
)


# In[48]:


# slitless is 2048 x 2048
# reference is 10000 x 10000

#catalog_file = simdata_path + "MOT_SCA1_roll_0_dither_0x_0y_cut_zcut_cat14-18.fits"
#catalog_file = simdata_path + "Euclid_total_ref.cat"
catalog_file = simdata_path + "Euclid-NISP_H_ref.cat"

print(catalog_file)
#primer = Table.read(catalog_file, format='fits') 
cat = Table.read(catalog_file, format='ascii.sextractor') 

#print(cat.colnames)
#print(len(cat))


########################################################
# find a source in all of the Euclid detector footprints
########################################################
# bright source near the center of DET11
#a0, d0 = 213.5310504, 57.4451151


a0,d0 = 213.5306062, 57.4122091

dr = np.sqrt((cat['X_WORLD']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (cat['Y_WORLD']-d0)**2)*3600.
id = cat['NUMBER'][np.argmin(dr)]
print(id)
print(np.argmin(dr))
obj_mag = cat['MAG_AUTO'][np.argmin(dr)]
print('ID:%d, mag=%.2f, dr=%.2f"' % (id, obj_mag, np.min(dr)))

#############################################################
# conversion between the reference image and slitless spectra
#############################################################

#ext = 1
#slitless_files = glob.glob(simdata_path + "Euclid_FRAME1_DET11_slitless.fits")

#slitless_files.sort()
#print(slitless_files)

fix_mosa = 1
det = '11'

hdu = pyfits.open(slitless_file)
head_slit = hdu['DET%s.SCI' % (det)].header
print(head_slit)
naxis1 = head_slit["NAXIS1"]
naxis2 = head_slit["NAXIS2"]


if fix_mosa:
    head_wcs = euclid_wcs( head_slit )
else:
    head_wcs = head_slit 
print(head_wcs)

wcs_slit = WCS(head_wcs)
print(wcs_slit)

A = wcs_slit.calc_footprint(axes=(naxis1,naxis2))
print(A)

p0 = [a0, d0]
print(p0)
inside = insidePolygon(A, p0)
print(inside)

if inside:
    pixcoords = wcs_slit.wcs_world2pix([p0],1)
    print(pixcoords)
    xc, yc = pixcoords[0]


#euc_dat.find_source(ra, dec)
#euc_dat.get_spectra_ids()
model_wav, model_flux = euc_dat.get_model_spectra(a0, d0)
euc_dat.plot_model_spectra()


# In[87]:


#aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
aXe.get_beam_invert_cutout(xc, yc, aper, beam='A', offset=0.0, dy0=-300.0, dy1=300.0)

#aXe.extract_cutout()
aXe.extract_cutout(disp_axis=1, spatial_axis=0)

aXe.plot_cutout(vmin=-10., vmax=30., verb=1, flip_axes=1)
aXe.plot_trace(vmin=-10, vmax=30)
#aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)
aXe.get_thumbnail(a0, d0)

aXe.plot_thumbnail(cat=cat, vmin=-1, vmax=1)
#aXe.plot_sources(vmin=-1, vmax=1)
#aXe.plot_sources(
#    x=cat["X_WORLD"], 
#    y=cat["Y_WORLD"], 
#    labels=cat["MAG_AUTO"], 
#    vmin=-1, 
#    vmax=1,
#    plot_all=True,
#    mag_cut=17,
#)

aXe.plot_footprint()

#aXe.plot_sources_wcs(
#    x=subcat["X_WORLD"], 
#    y=subcat["Y_WORLD"], 
#    labels=subcat["MAG_AUTO"], 
#    vmin=-1, 
#    vmax=500,
#    plot_all=True,
#    mag_cut=17,
#)

aXe.plot_extraction(wav_units=False)

aXe.plot_extraction()

aXe.extract_cutout(disp_axis=1, spatial_axis=0, flip_axis=1)
aXe.plot_extraction()

y0 = -1e-17
y1 = 1e-16

#y1 = 1e-13

aXe.plot_extraction(flux_units=True, scale=1./gain, y0=y0, y1=y1)

#aXe.plot_extraction(flux_units=True, scale=1./gain,
#                    model=[model_wav, model_flux])
aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, scale=1./gain,
                    model=[model_wav, model_flux])


# In[23]:


#help(EuclidFilters)
ef = EuclidFilters()
ef.plot_filters()
print(ef.filters.keys())
model_fnu = flam_to_fnu(flam=model_flux, wav=model_wav)
ef.calc_spectra_mag("NISP_H", model_fnu, model_wav/10.)


