#!/usr/bin/env python
# coding: utf-8

# # Euclid Grizli Pipeline using grizli 1.6.0 (2023.0217)

# ### Euclid Parameters and Requirements
# https://sci.esa.int/web/euclid/-/euclid-nisp-instrument
# 
# https://www.euclid.caltech.edu/page/technical-details

# ### Wide Survey Imaging
# (20,000 deg^2 = 2pi sr)
# 
# RIZ <= 24.5 AB (10 sigma, extended source)
# 
# YJH <= 24 AB (5 sigma)
# 
# ### Deep Survey Imaging
# (40 deg^2)
# 
# RIZ <= 26.5 AB (10 sigma)
# 
# YJH <= 26 AB (5 sigma)
# 
# 
# https://sci.esa.int/c/portal/doc.cfm?fobjectid=46064
# 
# https://www.euclid.caltech.edu/page/62
# 
# https://ui.adsabs.harvard.edu/abs/2022A%26A...662A.112E/abstract

# ## To-Do

# # Table of Contents
# 
# 1. Setup
#     1. [Import Python modules](#Import-Python-modules)
#     2. [Install templates for redshift fitting](#Install-templates-for-redshift-fitting)
#     3. [Python Helper Functions](#Python-Helper-Functions) 
#     4. [Path to my simulation directories](#Path-to-my-simulation-directories)
# 2. [Find sources in the direct image](#Find-sources-in-the-direct-image)
# 3. [Read SExtractor Photometry of Direct Images](#Read-SExtractor-Photometry-of-Direct-Images)
# 4. [Euclid object simulation](#Euclid-object-simulation)
# 5. [Check simulation](#Check-simulation)
# 6. [Show direct images and slitless spectra](#Show-direct-images-and-slitless-spectra)
# 7. [Loop over all objects and fit their redshifts](#Loop-over-all-objects-and-fit-their-redshifts)
# 8. [Inspect Redshift Fitting Results](#Inspect-Redshift-Fitting-Results)
# 9. [Extract a single 2D spectrum](#Extract-a-single-2D-spectrum)
# 10. [1D Spectral Extraction](#1D-Spectral-Extraction)
# 11. [Display Redshift Fit](#Display-Redshift-Fit)
# 
# Appendix - Old
# 1. [aXeSIM predictions based on conf file](#aXeSIM-predictions-based-on-conf-file)
# 2. [Show 2D beam](#Show-2D-beam)
# 3. [Simple SN calculations based on the spcontetc](#Simple-SN-calculations-based-on-the-spcontetc)
# 4. [Simple SN calculations based on the pzcaletc](#Simple-SN-calculations-based-on-the-pzcaletc)
# 5. [Simple SN calculations based on the apttables2021](#Simple-SN-calculations-based-on-the-apttables2021)
# 6. [Roman and Euclid Sensitivity Function](#Roman-and-Euclid-Sensitivity-Function)
# 7. [Velocity resolution](#Velocity-resolution)
# 8. [Fit redshift to source](#Fit-redshift-to-source)
# 9. [Coordinates Check](#Coordinates-check)
# 10. [SED Check](#SED-check)
# 
# [top](#Table-of-Contents)

# # Setup

# In[ ]:


#get_ipython().run_line_magic('matplotlib', 'inline')


# ## Import Python modules
# [top](#Table-of-Contents)

# In[ ]:


# Only use if you are editing the Python code while running the Jupyter notebook.
#import importlib
#importlib.reload(grizli_functions)


# In[ ]:

import time
T_START = time.time()




import grizli_functions
#from grizli_functions import wcs_pixel_scale, check_sims, create_circular_mask
from grizli_functions import add_noise, wcs_pixel_scale, check_sims, display_grizli
from grizli_functions import fake_euclid_direct, fake_euclid_ref, total_image
from grizli_functions import euclid_det, read_slitless_headers 
from grizli_functions import write_individual_slitless, plot_slitless, map_src_to_det
print(grizli_functions.__file__)
# import jwst failure is ok!


# In[ ]:


import glob, os, sys
from collections import OrderedDict

import matplotlib as mpl    
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

from IPython.display import Image

mpl.rcParams['figure.figsize'] = (10.0, 6.0)
mpl.rcParams['font.size'] = 14
mpl.rcParams['savefig.dpi'] = 72

import numpy as np
from scipy import integrate
#from math import cos, sin, atan2, pi
from math import sqrt, log


import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord

import astropy.io.fits as pyfits
from astropy import wcs
from astropy.table import Table, unique, join, vstack
from astropy.modeling import models

import drizzlepac
import photutils

import grizli
import grizli.model
import grizli.multifit
from grizli import utils, multifit, fitting
import grizli.fake_image
from grizli.pipeline import auto_script
from grizli import prep

print('\n Python version: ', sys.version)
print('\n Grizli version: ', grizli.__version__)
print('\n Astropy version: ', astropy.__version__)


# In[ ]:

####################################
# paramters
####################################
plot = 0


####################################
# constants
####################################
Jy = 1.0E-23        # erg/s/cm^2/Hz
mJy = 1e-3          # Jy
uJy = 1e-6          # Jy
Mpc = 3.086e24      # cm
Ang = 1E-8          # cm
mu = 1E-4           # cm

c_km = 2.9979E5     # km/s
c = 2.9979E10       # cm/s
h = 6.626068E-27    # cm^2*g/s           # erg s 
k = 1.3806503E-16   # cm^2*g/(s^2*K)


# ## Install templates for redshift fitting
# [top](#Table-of-Contents)
# 
# Run only once for the install

# #### WFC3 and ACS calibs

# In[ ]:


#grizli.utils.fetch_default_calibs()


# #### WFC3 PSF and Pickles stars

# In[ ]:


#grizli.utils.fetch_config_files()


# #### Templates used in fitting

# In[ ]:


#grizli.utils.symlink_templates(force=False)


# # Python Helper Functions
# [top](#Table-of-Contents)

# In[ ]:


emlines = [["OVI",         1038.0],         # 0
           ["Ly$\\alpha$", 1215.67],        # 1
           ["CIV",     1550.0],             # 2
           ["CIII]",   1909.],              # 3
           ["CII]",    2327.],              # 4
           ["MgII",    2796.4],             # 5
           ["MgII",    2803.5],             # 6
           ["NeV",     3326.],              # 7
           ["[OII]",   3727.],  # O2        # 8
           ["[NeIII]", 3868.7],             # 9
           ["H$\gamma$",  4340.5],  # Hg    # 10
           ["[OIII]",  4363.0],  # O31      # 11
           ["H$\\beta$",   4861.3],  # Hb   # 12
           ["[OIII]",  4959.0],  # O32      # 13
           ["[OIII]",  5007.0],  # O33      # 14
           ["[NII]",   6548.1],             # 15
           ["H$\\alpha$",  6562.8],  # Ha   # 16
           ["[NII]",   6583.0],             # 17
           ["[SII]",   6717.0],             # 18
           ["[SII]",   6731.0],             # 19
           ["P$\\delta$", 10049.8],  # Pd   # 20
           ["P$\\gamma$", 10938.0],  # Pg   # 21
           ["P$\\beta$",  12818.1],  # Pb   # 22
           ["P$\\alpha$", 18750.1],  # Pa   # 23 
           ["Br$\\delta$", 19440.0],  # Br-d (wikipedia, not exact)
           ["Br$\\gamma$", 21660.0],  # Br-g (wikipedia, not exact)
           ["Br$\\beta$",  26250.0],  # Br-b (wikipedia, not exact)
           ["Br$\\alpha$", 40510.0],  # Br-a (wikipedia, not exact) 
          ]

# http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/hydspec.html
# http://articles.adsabs.harvard.edu//full/1934ApJ....80...19M/0000022.000.html


# ## (TESTING) Reload Core Grizli Python Functions - SKIP
# 
# Only use if you are editing the Python code while running the Jupyter notebook.

# In[ ]:


#import importlib
#importlib.reload(grizli.grismconf)


# In[ ]:


#import importlib
#importlib.reload(grizli.pipeline.auto_script)


# In[ ]:


#import importlib
#importlib.reload(grizli.model)


# In[ ]:


#import importlib
#importlib.reload(grizli.multifit)


# In[ ]:


#import importlib
#importlib.reload(grizli.fitting)


# In[ ]:


#import importlib
#importlib.reload(grizli.utils)


# # Start processing here!

# ## Path to my simulation directories
# [top](#Table-of-Contents)

# In[ ]:


#os.chdir('../')
#os.chdir('/Users/gwalth/data/Roman/grizli/sims/')
os.chdir('/Users/gwalth/data/Roman/grizli/sims/Euclid')

#os.chdir('/local/RomanSims/grizli/sims/') # cygnusd
HOME_PATH = os.getcwd()
print('HOME_PATH = ', HOME_PATH)
#root = "SIM_10_18_22"
#root = "SIM_12_23_22"
root = "TestPoints"


# ## Input files

# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))


# In[ ]:


#slitless_files = glob.glob('EUC_SIM_*.fits')
#slitless_files.sort()
#print(len(slitless_files))

# First test
slitless_files = ['NISPS_TIPS_TestPoints_highSNR_mod1_14324_2023_05_26_frame1.fits']
catalog_files = ['TestPoints_highSNR_mod1_14324.fits']

print(slitless_files)
print(catalog_files)


# In[ ]:


#slitless_files = glob.glob("EUC_SIM_NISR*.fits")
#catalog_files = glob.glob("CATALOG*.fits")
#print(slitless_files)
#print(catalog_files)


# ## Directory Structure
# 
# I was structing it similar to Grizli with Prep, RAW and Extraction directories.  If this were real mission data, the stage that we recieved from Anihita would have been drizzled images and spectra which would go into the Prep directories.
# 
# This is just showing that we have the right directories and we can find all of the files.

# In[ ]:


# Clean
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

os.system('rm *_direct.fits')
os.system('rm *_slitless.fits')
os.system('rm *_wcs.fits')
os.system('rm *_final*.fits')
os.system('rm *_final.cat')
os.system('rm Euclid_GrismFLT.pickle')


# In[ ]:


# detector header names
all_det = euclid_det()


# In[ ]:


## Write individual files for each extension of the slitless spectra
# Grizli is easier to manage when writing out all of the files. 
# At some point we'll want to read the data extensions directly into Grizli, 
# this is currently a kludge.
all_slitless = [write_individual_slitless(sf) for sf in slitless_files]
print(all_slitless)


# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

if plot:
    file = "Euclid_DET11_slitless.fits"
    pf = pyfits.open(file)

    data = pf['SCI'].data - 1024.
    X = data.flatten()

    data = pf['ERR'].data
    Y = data.flatten()

    fig = plt.figure()

    ax1 = fig.add_subplot(121)
    ax1.hist(X,bins=20, range=(0,2000))
    ax1.set_yscale("log")
    ax1.set_xlabel("Counts")
    ax1.set_title("SCI")

    ax2 = fig.add_subplot(122)
    ax2.hist(Y,bins=20, range=(0,100))
    ax2.set_yscale("log")
    ax2.set_xlabel("Counts")
    ax2.set_title("ERR")

    plt.show()


# In[ ]:


# Plot the slitless extensions
if plot:
    plot_slitless(slitless_files[0], vmin=500, vmax=1700, verb=0)


# In[ ]:


## Read slitless headers and plot the image coordinates of the detectors relative to each other
heads = read_slitless_headers(slitless_files[0], verb=0, plot=0)
print(heads)


# ## Read the source catalog

# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

print(catalog_files)
primer = Table.read(catalog_files[0]) 
print(primer.colnames)
print(len(primer))


# In[ ]:


filt = primer['RA'] < 200.
print(filt)
print(primer[filt])

index = np.arange(len(primer))
print(index[filt])
primer.remove_row(index[filt][0])

#filt = primer['VIS'] == -99
#print(primer[filt])

#filt = primer['NIR_Y'] == -99
#print(primer[filt])

#filt = primer['NIR_J'] == -99
#print(primer[filt])

#filt = primer['NIR_H'] == -99
#print(primer[filt])


# In[ ]:


print([col for col in primer.colnames if "TU_" in col])
Euclid_bands = ['VIS','NIR_Y','NIR_J','NIR_H']
Euclid_bands_flux = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG'] 


# In[ ]:


for bk,fk in zip(Euclid_bands,Euclid_bands_flux):
    fnu_Jy = primer[fk] # Jy    
    mab = -2.5*np.log10(fnu_Jy) + 8.90   # Jy --> AB mag            
    mab[np.isinf(mab)]=-99.
    primer[bk] = mab    


# In[ ]:


#primer[:10].show_in_notebook()


# In[ ]:


#primer[Euclid_bands][:10].show_in_notebook()


# In[ ]:


#primer[Euclid_bands_flux][:10].show_in_notebook()


# ## Plot the magnitude histogram of sources

# In[ ]:


primer['TU_FNU_VIS_MAG'] # Jy
primer['TU_FNU_Y_NISP_MAG'] # Jy
primer['TU_FNU_J_NISP_MAG'] # Jy
primer['TU_FNU_H_NISP_MAG'] # Jy
primer['VIS'] # mag
primer['NIR_Y'] # mag
primer['NIR_J'] # mag
primer['NIR_H'] # mag

if plot:
    fig = plt.figure(figsize=(7,5))

    ax1 = fig.add_subplot(111)
    ax1.hist(primer['VIS'], bins=50, range=(15,30), alpha=0.6, label="VIS")
    ax1.hist(primer['NIR_Y'], bins=50, range=(15,30), alpha=0.6, label="NIR_Y")
    ax1.hist(primer['NIR_J'], bins=50, range=(15,30), alpha=0.6, label="NIR_J")
    ax1.hist(primer['NIR_H'], bins=50, range=(15,30), alpha=0.6, label="NIR_H")
    ax1.set_xlabel("mag")
    ax1.legend()

    plt.show()


# ## Plot the RA and Dec distribution across the detectors

# In[ ]:

if plot:
    fig = plt.figure()

    ax1 = fig.add_subplot(111)
    ax1.scatter(primer['RA'],primer['DEC'],s=0.01, alpha=1.0, edgecolors="none")
    ax1.set_aspect(1.)
    ax1.set_xlabel("RA [deg]")
    ax1.set_ylabel("Dec [deg]")

    plt.show()


# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))
direct_thumbnail_files = glob.glob('NIS_catalog_file_??.thm.beamA.fits')
direct_thumbnail_files.sort()
print(direct_thumbnail_files)
print(len(direct_thumbnail_files))


# In[ ]:


#os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))
## Map the sources the sources from the catalog to the detectors 
#det_tbl, det_dict = map_src_to_det(plot=0)
#det_dict = map_src_to_det(primer, heads, plot=0)


# ## Print each of the detectors WCS header information

# In[ ]:




#os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
#wcs_temp = "Euclid_DET%s_wcs.fits"

#for i in range(0,4): 
#    for j in range(0,4):
#        print(i+1,j+1)
#        num = '%i%i' % (i+1,j+1)
#        hdu = pyfits.open(wcs_temp % (num), ignore_missing_simple=True)    
#        head = hdu[0].header
#        #print(head)
#        print(head['CRVAL1'],head['CRVAL2'],head['CRPIX1'],head['CRPIX2'])


# In[ ]:


#primer["TU_FNU_H_NISP_MAG","NIR_H"][:10].show_in_notebook()


# In[ ]:


#    # effective wavelength (pivot wavelength)
#    # http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/
#    numerator = integrate.trapz(tran,wav)
#    denominator = integrate.trapz(tran/wav**2,wav)
#    effwav_sq = numerator/denominator
#    effwav = np.sqrt(effwav_sq)



# ## Create input table for the reference image

# In[ ]:


thumb_temp = 'NIS_catalog_file_%s.thm.beamA.fits'
os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))

source_ids = []
thumbs = []

t0 = time.time()

k = 0
for i in range(0,4): 
    for j in range(0,4):
        print(i+1,j+1)
        num = '%i%i' % (i+1,j+1)
        hdus = pyfits.open(thumb_temp % (num))
        #source_ids_0 = [h.header['EXTNAME'] for h in hdus[1:]]
        tmp_source_ids = []
        for h in hdus[1:]:
            tmp_source_ids.append(h.header['EXTNAME'])
        print("N =", len(tmp_source_ids))
        source_ids += tmp_source_ids
        
        #thumbs += [h.data for h in hdus[1:]]
        for h in hdus[1:]:
          
            thumbs.append(h.data)

            sys.stdout.write("Progress: %d   \r" % (k+1))
            sys.stdout.flush()

            k += 1
        print()

t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))

        
#print(source_ids)
print(len(source_ids))

#print(set(source_ids))
print(len(set(source_ids)))

print(sys.getsizeof([]),"bytes")
print(sys.getsizeof([1]),"bytes") 
print("source_ids =",sys.getsizeof(source_ids),"bytes")   
print("thumbs =",sys.getsizeof(thumbs),"bytes")   

src_tbl = Table([source_ids,thumbs], names=['SOURCE_ID','THUMB'])
all_tbl = join(src_tbl, primer, keys='SOURCE_ID')
final_tbl = unique(all_tbl, keys='SOURCE_ID') # look for duplicate entries

print()
print()
print()
print("Size of table =", len(all_tbl))
print("Size of table (no duplicates) =", len(final_tbl))
print()
print()
print()
#world = [[ra,dec] for ra,dec in det_tbl['RA','DEC']]
#world = [[ra,dec] for ra,dec in det_tbl['RA','DEC']]
ra = final_tbl['RA']
dec = final_tbl['DEC']

ra_avg = np.mean(ra)
dec_avg = np.mean(dec)
print("RA (avg) =", ra_avg)
print("Dec (avg) =", dec_avg)


# In[ ]:


print()
print()
print()
print("source_ids =",sys.getsizeof(source_ids),"Bytes")   
print("src_tbl =",sys.getsizeof(src_tbl),"Bytes")   
print("all_tbl =",sys.getsizeof(all_tbl),"Bytes")   
print("final_tbl =",sys.getsizeof(final_tbl),"Bytes")   
print()
print()
print()

# In[ ]:


#all_tbl[:10].show_in_notebook()


# In[ ]:


#print(primer)
print(primer.colnames)


# ## Create reference images

# In[ ]:


# Mauri et al. 2020
# (ropper et al. 2016
os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))

# noise characteristics of "direct image"
# Scaramella et al. 2022 - Euclid preparation I. The Euclid Wide Survey
# RGS000, RGS180, RGS000_rot, RGS180_rot
spec_exptime = 574 # seconds
spec_gain = 2.0


# In[ ]:


nexp = 4 

# VIS 
readnoise = 4.5  # e-
background = 1.2 # e-/s
#pixel_size = 12. # um
dir_gain = 1.0   # e-/photon

# VIS
flux_key = "TU_FNU_VIS_MAG"
mag_key = "VIS"
dir_exptime = 570 # seconds
wav_cen = 7102.613 # calculated from pivot
wav_width = 9000. - 5600.
eff_tot = 0.70
#eff_tot = 0.74
output = "Euclid_VIS_ref.fits"
filt = "VIS"
instr = "VIS"

t0 = time.time()

print("Creating VIS reference image")
test_fluxes_vis, test_mags_vis = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, 
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output, 
                                         filt=filt, instr = instr)



t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))
# In[ ]:


t0 = time.time()


#dir_exptime = 448
#nexp = 1
nexp = 4

# NISP
readnoise = 6.   # e-
background = 1.2 # e-/s

#background = 0 # e-/s
#readnoise = 0
dir_gain = 1.0   # e-/photon

# Y
flux_key = "TU_FNU_Y_NISP_MAG"
mag_key = "NIR_Y"
dir_exptime = 112. # seconds
wav_cen = 10809. # Ang
wav_width = 2627. # Ang
eff_tot = 0.772
output = "Euclid_NISP_Y_ref.fits"
filt = "NISP_Y"
instr = "NISP"

print("Creating NISP_Y reference image")
test_fluxes_nisp_y, test_mags_nisp_y = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, 
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output,
                                         filt=filt, instr = instr)

t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))
# In[ ]:


t0 = time.time()
# J
flux_key = "TU_FNU_J_NISP_MAG"
mag_key = "NIR_J"
dir_exptime = 112. # seconds
wav_cen = 13673. # Ang
wav_width = 3994. # Ang
eff_tot = 0.790
output = "Euclid_NISP_J_ref.fits"
filt = "NISP_J"
instr = "NISP"

print("Creating NISP_J reference image")
test_fluxes_nisp_j, test_mags_nisp_j = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, 
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output,
                                         filt=filt, instr = instr)
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))
# In[ ]:


t0 = time.time()
# H
flux_key = "TU_FNU_H_NISP_MAG"
mag_key = "NIR_H"
wav_cen = 17714. # Ang
wav_width = 4999. # Ang
dir_exptime = 112. # seconds
eff_tot = 0.782
output = "Euclid_NISP_H_ref.fits"
filt = "NISP_H"
instr = "NISP"

print("Creating NISP_H reference image")
test_fluxes_nisp_h, test_mags_nisp_h = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, 
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output,
                                         filt=filt, instr = instr)
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))
# In[ ]:



#print(all_direct)
# should label direct with "_direct"


# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

t0 = time.time()
ref_files = ["Euclid_VIS_ref.fits", "Euclid_NISP_Y_ref.fits", "Euclid_NISP_J_ref.fits",
             "Euclid_NISP_H_ref.fits"]
total_image(ref_files, output="Euclid_total_ref.fits", img_ext='REF') 
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))




# In[ ]:


# remove segmentation FITS to redo "next" step
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
os.system('rm *_seg.fits')


# In[ ]:


## Make SExtractor catalog
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

# Euclid (Schirmer et al. 2022, XVIII. The NISP photometric system)
#mag_zero = 25.04 # +- 0.05 mag   Y 
#mag_zero = 25.26 # +- 0.05 mag   J
#mag_zero = 25.21 # +- 0.05 mag   H

#           VIS   Y      J      H      total
mag_zero = [25.6, 25.04, 25.26, 25.21, 26.0]

det_img = "Euclid_total_ref.fits"

all_cat = []
all_seg = []
all_bkg = []

all_ref_files = ref_files + [det_img]

for i,ref in enumerate(all_ref_files):

    prefix = ref.replace(".fits","")
    
    sex = "Euclid.sex"
             
    cat = prefix + ".cat" 
    #wht = prefix + "_wht.fits"
    seg = prefix + "_seg.fits"
    bkg = prefix + "_bkg.fits"
    #aper = prefix + "_aper.fits"
    
    all_cat.append(cat)
    all_seg.append(seg)
    all_bkg.append(bkg)
        
    if not os.path.exists(seg):
    
        det_wht = det_img + "[2]"
        wht = ref + "[2]"
        det_ext = det_img + "[1]"
        ref_ext = ref + "[1]"

        
        checkimage_name = seg + "," + bkg

        # detection image
        sex_str = 'sex ' + det_ext + "," + ref_ext + ' -c ' + sex + \
                  ' -WEIGHT_IMAGE ' + det_wht + "," + wht + \
                  ' -CHECKIMAGE_NAME ' + checkimage_name + ' -CATALOG_NAME ' + cat + \
                  ' -MAG_ZEROPOINT %.2f' % (mag_zero[i])
        
        # single image
        #sex_str = 'sex ' + ref_ext + ' -c ' + sex + ' -WEIGHT_IMAGE ' + wht + \
        #          ' -CHECKIMAGE_NAME ' + checkimage_name + ' -CATALOG_NAME ' + cat + \
        #          ' -MAG_ZEROPOINT %.2f' % (mag_zero[i])
        print(sex_str)
        os.system(sex_str)
            
    else:
        print("Skipping...")
        

# awk '{ printf "circle(%f, %f, 0.00007) # text={%.3f}\n", $4, $5, $42 }' GRS_FOV1_roll0_dx0_dy0_SCA1_direct.cat > GRS_FOV1_roll0_dx0_dy0_SCA1_direct.reg
        


# In[ ]:

if plot:
    test_fluxes = [test_fluxes_vis, test_fluxes_nisp_y, test_fluxes_nisp_j, test_fluxes_nisp_h]
    test_mags = [test_mags_vis, test_mags_nisp_y, test_mags_nisp_j, test_mags_nisp_h]
    
    for i,zp in enumerate(mag_zero[:-1]):
        
        mag = -2.5 * np.log10(test_fluxes[i]) + zp
    
    
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
        ax.scatter(mag,test_mags[i])
        ax.plot([9,25],[9,25],c="tab:red")
        
        #ax.set_title(cat)
        ax.set_xlabel("ZP derived magnitudes")
        ax.set_ylabel("input magnitudes")
    
        plt.show()
    
    
    

all_direct = all_ref_files


# In[ ]:


dir()
print(all_direct)
print(all_slitless)


# ## Add all of the header metadata needed for Grizli

# In[ ]:


print(dir_exptime)
print(spec_exptime)
print(nexp)


# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

suffix = "_final"

all_final_direct = []
all_final_slitless = []


# ## Subtract 1024 from science and sqrt(error), and adjust header keys

# In[ ]:


###############
### SPECTRA ###
###############
for slitless in all_slitless[0]:
    
    new_slitless = slitless.replace(".fits",suffix+".fits")
    all_final_slitless.append(new_slitless)
        
    hdu = pyfits.open(slitless)
    hdu.info()

    ext = 0
    #hdu[ext].header['INSTRUME'] = 'NISP' 
    hdu[ext].header['INSTRUME'] = 'NISP-GLWv1' 
    # v2:
    # - optical model looks the same for each detector
    # - sensitivity for each detector (are they different?)
    hdu[ext].header['FILTER'] = 'RED'
    hdu[ext].header['EXPTIME'] = spec_exptime

    ext = 1
    #hdu[ext].header['INSTRUME'] = 'NISP'
    hdu[ext].header['INSTRUME'] = 'NISP-GLWv1'
    hdu[ext].header['FILTER'] = 'RED'
    hdu[ext].header['EXTVER'] = ext
    hdu[ext].header['EXPTIME'] = spec_exptime

    sci = hdu[ext].data
    #hdu[ext].data = sci/spec_exptime/gain
    hdu[ext].data = (sci-1024.)/spec_exptime/spec_gain


    ext = 2
    chi2 = hdu[ext].data
    hdu[ext].data = np.sqrt(chi2)/spec_exptime/spec_gain

    hdu.writeto(new_slitless, overwrite=True, output_verify='fix')
    print("Writing",new_slitless)


# In[ ]:


print(all_final_direct)
print(all_final_slitless)


# ## Check the pixel scale of the images

# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

#for i in range(4):
#    for j in range(4):
#        num = '%i%i' % (i+1,j+1)
#        new_direct = "Euclid_DET%s%s.fits" % (num,suffix)
#        new_slitless = "Euclid_DET%s_slitless%s.fits"  % (num,suffix)

for direct in all_final_direct:        
    wcs_pixel_scale(direct)

for slitless in all_final_slitless:
    wcs_pixel_scale(slitless)


# In[ ]:


all_phot = []

for cat in all_cat:
    #print(cat)
    phot = Table.read(cat, format='ascii.sextractor') # ref_cat in multimission
    all_phot.append(phot)
    print(cat," number_of_sources =",len(phot))
print()
print(phot.colnames)


# In[ ]:

if plot:
    fig = plt.figure(figsize=(12,12))
    
    k = 0
    ax = []
    for i,phot in enumerate(all_phot):
        
        cat = all_cat[i]
    
        #print(phot["MAG_AUTO"])
        #print(phot["FLUX_AUTO"])
    
        if len(phot) > 0:    
            print(np.min(phot["MAG_AUTO"]),np.max(phot["MAG_AUTO"]))
            #print(np.min(phot["FLUX_AUTO"]),np.max(phot["FLUX_AUTO"]))
    
        
        ax = fig.add_subplot(4,4,k+1)
        ax.hist(phot["MAG_AUTO"],range=(10,32),bins=44)
        #ax.hist(phot["MAG_AUTO"],bins=20)
        #ax.hist(phot["FLUX_AUTO"],bins=20)
    
    
        ax.set_xlabel(cat)
        ax.set_ylabel("N")
        ax.set_xlim(10,32)
        
        k += 1
        
    plt.show()
    
    
# In[ ]:


for phot in all_phot:
    print(phot)


# In[ ]:

if plot:
    plot_labels = 1
    mag_cut = 17.
    
    fig = plt.figure(figsize=(12,12))
    
    k = 0
    #for direct,phot in zip(all_direct,all_phot):
    for direct,phot in zip(all_ref_files,all_phot):
    
    
        hdu = pyfits.open(direct)
        img = hdu[1].data
        
        ax1 = fig.add_subplot(4,4,k+1)
        ax1.imshow(img,origin="lower",vmin=-0.1,vmax=0.1)
        ax1.scatter(phot["X_IMAGE"],phot["Y_IMAGE"],fc="None",ec="r",s=60,lw=0.1, alpha=0.1)
    
        if plot_labels:
            for j,src in enumerate(phot['MAG_AUTO']):
            #for j,src in enumerate(phot['NUMBER']):
            #for j,src in enumerate(phot['FLUX_AUTO']):
                if src < mag_cut:
                    ax1.text(phot['X_IMAGE'][j],phot['Y_IMAGE'][j],src)
                
        
        k += 1
    plt.show()
    
    



verb = 1
search_rad = 0.6

if plot:
    fig = plt.figure(figsize=(12,12))

all_phot_matched_clean = []
    
k = 0
for i,phot in enumerate(all_phot):
    
    #print(all_final_direct[i])
    print("Total number of sources found =",len(phot))
    print("Search primer and all_tbl for RA/DEC matches")
    
    if len(phot) > 0:
        c_prime = SkyCoord(ra=primer["RA_MAG"]*u.degree, dec=primer["DEC_MAG"]*u.degree)
        c_phot = SkyCoord(ra=phot["X_WORLD"], dec=phot["Y_WORLD"])

        #idx, d2d, d3d = c_prime.match_to_catalog_sky(c_phot)
        idx, d2d, d3d = c_phot.match_to_catalog_sky(c_prime)

        filt = d2d < 1*u.arcsec

        if plot:
            p1 = fig.add_subplot(4,4,k+1)
            p1.hist(d2d[filt].value*3600.,bins=25)
        
        #primer['idx'] = idx
        #primer['d2d'] = d2d

        phot['idx'] = idx
        phot['d2d'] = d2d
        #print(primer.colnames)
        
        #phot['idx'] = np.arange(len(phot))
        primer['idx'] = np.arange(len(primer))
        #print(all_tbl.colnames)
        #print(all_tbl['idx'])
        
        print("Join all_tbl and primer tables")
        match_tbl = join(phot, primer, keys='idx')
        #print(match_tbl)
        #print(match_tbl.colnames)
        
        print("Select only sources < %.2f arcsec" % (search_rad))
        #filt = match_tbl['d2d'] < 1.0*u.arcsec
        clean_filt = match_tbl['d2d'] < search_rad*u.arcsec
        match_clean_tbl = match_tbl[clean_filt]
        print(list(match_clean_tbl['idx']))
        print("Number of matches =",len(match_clean_tbl))
        print()
        
        if not i: print(match_clean_tbl.colnames)
        
    else:
        match_clean_tbl = []
    
    all_phot_matched_clean.append(match_clean_tbl)
    
    k+=1

if plot:
    plt.show()
        


# In[ ]:


if plot:
    plot_labels = 0
    
    fig = plt.figure(figsize=(12,12))
    
    k = 0    
    #for direct,match_clean_tbl in zip(all_direct,all_phot_matched_clean):
    for direct,match_clean_tbl in zip(all_ref_files,all_phot_matched_clean):
    
    
        hdu = pyfits.open(direct)
        img = hdu[1].data
    
        ax1 = fig.add_subplot(4,4,k+1)
        ax1.imshow(img,origin="lower",vmin=-0.1,vmax=0.1)
        #p1.scatter(phot["X_IMAGE"],phot["Y_IMAGE"],fc="None",ec="r",s=60)
        if len(match_clean_tbl) > 0:
            ax1.scatter(match_clean_tbl['X_IMAGE'],match_clean_tbl['Y_IMAGE'],lw=0.3,fc="None",ec="r",s=60)
        
        if plot_labels and len(match_clean_tbl) > 0:
            for j,src in enumerate(match_clean_tbl['SOURCE_ID']):
                #print(j)
                ax1.text(match_clean_tbl['X_IMAGE'][j],match_clean_tbl['Y_IMAGE'][j],src)
        k += 1
    plt.show()
    
    
# In[ ]:


print([col for col in primer.colnames if "TU_" in col])
Euclid_bands = ['VIS','NIR_Y','NIR_J','NIR_H']
Euclid_bands_flux = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG']


# In[ ]:


print("N =",len(all_phot_matched_clean))
print()

verb = 0

for i,match_clean_tbl in enumerate(all_phot_matched_clean[:-1]):  
    zps = []
    mags = []
    
    if len(match_clean_tbl) > 0:
        print("i =",i)
                
        for k,row in enumerate(match_clean_tbl):
            
            #ind = np.argmax(match_clean_tbl['FLUX_AUTO'])

            #print("ind =",ind)
            dn = row['FLUX_AUTO']
            mag = row[Euclid_bands[i]]
            flux = row[Euclid_bands_flux[i]]

            ## mag = -2.5 * log10(dn/exptime) + zp
            zp = mag + 2.5*np.log10(dn)
            
            
            mags.append(mag)
            zps.append(zp)
            if verb:
                print(k)
                print("dn =", dn)
                print("mag =", mag)
                print("flux =", flux)
                print("zp =", zp)
                print()

    print("Magnitudes:")
    print("min = %.3f" % (np.nanmin(mags)))
    print("max = %.3f" % (np.nanmax(mags)))
    print()

    print("Zeropoint:")
    print("mean   = %.3f" % (np.nanmean(zps)))
    print("median = %.3f" % (np.nanmedian(zps)))
    print("std    = %.3f" % (np.nanstd(zps)))
    print("min    = %.3f" % (np.nanmin(zps)))
    print("max    = %.3f" % (np.nanmax(zps)))


    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(mags,zps)
        plt.show()


# In[ ]:


rms = 0.15
zp = 25.13
mag = - 2.5*np.log10(5*rms) + zp
print(mag)


# ## Euclid object simulation
# [top](#Table-of-Contents)

# In[ ]:


#import importlib
#importlib.reload(grizli.model)


# In[ ]:


#import grizli.model


# In[ ]:


print(all_final_slitless)
print(all_ref_files)
print(all_seg)
#print(all_phot_matched_clean)
print(all_phot_matched_clean[-1]['NUMBER','X_IMAGE','Y_IMAGE','MAG_AUTO'])

T_END = time.time()

print()
print("Finished in %.1f seconds" % (T_END-T_START))
sys.exit()

# In[ ]:


print(all_ref_files[3])
print(all_seg[-1])
print(all_phot_matched_clean[-1])


# In[ ]:


#import importlib
#importlib.reload(grizli.grismconf)


# In[ ]:


import grizli.grismconf


# In[ ]:


###################
# Reference image #
###################

T0 = time.time()

os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

# allow simulation of objects at the edges
pad=200 # pixels
#pad = 800 # I think this may be optimal given the spectra size (only seems to add 10-20 spectra)

mag_limit = 30 

all_euclid = []

#for i in range(len(all_final_slitless)):
for i in range(len(all_final_slitless)):
#for i in [0]:

    
    t0 = time.time()

    Euclid = grizli.model.GrismFLT(grism_file=all_final_slitless[i], verbose=True, pad=pad,  
                                   ref_file=all_ref_files[3], ref_ext=1, # NISP_H
                                   seg_file=all_seg[-1], # total
                                   #shrink_segimage=False)
                                   shrink_segimage=True)

    Euclid_cat = Euclid.blot_catalog(all_phot_matched_clean[-1], sextractor=True) 
    #Euclid_cat = Euclid.blot_catalog(all_phot[i], sextractor=True) 
    Euclid.catalog = Euclid_cat

    mask = Euclid_cat['MAG_AUTO'] < mag_limit
    print('N=%d' %(mask.sum()))
    #Euclid.compute_full_model(verbose=True)
    #Euclid.compute_full_model(compute_beams=['A'], mask=mask, verbose=False)
    Euclid.compute_full_model(ids=Euclid_cat['NUMBER'][mask], mags=Euclid_cat['MAG_AUTO'][mask], verbose=True)
    #Euclid.compute_full_model(ids=Euclid_cat['NUMBER'][mask], mags=Euclid_cat['MAG_AUTO'][mask], verbose=True)

    all_euclid.append(Euclid)
    
    t1 = time.time()
    
    print("Detector finished in %.1f seconds" % (t1-t0))
    
    
T1 = time.time()

print()
print("Finished in %.1f seconds" % (T1-T0))


# In[ ]:


import importlib
importlib.reload(grizli.model)


# In[ ]:


import grizli.model


# In[ ]:


print(all_euclid)


# In[ ]:


t0 = time.time()

print(all_euclid)

import pickle

# ~ 5GB file
with open('Euclid_GrismFLT.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(all_euclid, f, pickle.HIGHEST_PROTOCOL)
    
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))


# In[ ]:


t0 = time.time()

import pickle

with open('Euclid_GrismFLT.pickle', 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    all_euclid = pickle.load(f)
    
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))


# In[ ]:


print(all_euclid)
print(type(all_euclid))
print(len(all_euclid))
#print(all_euclid[0].__dict__)


# ## Check simulation
# [top](#Table-of-Contents)
# 
# Possible problem with the aXeSIM conf file.  This function checks if the OrderDict exists within GrismFLT class.

# In[ ]:


mag_limit = 30


# In[ ]:


all_euclid_srcs = []

for i,euclid in enumerate(all_euclid):
    #Euclid_all,Euclid_magcut,Euclid_extract = check_sims(euclid, mag_limit)
    #print(euclid.catalog["MAG_AUTO"])
    #print(euclid.object_dispersers)
    #print(len(euclid.object_dispersers))
    #print(euclid.object_dispersers[93])
    #############################################
    # These are the actual Source Extractor IDs #
    #############################################
    #for id in euclid.object_dispersers:
        #print(id)
        #print(id, euclid.object_dispersers[id])
        #print(len(euclid.object_dispersers[id]))
    #print([(j,id) for j,id in enumerate(euclid.object_dispersers)])
    #sim.object_dispersers.index(id)
    #print(all_det[i])
    all_euclid_srcs.append(check_sims2(euclid, mag_limit))


# In[ ]:


print(all_euclid_srcs)


# In[ ]:


# Euclid_extract
det11 = all_euclid_srcs[0][2]

fig = plt.figure(figsize=(4,4))
ax1 = fig.add_subplot(111)
ax1.hist(det11["MAG_AUTO"],bins=15)
plt.show()

filt = det11["MAG_AUTO"] < 17.0
print(det11["NUMBER","ra","dec","MAG_AUTO"][filt])


# In[ ]:


r0,d0 = 228.68820022, 6.3239380205 # mag=16.7511


# ## Extract a single 2D spectrum 
# Science, model and contamination based on RA and Dec
# 
# [top](#Table-of-Contents)

# In[ ]:


det_ind = 0 # DET11
Euclid = all_euclid[det_ind]
phot = all_phot[det_ind]

## Find the object ID near coordinates (r0,d0), 
## SExtractor IDs might not be constant across platforms

dr = np.sqrt((phot['X_WORLD']-r0)**2*np.cos(d0/180*np.pi)**2 + 
             (phot['Y_WORLD']-d0)**2)*3600.
id = phot['NUMBER'][np.argmin(dr)]
obj_mag = phot['MAG_AUTO'][np.argmin(dr)]
print('ID:%d, mag=%.2f, dr=%.2f"' %(id, obj_mag, np.min(dr)))

beams = OrderedDict()

ix = Euclid.catalog['id'] == id
x0, y0 = Euclid.catalog['x_flt'][ix][0], Euclid.catalog['y_flt'][ix][0]
print(Euclid.direct.instrument, x0, y0)
#print(Roman.wcs.pscale)
#dim = 18*0.135/sim.flt_wcs.pscale 
#beam = grizli.model.BeamCutout(id=id, x=x0, y=y0, 
#                               cutout_dimensions=np.cast[int]((dim, dim)), 
#                               conf=sim.conf, GrismFLT=sim)

print(Euclid.object_dispersers[id])

is_cgs, spectrum_1d, b = Euclid.object_dispersers[id]
#print(b)

cutout = grizli.model.BeamCutout(Euclid, b['A'], min_sens=0,) # min_mask=0) 

cutout.beam.compute_model()  
cutout.contam = cutout.beam.cutout_from_full_image(Euclid.model)
if id in Euclid.object_dispersers:
    cutout.contam -= cutout.beam.model

print(dir(cutout.beam))
#print(cutout.contam)
#print(cutout.beam.model)
#print(dir(cutout.beam))
print(cutout.beam.seg)
print(cutout.beam.seg.shape)

img = cutout.grism.data['SCI']*1
X = img.flatten()
std = np.std(X)
med = np.median(X)

sig = 1.
vmin = med-sig*std
vmax = med+sig*std

fig = plt.figure()
p1 = fig.add_subplot(111)
p1.imshow(cutout.grism.data['SCI']*1, origin='lower',cmap='gray_r',vmin=vmin,vmax=vmax)
p1.set_xlabel("X [pixels]")
p1.set_ylabel("Y [pixels]")

#p1.imshow(cutout.grism.data['ERR']*1)#,vmin=-0.1,vmax=0.1, origin='lower',cmap='gray_r')
    
beams[Euclid.grism.instrument] = cutout

cutout.write_fits() # still learning about the output


# ## 1D Spectral Extraction
# [top](#Table-of-Contents)

# In[ ]:


print(beams.keys())
key = "NISP-GLWv1"
print(beams[key].grism.data.keys())


# In[ ]:


from scipy import integrate

### Plot 1D spectra

fig = plt.figure(figsize=(6,4))
p1 = fig.add_subplot(111)
key = "NISP-GLWv1"
w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], bin=0)
#w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], 
#                                          ivar=1./(beams[key].grism.data['ERR'])**2,bin=0)

#sf = f/np.nanmax(f)
sf = f
std = np.nanstd(sf)
med = np.nanmedian(sf)

sig = 2.
y0 = med-sig*std
y1 = med+sig*std

p1.text(0.05,0.9,"ID = %i" % (id),transform=p1.transAxes)
p1.plot(w/1.e4, sf, c="k") # linestyle='steps-mid')

#p1.plot(w/1.e4, e, c="r") # linestyle='steps-mid')
p1.set_ylabel("Flux [Arbitrary]")
p1.set_xlabel("Wavelength [um]")

p1.set_ylim(y0,y1)


# In[ ]:


def ab2flux(mab,eff_wav):
    c = 2.9979E10 # cm/s
    Ang = 1E-8    # cm
    Jy = 1E-23    # erg/s/cm^2/Hz
    
    # mab = -2.5*np.log10(fnu) - 48.6 
    fnu = 10**(-0.4*(mab + 48.6))          # erg/s/cm^2/Hz
    flambda = fnu*(c/(eff_wav*Ang)**2)*Ang # erg/s/cm^2/Ang
   
    #print("%.2f AB" % (mab))
    #print("%.2e erg/s/cm^2/Hz" % (fnu))
    #print("%.2e Jy" % (fnu/Jy))
    #print("%.4f uJy" % (fnu/(1e-6*Jy)))
    #print("%.1f nJy" % (fnu/(1e-9*Jy)))
    #print()
    #print("%.2e erg/s/cm^2/Ang" % (flambda))
    #print()
    
    return flambda


# In[ ]:


print(os.environ['GRIZLI'])


# In[ ]:


######################
#HOME_PATH = "/Users/gwalth/data/Roman/grizli/"
GRIZLI_PATH = os.environ['GRIZLI']
sens_file = GRIZLI_PATH + '/CONF/Euclid/CONF11/SENS_A.fits'
sens_tbl = Table.read(sens_file)
print(sens_tbl.colnames)

from scipy import interpolate

R = interpolate.interp1d(sens_tbl['WAVELENGTH'],sens_tbl['SENSITIVITY'])
w0,w1 = sens_tbl['WAVELENGTH'][0],sens_tbl['WAVELENGTH'][-1]

w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], bin=0)
#w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], 
#                                          ivar=1./(beams[key].grism.data['ERR'])**2,bin=0)

i0 = np.argmin(np.abs(w-w0))
i1 = np.argmin(np.abs(w-w1))

w = w[i0:i1]
f = f[i0:i1]
e = e[i0:i1]
    
response = [R(w_el) for w_el in w]
    
flux2 = f/response
noise2 = e/response

fig = plt.figure()
p1 = fig.add_subplot(111)
# plot spectrum
p1.plot(w/10., flux2, c="orange", label="extracted spectra") # linestyle='steps-mid')
# factor of 10 off, could be WAV_AB being wrong (i.e. 1600 instead of 160)

p1.set_xlim(1190,1910.0)
p1.set_ylim(1e-15,1e-13)
p1.set_yscale("log")
p1.set_xlabel("Wavelength [nm]")
#p1.set_ylabel("e-/s per erg/s/cm$^2$/Ang") # according to axe_manual

p1.legend(loc=1)
plt.show()


# In[ ]:


from scipy import integrate

### Plot 1D spectra

fig = plt.figure(figsize=(6,4))
p1 = fig.add_subplot(111)
key = "NISP-GLWv1"
w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], bin=0)
#w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], 
#                                          ivar=1./(beams[key].grism.data['ERR'])**2,bin=0)

#sf = f/np.nanmax(f)
sf = f
std = np.nanstd(sf)
med = np.nanmedian(sf)

sig = 2.
y0 = med-sig*std
y1 = med+sig*std

p1.text(0.05,0.9,"ID = %i" % (id),transform=p1.transAxes)
p1.plot(w/1.e4, sf, c="k") # linestyle='steps-mid')

#p1.plot(w/1.e4, e, c="r") # linestyle='steps-mid')
p1.set_ylabel("Flux [Arbitrary]")
p1.set_xlabel("Wavelength [um]")

p1.set_ylim(y0,y1)


# In[ ]:





# In[ ]:





# ## Loop over all objects and fit their redshifts
# [top](#Table-of-Contents)

# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Extractions'))


# In[ ]:


total_Nl = 0

all_Nl = []
for i,euclid_srcs in enumerate(all_euclid_srcs):
    euclid_all,euclid_magcut,euclid_extract = euclid_srcs
    Nl = [id for id in euclid_extract['NUMBER']]
    print(all_det[i],"number_of_sources =",len(Nl))
    all_Nl.append(Nl)
    
    total_Nl += len(Nl)

print()
print("Total sources for all detectors =",  total_Nl)
print(all_Nl)


# In[ ]:


import importlib
importlib.reload(grizli.fitting)


# In[ ]:


import importlib
importlib.reload(grizli.multifit)


# In[ ]:


# subset
all_Nl = [all_Nl[0]]


# In[ ]:


T0 = time.time()

os.chdir(os.path.join(HOME_PATH, root, 'Extractions'))

#fwhm = 395 # km/s
fwhm = 400 # km/s

#t0 = utils.load_templates(fwhm=fwhm, line_complexes=True, fsps_templates=True) # redshift fits, fixed line ratios
#t1 = utils.load_templates(fwhm=fwhm, line_complexes=False, fsps_templates=True) # final fits

# Fitting templates

# First is set with combined emission line complexes for the redshift fit 
# (don't allow infinite freedom) of the line ratios / fluxes
temp0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                     full_line_list=None,  continuum_list=None, 
                                     fsps_templates=True)

# Second set has individual line templates for fitting the line fluxes
temp1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                     full_line_list=None, continuum_list=None, 
                                     fsps_templates=True)

#print(temp0)
#print(len(temp0))
#print(temp1)
#print(len(temp1))

for i,Nl in enumerate(all_Nl):
    
    print("Fitting redshifts for %s" % (all_det[i]))
    euclid = all_euclid[i]
    
    group_name = root + "_" + all_det[i]
    
    t0 = time.time()
    
    if not os.path.exists(all_det[i]):
        os.mkdir(all_det[i])    
    os.chdir(all_det[i])

    for j,id in enumerate(Nl):

        print("id =",id)
        print("%i of %i" % (j+1,len(Nl)))
    
        #beams = OrderedDict()

        is_cgs, spectrum_1d, b = euclid.object_dispersers[id]
        cutout = grizli.model.BeamCutout(euclid, b['A'], min_sens=0,) # min_mask=0) 

        cutout.beam.compute_model()  
        cutout.contam = cutout.beam.cutout_from_full_image(euclid.model)
        if id in euclid.object_dispersers:
            cutout.contam -= cutout.beam.model

        hdu = cutout.write_fits(get_hdu=True)
        ext = 0
        hdu[ext].header['EXPTIME'] = hdu['SCI'].header['EXPTIME']

        beam = 'beam_%05d.grism.A.fits' % (id)
        hdu.writeto(beam,overwrite=True)

        mb = multifit.MultiBeam([beam], fcontam=0.2, group_name=group_name, psf=False, min_sens=0.05)
        mb.write_master_fits()
    
        # kludge
        os.remove(beam)
        ###########################
        
        fitting.run_all(id, temp0, temp1, fit_only_beams=True, fwhm=fwhm, zr=[0.05, 3.0], 
                        dz=[0.004, 0.0002], fitter=['nnls', 'bounded'], group_name=group_name)
    

    os.chdir("..")
    
    t1 = time.time()
    print()
    print("Finished %s in %.1f seconds" % (all_det[i],t1-t0))
    print()
    
    
T1 = time.time()

print()
print("Finished in %.1f seconds" % (T1-T0))


# In[ ]:


print(os.getcwd())


# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Extractions'))


# In[ ]:


#id = 29
#det = 11

full = "DET%i/%s_DET%i_%05d.full.fits" % (det,root,det,id)
row = "DET%i/%s_DET%i_%05d.row.fits" % (det,root,det,id)
oned = "DET%i/%s_DET%i_%05d.1D.fits" % (det,root,det,id)


hdu = pyfits.open(full)
hdu.info()
head = hdu[0].header
ra = head["RA"]
dec = head["DEC"]
print(ra,dec)

hdu = pyfits.open(row)
hdu.info()
#hdu[0].header

hdu = pyfits.open(oned)
hdu.info()
hdu[0].header


# In[ ]:


prefix = '{0}_DET{1}_{2:05d}'.format(root,det,id)
#display_grizli(prefix, path="DET%i" % (det), dispersers=["RED"], 
#               w0=1.18, w1=1.92)
display_grizli(prefix, path="DET%i" % (det), dispersers=["RED"], 
               w0=1.18, w1=1.92, ) #norm=1e-19)


# In[ ]:


#print(primer)
print(primer.colnames)
#print([col for col in primer.colnames if "TU_" in col])
#bands = [col for col in primer.colnames if "_MAG" in col]
bands = [col for col in primer.colnames if "TU_" in col]


Euclid_bands = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG']
print(ra,dec)

dr = np.sqrt((primer['RA']-ra)**2*np.cos(dec/180*np.pi)**2 + 
             (primer['DEC']-dec)**2)*3600.

ind = np.argmin(dr)
source_id = primer['SOURCE_ID'][ind]
obj_mag = primer['TU_FNU_H_NISP_MAG'][ind]
print('SOURCE_ID:%d, H_mag=%.2f, dr=%.2f"' % (source_id, obj_mag, np.min(dr)))
print(primer['Z_OBS'][ind])
fnu = primer['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG'][ind]

fnu = np.array(list(fnu))
#print(dir(fnu))
#print(primer[bands][ind])
print(fnu)

Jy = 1e-23 # erg/s/cm^2/Hz

m = -2.5*np.log10(fnu*Jy)-48.60
print(m)

primer_small = primer['SOURCE_ID','RA','DEC','TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 
                      'TU_FNU_H_NISP_MAG']

tbl = Table(rows=primer_small[ind])
#print(tbl.colnames)
print(tbl)
print()

# rename columns
tbl.rename_column("SOURCE_ID","id")   
tbl.rename_column("RA","ra")   
tbl.rename_column("DEC","dec")   

for col in tbl.colnames:
    if "TU_FNU_" in col:
        print(col)
        new_col = col.replace("TU_FNU_","").replace("_MAG","_FLUX")
        print(new_col)
        tbl.rename_column(col,new_col)    
print()
    
# adjust flux and add errors
for col in tbl.colnames:
    if "_FLUX" in col:
        print(col)
        tbl[col] /= 1e-6
        #row[col]
        new_col = col + "ERR"
        tbl[new_col] = 0.5

tbl["z_spec"] = -1.0
        
#print(tbl)

euclid_phot_file = "Euclid_phot.fits"
tbl.write(euclid_phot_file, overwrite=True)            
tbl.show_in_notebook()


# In[ ]:


# Requires eazy-py:  https://github.com/gbrammer/eazy-py
import eazy
print('\n Eazy-py version: ', eazy.__version__)


# In[ ]:


# Preparation for eazy-py
eazy.symlink_eazy_inputs()


# In[ ]:


### Initialize **eazy.photoz** object

params = {}

translate_file = 'Euclid_phot.translate'
params['CATALOG_FILE'] = 'Euclid_phot.fits'
params['MAIN_OUTPUT_FILE'] = 'Euclid_phot.eazypy'
params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'

# Galactic extinction
params['MW_EBV'] = 0.0
params['Z_STEP'] = 0.002
params['Z_MAX'] = 3.
params['FIX_ZSPEC']='n'
#params['PRIOR_FILTER'] = 205

ez = eazy.photoz.PhotoZ(param_file=None, translate_file=translate_file, 
                        zeropoint_file=None, params=params, 
                        load_prior=True, load_products=False)


# In[ ]:


from grizli.pipeline import photoz
## Grism fitting arguments created in Grizli-Pipeline
#args = np.load('fit_args.npy', allow_pickle=True)[0]

## First-pass redshift templates, similar to the eazy templates but 
## with separate emission lines
#t0 = args['t0'] # read earlier

#############
## Make a helper object for generating photometry in a format that grizli 
## understands. 

## Passing the parameters precomputes a function to quickly interpolate
## the templates through the broad-band filters.  It's not required, 
## but makes the fitting much faster.
## 
## `zgrid` defaults to ez.zgrid, be explicit here to show you can 
## change it. 
print(ez.zgrid)

phot_obj = photoz.EazyPhot(ez, grizli_templates=temp0, zgrid=ez.zgrid) 


# In[ ]:


ez.fit_parallel(n_proc=4)
ez.error_residuals()

print('Get physical parameters')
ez.standard_output()


# In[ ]:


# Show SEDs with best-fit templates and p(z)
ez.show_fit(0, id_is_idx=True)


# In[ ]:


### Spline templates for dummy grism continuum fits
wspline = np.arange(4200, 2.5e4)
Rspline = 50
df_spl = len(utils.log_zgrid(zr=[wspline[0], wspline[-1]], dz=1./Rspline))
tspline = utils.bspline_templates(wspline, df=df_spl+2, log=True, clip=0.0001)


# In[ ]:


print(source_id) # primer
print(id)        # sims.catalog

catalog = all_euclid[0].catalog["NUMBER"]
idx = list(catalog).index(id)
print(idx)
#print(list(catalog[indices]))
#print(all_euclid[0].catalog[indices])


# In[ ]:


## This isn't necessary for general fitting, but 
## load the grism spectrum here for demonstrating the grism/photometry scaling
group_name = root + "_DET%i" % (det)
beams_file = "%s_DET%i_%05d.beams.fits" % (root,det,id)
mb = multifit.MultiBeam(beams_file, fcontam=0.2, group_name=group_name)


# In[ ]:


# Generate the `phot` dictionary
phot, ii, dd = phot_obj.get_phot_dict(mb.ra, mb.dec)
label = "Euclid Catalog ID: {0}, dr={1:.2f}, zphot={2:.3f}"
print(label.format(ez.cat['id'][ii], dd, ez.zbest[ii]))

print('\n`phot` keys:', list(phot.keys()))
for k in phot:
    print('\n'+k+':\n', phot[k])
    
# Initialize photometry for the MultiBeam object
mb.set_photometry(**phot)


# In[ ]:


# parametric template fit to get reasonable background
sfit = mb.template_at_z(templates=tspline, fit_background=True, 
                        include_photometry=False)
fig = mb.oned_figure(tfit=sfit)

ax = fig.axes[0]
ax.errorbar(mb.photom_pivot/1.e4, mb.photom_flam/1.e-19, 
            mb.photom_eflam/1.e-19, 
            marker='s', color='k', alpha=0.4, linestyle='None',
            label='Euclid photometry')

ax.legend(loc='upper left', fontsize=8)

ax.set_ylim(0,100)
ax.set_xlim(0.6,1.9)

print(mb.photom_pivot/1.e4)
print(mb.photom_flam/1.e-19)


# In[ ]:


## First example:  no rescaling
z_phot = ez.zbest[0]

# Reset scale parameter
if hasattr(mb,'pscale'):
    delattr(mb, 'pscale')
    
#t1 = args['t1']
tfit = mb.template_at_z(z=z_phot)
print('No rescaling, chi-squared={0:.1f}'.format(tfit['chi2']))
fig = fitting.full_sed_plot(mb, tfit, zfit=None, bin=4)


# In[ ]:


# Reset scale parameter
if hasattr(mb,'pscale'):
    delattr(mb, 'pscale')

# Template rescaling, simple multiplicative factor
scl = mb.scale_to_photometry(order=0)

# has funny units of polynomial coefficients times 10**power, 
# see `grizli.fitting.GroupFitter.compute_scale_array`
# Scale value is the inverse, so, e.g., 
# scl.x = [8.89] means scale the grism spectrum by 10/8.89=1.12
print(scl.x) 

mb.pscale = scl.x 

# Redo template fit
tfit = mb.template_at_z(z=z_phot)
print(z_phot)
print('Simple scaling, chi-squared={0:.1f}'.format(tfit['chi2']))
fig = fitting.full_sed_plot(mb, tfit, zfit=None, bin=4)


# In[ ]:


# Reset scale parameter
if hasattr(mb,'pscale'):
    delattr(mb, 'pscale')
    mb.compute_model()

# Template rescaling, linear fit
scl = mb.scale_to_photometry(order=0)

# has funny units of polynomial coefficients times 10**power, 
# see `grizli.fitting.GroupFitter.compute_scale_array`
# Scale value is the inverse, so, e.g., 
# scl.x = [8.89] means scale the grism spectrum by 10/8.89=1.12
print(scl.x) 

mb.pscale = scl.x 

# Redo template fit
tfit = mb.template_at_z(z=z_phot)
print('Simple scaling, chi-squared={0:.1f}'.format(tfit['chi2']))
fig = fitting.full_sed_plot(mb, tfit, zfit=None, bin=4)


# In[ ]:


if hasattr(mb,'pscale'):
    delattr(mb, 'pscale')
    mb.compute_model()

# Now run the full redshift fit script with the photometry, which will also do the scaling
order=0
#fitting.run_all_parallel(id, phot=phot, verbose=False, 
#                         scale_photometry=order+1, zr=[0.05, 3.0])

fitting.run_all(id, temp0, temp1, phot=phot, scale_photometry=order+1, fit_only_beams=True, fwhm=fwhm, 
                zr=[0.05, 3.0], dz=[0.008, 0.0004], fitter=['nnls', 'bounded'], group_name=group_name)

#fitting.run_all(id, temp0, temp1, fit_only_beams=True, fwhm=fwhm, zr=[0.05, 3.0], 
#                        dz=[0.004, 0.0002], fitter=['nnls', 'bounded'], group_name=group_name)
    


# In[ ]:


zfit = pyfits.open('{0}_DET{1}_{2:05d}.full.fits'.format(root, det, id))
z_grism = zfit['ZFIT_STACK'].header['Z_MAP']
print('Best redshift: {0:.4f}'.format(z_grism))

# Compare PDFs
pztab = utils.GTable.gread(zfit['ZFIT_STACK'])
plt.plot(pztab['zgrid'], pztab['pdf'], label='grism+Euclid')

plt.plot(ez.zgrid, np.exp(ez.lnp[0,:]), label='photo-z')

plt.semilogy()
plt.xlim(z_grism-0.05, z_grism+0.05); plt.ylim(1.e-10, 1000)
plt.xlabel(r'$z$'); plt.ylabel(r'$p(z)$')
plt.grid()
plt.legend()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ## Extract a single 2D spectrum 
# Science, model and contamination based on RA and Dec
# 
# [top](#Table-of-Contents)

# In[ ]:


det_ind = 0 # DET11
Euclid = all_euclid[det_ind]
phot = all_phot[det_ind]

## Find the object ID near coordinates (r0,d0), 
## SExtractor IDs might not be constant across platforms

#r0,d0 = 228.6882041, 6.2969624 # cont src
#r0,d0 = 228.6882147, 6.3239512 #  emission line src
#r0,d0 = 228.6881854, 6.3209442
#r0,d0 = 228.6882128, 6.3028898 
#r0,d0 = 228.6881973, 6.3179116 
r0,d0 = 228.6881976, 6.3149738 


dr = np.sqrt((phot['X_WORLD']-r0)**2*np.cos(d0/180*np.pi)**2 + 
             (phot['Y_WORLD']-d0)**2)*3600.
id = phot['NUMBER'][np.argmin(dr)]
obj_mag = phot['MAG_AUTO'][np.argmin(dr)]
print('ID:%d, mag=%.2f, dr=%.2f"' %(id, obj_mag, np.min(dr)))

beams = OrderedDict()

ix = Euclid.catalog['id'] == id
x0, y0 = Euclid.catalog['x_flt'][ix][0], Euclid.catalog['y_flt'][ix][0]
print(Euclid.direct.instrument, x0, y0)
#print(Roman.wcs.pscale)
#dim = 18*0.135/sim.flt_wcs.pscale 
#beam = grizli.model.BeamCutout(id=id, x=x0, y=y0, 
#                               cutout_dimensions=np.cast[int]((dim, dim)), 
#                               conf=sim.conf, GrismFLT=sim)

print(Euclid.object_dispersers[id])

is_cgs, spectrum_1d, b = Euclid.object_dispersers[id]
#print(b)

cutout = grizli.model.BeamCutout(Euclid, b['A'], min_sens=0,) # min_mask=0) 

cutout.beam.compute_model()  
cutout.contam = cutout.beam.cutout_from_full_image(Euclid.model)
if id in Euclid.object_dispersers:
    cutout.contam -= cutout.beam.model

print(dir(cutout.beam))
#print(cutout.contam)
#print(cutout.beam.model)
#print(dir(cutout.beam))
print(cutout.beam.seg)
print(cutout.beam.seg.shape)

img = cutout.grism.data['SCI']*1
X = img.flatten()
std = np.std(X)
med = np.median(X)

sig = 1.
vmin = med-sig*std
vmax = med+sig*std

fig = plt.figure()
p1 = fig.add_subplot(111)
p1.imshow(cutout.grism.data['SCI']*1, origin='lower',cmap='gray_r',vmin=vmin,vmax=vmax)
p1.set_xlabel("X [pixels]")
p1.set_ylabel("Y [pixels]")

#p1.imshow(cutout.grism.data['ERR']*1)#,vmin=-0.1,vmax=0.1, origin='lower',cmap='gray_r')
    
beams[Euclid.grism.instrument] = cutout

cutout.write_fits() # still learning about the output


# ## 1D Spectral Extraction
# [top](#Table-of-Contents)

# In[ ]:


from scipy import integrate

### Plot 1D spectra
fig = plt.figure(figsize=(10,6))
p1 = fig.add_subplot(111)
key = "NISP-GLWv1"
w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], bin=0)
#w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], 
#                                          ivar=1./(beams[key].grism.data['ERR'])**2,bin=0)


#sf = f/np.nanmax(f)
sf = f
std = np.nanstd(sf)
med = np.nanmedian(sf)

sig = 2.
y0 = med-sig*std
y1 = med+sig*std

p1.text(0.05,0.9,"ID = %i" % (id),transform=p1.transAxes)
p1.plot(w/1.e4, sf, c="k") # linestyle='steps-mid')

#p1.plot(w/1.e4, e, c="r") # linestyle='steps-mid')
p1.set_ylabel("Flux [Arbitrary]")
p1.set_xlabel("Wavelength [um]")

p1.set_ylim(y0,y1)


# ## Display Redshift Fit
# [top](#Table-of-Contents)

# In[ ]:


det_ind = 0
os.chdir(os.path.join(HOME_PATH, root, 'Extractions',all_det[det_ind]))
group_name = root + "_" + all_det[det_ind]


# In[ ]:


display_grizli(group_name, id, w0=1.15, w1=1.95, labels=1)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


#########################
# copy sextractor files #
#########################

# cd ~/data/Roman/grizli/my_roman_sims/Prep

# rsync -avz Roman.param cygnusc.ipac.caltech.edu:/home/gwalth/data/Roman/grizli/sims/sim_v2/Prep
# rsync -avz Roman.sex   cygnusc.ipac.caltech.edu:/home/gwalth/data/Roman/grizli/sims/sim_v2/Prep
# rsync -avz default.nnw cygnusc.ipac.caltech.edu:/home/gwalth/data/Roman/grizli/sims/sim_v2/Prep


# In[ ]:


##############################
# compare branch differences #
##############################

# git diff roman_sims_v1_gwalth..master -- *.py

# cd ~/python/grizli_1.3.2/grizli

# rsync -avz grismconf.py cygnusc.ipac.caltech.edu:/home/gwalth/python/src/grizli/grizli
# rsync -avz model.py     cygnusc.ipac.caltech.edu:/home/gwalth/python/src/grizli/grizli
# rsync -avz multifit.py  cygnusc.ipac.caltech.edu:/home/gwalth/python/src/grizli/grizli
# rsync -avz utils.py     cygnusc.ipac.caltech.edu:/home/gwalth/python/src/grizli/grizli


# In[ ]:


#########################
# copy grism conf files #
#########################

# cd ~/data/Roman/grizli/grizli/CONF

# rsync -avz Roman.G150-v?-GLW.conf cygnusc.ipac.caltech.edu:/home/gwalth/data/Roman/grizli/grizli/CONF
# rsync -avz sens_0720_2020.fits    cygnusc.ipac.caltech.edu:/home/gwalth/data/Roman/grizli/grizli/CONF


# In[ ]:





# In[ ]:





# In[ ]:





# # Appendix - Old

# ## aXeSIM predictions based on conf file
# [top](#Table-of-Contents)

# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))


# In[ ]:


from grizli import grismconf


# In[ ]:


# This one gives the error
import os
os.chdir("/Users/gwalth/data/Roman/grizli/grizli/CONF")
conf = grismconf.load_grism_config("./Euclid.Gred.0.conf")

# This one works!
#conf = grismconf.load_grism_config("../../../grizli/CONF/Roman.det1.07242020.conf")
#conf = grismconf.load_grism_config("../../../grizli/CONF/Roman_bad/Roman.det1.07242020.conf")

dx = conf.dxlam["A"]
print(dx)
x0,y0 = 1024,1024

dy,lam = conf.get_beam_trace(x=x0,y=y0,dx=dx)
print(dy)
print(lam)


# In[ ]:


conf.show_beams()


# In[ ]:


# This one gives the error
import os
os.chdir("/Users/gwalth/data/Roman/grizli/grizli/CONF/CONF11")
conf = grismconf.load_grism_config("./NISP_RGS000_21.conf")
#conf = grismconf.load_grism_config("./NISP_RGS000_21.conf.beamA")

# This one works!
#conf = grismconf.load_grism_config("../../../grizli/CONF/Roman.det1.07242020.conf")
#conf = grismconf.load_grism_config("../../../grizli/CONF/Roman_bad/Roman.det1.07242020.conf")

dx = conf.dxlam["A"]
print(dx)
x0,y0 = 1024,1024

dy,lam = conf.get_beam_trace(x=x0,y=y0,dx=dx)
print(dy)
print(lam)


# In[ ]:


conf.show_beams()


# In[ ]:


print(os.getcwd())


# ## Show 2D beam
# [top](#Table-of-Contents)

# In[ ]:


## Spectrum with lines & noise
#spectrum_file = os.path.join(os.path.dirname(grizli.__file__), 'data/erb2010.dat')
#spectrum_file = '../../../grizli/templates/erb2010.dat'
#erb = np.loadtxt(spectrum_file, unpack=True)
#z = 2.0 # test redshift

## normalize spectrum to unity to use normalization defined in the direct image
#import pysynphot as S
#spec = S.ArraySpectrum(erb[0], erb[1], fluxunits='flam')
#spec = spec.redshift(z).renorm(1., 'flam', S.ObsBandpass('wfc3,ir,f140w'))
#spec.convert('flam') # bug in pysynphot, now units being converted automatically above? (11/10/16)

fig = plt.figure(figsize=(14,4))
for i, key in enumerate(beams.keys()):
    #     beams[key].compute_model(beams[key].thumb, id=beams[key].id, 
    #                              xspec=spec.wave, yspec=spec.flux)
    #print(key)
    #beams[key].beam.compute_model(spectrum_1d=[spec.wave, spec.flux]) 
    #beams[key].beam.compute_model()
    
    gdata = beams[key].grism.data['SCI']*1
    gmodel = beams[key].model
    gcontam = beams[key].contam
    
    print(gdata)
    print(gmodel)
    print(gcontam)
    
    print(np.min(gdata),np.max(gdata))
    print(np.min(gmodel),np.max(gmodel))
    print(np.min(gcontam),np.max(gcontam))
    
    print(gdata.shape)
    print(gmodel.shape)
    print(gcontam.shape)
    
    axl = fig.add_subplot(321+i*2)
    #axl.imshow(beams[key].model + beams[key].grism.data['SCI']*1, interpolation='Nearest', 
    #       origin='lower', vmin=-0.006, vmax=0.16, cmap='gray_r', aspect='auto')
    
    axl.imshow(beams[key].grism.data['SCI']*1, interpolation='Nearest', 
           origin='lower', vmin=-0.1,vmax=0.1, cmap='gray_r', aspect='auto')
    
    
    axr = fig.add_subplot(321+i*2+1)
    #axr.imshow(beams[key].contam + beams[key].grism.data['SCI'] + beams[key].model, 
    #           interpolation='Nearest', 
    #           origin='lower', vmin=-1, vmax=2, cmap='gray_r', aspect='auto')

    axr.imshow(beams[key].contam , 
               interpolation='Nearest', 
               origin='lower', vmin=-0.1,vmax=0.1, cmap='gray_r', aspect='auto')
    
    #axl.set_title('%s - %s' %(key, beams[key].grism.filter))
    #for ax in [axl, axr]:
    #    beams[key].beam.twod_axis_labels(wscale=1.e4, mpl_axis=ax)
    #    beams[key].beam.twod_xlim(1.3,1.75, wscale=1.e4, mpl_axis=ax)
    #    if i < 2:
    #        ax.set_xticklabels([])
    #    else:
    #        ax.set_xlabel(r'$\lambda$')

fig.tight_layout(pad=0.5)


# In[ ]:


def ab2flux(mab,eff_wav):
    c = 2.9979E10 # cm/s
    Ang = 1E-8    # cm
    Jy = 1E-23    # erg/s/cm^2/Hz
    
    # mab = -2.5*np.log10(fnu) - 48.6 
    fnu = 10**(-0.4*(mab + 48.6))          # erg/s/cm^2/Hz
    flambda = fnu*(c/(eff_wav*Ang)**2)*Ang # erg/s/cm^2/Ang
   
    #print("%.2f AB" % (mab))
    #print("%.2e erg/s/cm^2/Hz" % (fnu))
    #print("%.2e Jy" % (fnu/Jy))
    #print("%.4f uJy" % (fnu/(1e-6*Jy)))
    #print("%.1f nJy" % (fnu/(1e-9*Jy)))
    #print()
    #print("%.2e erg/s/cm^2/Ang" % (flambda))
    #print()
    
    return flambda


# In[ ]:


# determine flux in Ang for the source in the direct image

mab = phot[phot['NUMBER']==id]['MAG_AUTO'][0]

f_scale = ab2flux(mab,15800.)


# ## Simple SN calculations based on the spcontetc
# [top](#Table-of-Contents)

# In[ ]:


# Per pix S/N=1
mag_per_pix_w0 = 21.9       # 1.00 um
mag_per_pix_w1 = 21.9       # 1.93 um
mag_per_pix_deepest = 22.7  # 1.34 um


# In[ ]:


ab2flux(mag_per_pix_w0,10000.)
ab2flux(mag_per_pix_deepest,13400.)


# In[ ]:


print(os.getcwd())
spcontetc_1n = Table.read("../etc/spcontetc_1n_0.2.dat",format="ascii")
spcontetc_8n = Table.read("../etc/spcontetc_8n_0.2.dat",format="ascii")
#spcontetc_1n = Table.read("../etc/spcontetc_1n_0.5.dat",format="ascii")
#spcontetc_8n = Table.read("../etc/spcontetc_8n_0.5.dat",format="ascii")
print(spcontetc_1n)

#|  um  | arcsec |      |exp per Jy|per pix |per resl|1e4 km/s|
fig = plt.figure()
p = fig.add_subplot(111)
#for n in range(5,8):
for n in range(5,7):
    #p.plot(spcontetc_1n['col1'],spcontetc_1n['col%i' % (n)]) # AB
    #p.plot(spcontetc_8n['col1'],spcontetc_8n['col%i' % (n)]) # AB
    
    f1 = ab2flux(spcontetc_1n['col%i' % (n)],spcontetc_1n['col1']*1e4)
    f8 = ab2flux(spcontetc_8n['col%i' % (n)],spcontetc_8n['col1']*1e4)
    
    p.plot(spcontetc_1n['col1'],f1)  
    #p.plot(spcontetc_8n['col1'],f8) 

    

#p.set_ylim(0e-18,2e-18)
p.set_xlabel(r'$\lambda$ [micron]')
p.set_ylabel(r'Flux [erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$]')
#p.set_ylabel(r'Mag [AB]')

#p.invert_yaxis()
    
plt.show()


# Roughly estimated from spcontetc (>1.4 micron):
# 
# <span style="color:red">5 sigma limit = 1.75e-18 erg/s/cm^2/Ang</span>
# 
# <span style="color:red">1 sigma limit = 1.00e-18 erg/s/cm^2/Ang</span>
# 

# ## Simple SN calculations based on the pzcaletc
# [top](#Table-of-Contents)

# In[ ]:


print(os.getcwd())
pzcaletc_1n = Table.read("../etc/pzcaletc_1n.dat",format="ascii")
pzcaletc_8n = Table.read("../etc/pzcaletc_8n.dat",format="ascii")
print(pzcaletc_1n)

fig = plt.figure()
p = fig.add_subplot(111)
#for n in range(2,13):
for n in range(2,3):
    #p.plot(pzcaletc_1n['col1'],pzcaletc_1n['col%i' % (n)]) # W/m^2
    p.plot(pzcaletc_1n['col1'],pzcaletc_1n['col%i' % (n)]*1000.) # erg/s/cm^2
    p.plot(pzcaletc_8n['col1'],pzcaletc_8n['col%i' % (n)]*1000.) # erg/s/cm^2

    
#p.legend(fontsize=14)
#p.set_xlim(w1-0.05, w2+0.05)
#p.set_ylim(-1e-21,3e-20)
p.set_xlabel(r'$\lambda$ [micron]')
p.set_ylabel(r'Flux [erg s$^{-1}$ cm$^{-2}$]')
    
plt.show()


# ## Simple SN calculations based on the apttables2021
# [top](#Table-of-Contents)

# In[ ]:


print(os.getcwd())
apt_rad0_0_zod1_2 = Table.read("../apttables2021/grism/spec_half-light-rad0.0_zod1.2.dat",format="ascii",
                               names=('AB','1.05','1.20','1.40','1.60','1.80','2.00'))
apt_rad0_2_zod1_2 = Table.read("../apttables2021/grism/spec_half-light-rad0.2_zod1.2.dat",format="ascii",
                               names=('AB','1.05','1.20','1.40','1.60','1.80','2.00'))
apt_rad0_3_zod1_2 = Table.read("../apttables2021/grism/spec_half-light-rad0.3_zod1.2.dat",format="ascii",
                               names=('AB','1.05','1.20','1.40','1.60','1.80','2.00'))
#print(apt_rad0_0_zod1_2)
#print(apt_rad0_2_zod1_2)
#print(apt_rad0_3_zod1_2)

fig = plt.figure()
p = fig.add_subplot(111)

wlist = ['1.05','1.20','1.40','1.60','1.80','2.00']

#for w in wlist:
#    p.plot(apt_rad0_0_zod1_2[w],apt_rad0_0_zod1_2['AB'],label=w) 

p.scatter(apt_rad0_0_zod1_2['1.20'],apt_rad0_0_zod1_2['AB'],label='half-light radius 0.0')
p.scatter(apt_rad0_2_zod1_2['1.20'],apt_rad0_0_zod1_2['AB'],label='half-light radius 0.2')
p.scatter(apt_rad0_3_zod1_2['1.20'],apt_rad0_0_zod1_2['AB'],label='half-light radius 0.3')
    
p.legend(fontsize=14)
p.set_xscale("log")

#p.set_xlim(0,2000)
#p.set_ylim(20,22)
    
plt.show()


# In[ ]:


##############
# My crude ETC
# based on 
# https://roman.gsfc.nasa.gov/science/apttables2021/table-grism.html
##############

# notes:
# might need a fit cutoff to t = 1e7 

from scipy.optimize import leastsq

t = 301.

param = [1,0]
err = 1.0

# y = ax + b
line_fn = lambda p,x: p[1]*x+p[0]

err_fn = lambda p: (line_fn(p,x) - y)/err

wlist = ['1.05','1.20','1.40','1.60','1.80','2.00'] # microns

mag_ab_5sig = []

fig = plt.figure()
p = fig.add_subplot(111)

#x = np.log10(apt_rad0_0_zod1_2['1.20'])
y = apt_rad0_2_zod1_2['AB']
for w in wlist:
    print(w)
    x = np.log10(apt_rad0_2_zod1_2[w])

    sol = leastsq(err_fn,param,full_output=1)
    print(sol[0])

    mag_ab = line_fn(sol[0],np.log10(t))
    print(mag_ab)
    mag_ab_5sig.append(mag_ab)
    print()

    logt = np.arange(2.0,8.6,0.1)
    mag = line_fn(sol[0],logt)
    

    p.plot(10**logt,mag,c="k")
    p.scatter(10**x,y,label=w,s=50,alpha=0.7)
    
p.legend(fontsize=14)
p.set_xscale("log")

p.set_xlim(400,2.5e8)
p.set_ylim(20.8,27.2)

p.invert_yaxis()
plt.show()


# In[ ]:


fig = plt.figure(figsize=(10,4))
p1 = fig.add_subplot(121)

warr = 1e4*np.array([float(w0) for w0 in wlist]) # Angstroms


mag_ab_5sig = np.array(mag_ab_5sig)

p1.plot(warr,mag_ab_5sig)
p1.set_xlabel(r'$\lambda$ [Ang]')
p1.set_ylabel(r'mag [AB]')
p1.invert_yaxis()


#p.set_xlim(0,2000)
#p.set_ylim(20,22)
p2 = fig.add_subplot(122)

#flux_5sig = [ab2flux(mag_ab_5sig[i],w0) for i,w0 in enumerate(warr)]
flux_5sig = ab2flux(mag_ab_5sig,warr) 
p2.plot(warr,flux_5sig)

p2.set_xlabel(r'$\lambda$ [Ang]')
p2.set_ylabel(r'Flux [erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$]')
    
plt.show()


# Comparison between spcontetc and apttables2021

# In[ ]:


print(os.getcwd())
spcontetc_1n = Table.read("../etc/spcontetc_1n_0.2.dat",format="ascii")
spcontetc_8n = Table.read("../etc/spcontetc_8n_0.2.dat",format="ascii")
#spcontetc_1n = Table.read("../etc/spcontetc_1n_0.5.dat",format="ascii")
#spcontetc_8n = Table.read("../etc/spcontetc_8n_0.5.dat",format="ascii")
print(spcontetc_1n)

fig = plt.figure()
p = fig.add_subplot(111)
#for n in range(5,8):
for n in range(5,7):
    #p.plot(spcontetc_1n['col1'],spcontetc_1n['col%i' % (n)]) # AB
    #p.plot(spcontetc_8n['col1'],spcontetc_8n['col%i' % (n)]) # AB
    
    f1 = ab2flux(spcontetc_1n['col%i' % (n)],spcontetc_1n['col1']*1e4)
    f8 = ab2flux(spcontetc_8n['col%i' % (n)],spcontetc_8n['col1']*1e4)
    
    p.plot(spcontetc_1n['col1'],f1,label='spcontetc 0.2"')  
    #p.plot(spcontetc_8n['col1'],f8) 


p.plot(warr/1e4,flux_5sig,label="apttables2021")

#p.set_ylim(0e-18,2e-18)
p.set_xlabel(r'$\lambda$ [micron]')
p.set_ylabel(r'Flux [erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$]')
#p.set_ylabel(r'Mag [AB]')

#p.invert_yaxis()
p.legend()
    
plt.show()


# In[ ]:


def poly_n(p, x):
    p = np.array(p)
    y = np.zeros((x.shape[0],))
    for n in np.arange(p.shape[0]):
        y += p[n]*x**n
    return y

def solve_for_y(poly_coeffs, y):
    pc = poly_coeffs.copy()
    print(pc)
    pc[0] -= y
    print(pc)
    return np.roots(pc)


print(os.getcwd())
apt_SN_55sec_zod1_2 = Table.read("../apttables2021/imaging/ap_phot_2pix_0.22rad_SN_55sec_zod1.2x.dat",format="ascii",
                               names=('AB','F062','F087','F106','F129','F158','F184','F146','F213'))
apt_SN_75sec_zod1_2 = Table.read("../apttables2021/imaging/ap_phot_2pix_0.22rad_SN_75sec_zod1.2x.dat",format="ascii",
                               names=('AB','F062','F087','F106','F129','F158','F184','F146','F213'))
apt_SN_100sec_zod1_2 = Table.read("../apttables2021/imaging/ap_phot_2pix_0.22rad_SN_100sec_zod1.2x.dat",format="ascii",
                               names=('AB','F062','F087','F106','F129','F158','F184','F146','F213'))


img_table = Table()
img_table['AB'] = apt_SN_55sec_zod1_2['AB']
img_table['55'] = apt_SN_55sec_zod1_2['F158']
img_table['75'] = apt_SN_75sec_zod1_2['F158']
img_table['100'] = apt_SN_100sec_zod1_2['F158']


#print(img_table)

t = 141

# SNR ~ t^(0.5)
# time
x = np.sqrt(np.array([55.,75.,100.]))
param = [1, 0]
err = 1.0

line_fn = lambda p,x: p[1]*x+p[0]

err_fn = lambda p: (line_fn(p,x) - y)/err

fig = plt.figure()
p = fig.add_subplot(111)

SN = [] # at t

###################
# fit for SNR and t
###################
for i,y in enumerate(img_table['55','75','100']):

    # SNR
    y = np.array(list(y))
    
    sol = leastsq(err_fn,param,full_output=1)
    print(sol[0])
    
    t0 = np.arange(50,200,0.1)
    snr = line_fn(sol[0],np.sqrt(t0))
    
    mag = np.ones(snr.shape)*img_table['AB'][i]

    p.plot(mag,snr,c="k")
    
    SN.append(line_fn(sol[0],np.sqrt(t)))

img_table['141'] = SN



param = [1,1,1,1]
# mag
x = img_table['AB']
#print(x)
err = np.ones(x.shape)
#err[:5] = 0.01
#err[-5:] = 0.01

err_fn2 = lambda p: (poly_n(p,x) - y)/err
#####################
# fit for SNR and mag
#####################

dict = {}

for i,t0 in enumerate(['55','75','100','141']):
    
    # SNR
    y = 2.5*np.log10(img_table[t0])
    
    sol = leastsq(err_fn2,param,full_output=1)
    print(sol[0])
    
    mag = np.arange(16.5,26.5,0.1)
    snr = 10**(0.4*poly_n(sol[0],mag))
    
    p.plot(mag,snr,c="k")
    print(len(sol[0]))
    
    dict['sol%s' % t0] = sol[0]

#################
# solving for s/n
#################
print(dict)
sol0 = dict['sol141']
print(solve_for_y(sol0,2.5*np.log10(10.0)))


p.scatter(img_table['AB'],img_table['55'],label='55 sec')
p.scatter(img_table['AB'],img_table['75'],label='75 sec')
p.scatter(img_table['AB'],img_table['100'],label='100 sec')
p.scatter(img_table['AB'],SN,label='141 sec')


p.plot([17,27],[10,10],"-",c="k",alpha=0.5)
p.plot([17,27],[5,5],"--",c="k",alpha=0.5)
p.plot([17,27],[3,3],"-.",c="k",alpha=0.5)

p.legend(fontsize=14)
p.set_yscale("log")

#p.set_xlim(0,2000)
#p.set_ylim(20,22)
#p.invert_yaxis()
    
plt.show()


# In[ ]:


print(os.getcwd())
apt_SN_55sec_zod1_2 = Table.read("../apttables2021/imaging/ap_phot_2pix_0.22rad_SN_55sec_zod1.2x.dat",format="ascii",
                               names=('AB','F062','F087','F106','F129','F158','F184','F146','F213'))
apt_SN_75sec_zod1_2 = Table.read("../apttables2021/imaging/ap_phot_2pix_0.22rad_SN_75sec_zod1.2x.dat",format="ascii",
                               names=('AB','F062','F087','F106','F129','F158','F184','F146','F213'))
apt_SN_100sec_zod1_2 = Table.read("../apttables2021/imaging/ap_phot_2pix_0.22rad_SN_100sec_zod1.2x.dat",format="ascii",
                               names=('AB','F062','F087','F106','F129','F158','F184','F146','F213'))
#print(apt_rad0_0_zod1_2)
#print(apt_rad0_2_zod1_2)
#print(apt_rad0_3_zod1_2)

fig = plt.figure()
p = fig.add_subplot(111)

flist = ['F062','F087','F106','F129','F158','F184','F146','F213']

#for f in flist:
#    p.scatter(apt_SN_55sec_zod1_2['AB'],apt_SN_55sec_zod1_2[f],label=f) 
    #p.scatter(apt_SN_75sec_zod1_2['AB'],apt_SN_75sec_zod1_2[f],label=f)
    #p.scatter(apt_SN_100sec_zod1_2['AB'],apt_SN_100sec_zod1_2[f],label=f) 

p.scatter(apt_SN_55sec_zod1_2['AB'],apt_SN_55sec_zod1_2['F158'],label='55 sec')
p.scatter(apt_SN_75sec_zod1_2['AB'],apt_SN_75sec_zod1_2['F158'],label='75 sec')
p.scatter(apt_SN_100sec_zod1_2['AB'],apt_SN_100sec_zod1_2['F158'],label='100 sec')


p.plot([17,27],[10,10],"-",c="k",alpha=0.5)
p.plot([17,27],[5,5],"--",c="k",alpha=0.5)
p.plot([17,27],[3,3],"-.",c="k",alpha=0.5)

p.legend(fontsize=14)
p.set_yscale("log")

#p.set_xlim(0,2000)
#p.set_ylim(20,22)
#p.invert_yaxis()
    
plt.show()


# ## Roman and Euclid Sensitivity Function
# 
# [top](#Table-of-Contents)

# In[ ]:


import os
print(os.getcwd())
print(HOME_PATH)
#os.chdir(os.path.join(HOME_PATH, root, 'Extraction'))


# ### Gabe's Roman sensitivity function

# In[ ]:


sens_file = HOME_PATH + '/../grizli/CONF/Roman.G150.v1.6.sens.fits' # Gabe's
sens_cat1 = Table.read(sens_file)
print(sens_cat1.colnames)


# ### Anahita's Roman sensitivity function

# In[ ]:


sens_file = HOME_PATH + '/../grizli/CONF/sens_0720_2020.fits'
sens_cat2 = Table.read(sens_file)
print(sens_cat2.colnames)


# ### AstroDeep Euclid sensitivity function

# In[ ]:


sens_file = HOME_PATH + '/../grizli/CONF/Euclid.Gred.1st.sens.0.fits'
sens_cat3 = Table.read(sens_file)
print(sens_cat3.colnames)


# ### TIPS Euclid sensitivity function

# In[ ]:


sens_file = HOME_PATH + '/../grizli/CONF/CONF11/SENS_A.fits'
sens_cat4 = Table.read(sens_file)
print(sens_cat4.colnames)


# In[ ]:


fig = plt.figure()

p1 = fig.add_subplot(111)
#p1.errorbar(sens_cat1['WAVELENGTH'],sens_cat1['SENSITIVITY'],yerr=sens_cat1['ERROR'])
#p1.errorbar(sens_cat2['WAVELENGTH'],sens_cat2['SENSITIVITY'],yerr=sens_cat2['ERROR'])
p1.plot(sens_cat1['WAVELENGTH'],sens_cat1['SENSITIVITY'],label="Roman (Pandeia)")
p1.plot(sens_cat2['WAVELENGTH'],sens_cat2['SENSITIVITY'],label="Roman (Anahita's)")
p1.plot(sens_cat3['WAVELENGTH'],sens_cat3['SENSITIVITY'],label="Euclid (AstroDeep)")
p1.plot(sens_cat4['WAVELENGTH'],sens_cat4['SENSITIVITY'],label="Euclid (TIPS)")
#p1.plot(pzcaletc_1n['col1']*1e4,1./(pzcaletc_1n['col%i' % (n)]*1000.)) # erg/s/cm^2
#p1.plot(pzcaletc_8n['col1']*1e4,1./(pzcaletc_8n['col%i' % (n)]*1000.)) # erg/s/cm^2
p1.set_xlim(8750,20750)
p1.set_ylabel("e-/s per erg/s/cm$^2$/Ang") # according to axe_manual
p1.legend()


# In[ ]:


from scipy import integrate

### Plot 1D spectra
fig = plt.figure(figsize=(6,18))
p1 = fig.add_subplot(311)
p2 = fig.add_subplot(312)
p3 = fig.add_subplot(313)
for i, key in enumerate(beams.keys()):
    print(key)
    print()
    w, f, e = beams[key].beam.optimal_extract(beams[key].model+beams[key].grism.data['SCI'], bin=0)
    #w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], 
    #                                          ivar=1./(beams[key].grism.data['ERR'])**2,bin=0)
    
    # normalize to the magnitude in the direct image
    print("sum = %.2e" % (np.nansum(f)))
    print("mean = %.2e" % (np.nanmean(f)))
    print("median = %.2e" % (np.nanmedian(f)))
    print("std = %.2e" % (np.nanstd(f)))
    print("min = %.2e" % (np.nanmin(f)))
    print("max = %.2e" % (np.nanmax(f)))
    
    
    print("S/N = %.2f" % (np.nanmax(f)/np.sqrt(np.nansum(e**2))))
    print()
    #

    
    # Method 1
    scale = f_scale/np.nansum(f)
    f = f * scale
    e = e * scale
    # Method 2
    #f = (f/np.nanmedian(f)) * f_scale
    # Method 3
    #skysub = f - np.nanmedian(f)
    #f = (skysub/np.nansum(skysub)) * f_scale
    
    print("sum = %.2e" % (np.nansum(f)))
    print("mean = %.2e" % (np.nanmean(f)))
    print("median = %.2e" % (np.nanmedian(f)))
    print("std = %.2e" % (np.nanstd(f)))
    print("min = %.2e" % (np.nanmin(f)))
    print("max = %.2e" % (np.nanmax(f)))
    print("S/N = %.2f" % (np.nanmax(f)/np.sqrt(np.nansum(e**2))))
    print()
    
    
    percent = [1,25,50,75,95,99,99.9]
    for per in percent:
        print("P(%s) = %.2e" % (per,np.nanpercentile(f, per)))
    print()
        
    N = len(f)
        
    print("Chunk  W0  W1  Sum  Mean  Median  STD  Min  Max  S/N")
    
    #chunks = 50. # depends on the width of the line?
    chunks = 70.
    for chunk in np.arange(chunks):
        #print("Chunk = %i" % (chunk+1))
        f_sect = f[int(chunk*N/chunks):int((chunk+1)*N/chunks)]
        e_sect = e[int(chunk*N/chunks):int((chunk+1)*N/chunks)]
        w_sect = w[int(chunk*N/chunks):int((chunk+1)*N/chunks)]
        w_sect0 = w_sect[0]/1.e4
        w_sect1 = w_sect[-1]/1.e4
        
        print("%3i  %.3f  %.3f  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %6.1f" % (chunk+1, w_sect0, w_sect1, np.nansum(f_sect), 
              np.nanmean(f_sect), np.nanmedian(f_sect), np.nanstd(f_sect), np.nanmin(f_sect), 
              np.nanmax(f_sect), np.nansum(f_sect)/np.sqrt(np.nansum(e_sect**2))))
        #print()
        
        
    #print("S/N =",np.nanmax(f)/np.nanstd(f))
    
    m_noi = np.nanmean(e)
    S_N = f/m_noi
    
    
    # plot spectrum
    p1.plot(w/1.e4, f, c="k") # linestyle='steps-mid')
    p1.plot(w/1.e4, e, c="r") # linestyle='steps-mid')
    p1.plot([w[0]/1e4,w[-1]/1e4],[m_noi,m_noi], c="g")

    # plot line
    p2.plot(w/1.e4, f, c="k") # linestyle='steps-mid')
    p2.plot(w/1.e4, e, c="r") # linestyle='steps-mid')
    #print(f)
    
    # S/N of line
    p3.plot(w/1.e4, S_N, c="b") # linestyle='steps-mid')
    
    # Line flux and S/N
    print()
    w1 = 1.24
    w2 = 1.26
    #w1 = 1.36
    #w2 = 1.37
    #w1 = 1.48
    #w2 = 1.52
    i1 = np.argmin(np.abs(w/1.e4-w1))
    i2 = np.argmin(np.abs(w/1.e4-w2))
    print(i1,i2)
    
    win = w2 - w1
    i0 = np.argmin(np.abs(w/1.e4-(w1-win)))
    i3 = np.argmin(np.abs(w/1.e4-(w2+win)))
    print(i0,i3)
    print()
    

    dw = w[1] - w[0]
    print(dw)
    
    
    line_flux = np.sum(f[i0:i1])*dw
    line_error = np.sqrt(np.sum((e[i0:i1]*dw)**2))
    print("Flux =", line_flux)
    print("Error =", line_error)
    print("S/N =", line_flux/line_error)
    print()
        
    
    line_flux = np.sum(f[i1:i2])*dw
    line_error = np.sqrt(np.sum((e[i1:i2]*dw)**2))
    print("Flux =", line_flux)
    print("Error =", line_error)
    print("S/N =", line_flux/line_error)
    
    
    print(integrate.trapz(f[i1:i2],w[i1:i2]))
    
    print()
    
    line_flux = np.sum(f[i2:i3])*dw
    line_error = np.sqrt(np.sum((e[i2:i3]*dw)**2))
    print("Flux =", line_flux)
    print("Error =", line_error)
    print("S/N =", line_flux/line_error)
    print()
        

    
    #y0 = np.nanmin(f)
    #y1 = np.nanmax(f)
    #print(y0,y1)
    
    #z = 1.726
    #lines = [4861.,5007.,6563.]
    #for line in lines:
    #    wobs = line/1e4*(1+z)
    #    p.plot([wobs,wobs],[y0,y1],"--",c="r")
    print()
    print("Size =", len(f))
    print("NaNs = ", np.sum(np.isnan(f)))
    #print(np.nanmedian(f))

#p.legend(fontsize=14)
p1.set_xlim(0.9, 2.05)
#p.set_ylim(0.0,0.025)
#p1.set_ylim(-2e-22,6e-22) # shows noise level
p1.set_xlabel(r'$\lambda$ [micron]')
p1.set_ylabel(r'F$_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')

#p.legend(fontsize=14)
p2.set_xlim(w1-0.05, w2+0.05)
p2.set_ylim(-1e-21,3e-20)
#p2.set_ylim(-1e-21,1e-21)
p2.set_xlabel(r'$\lambda$ [micron]')
p2.set_ylabel(r'F$_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')

p3.set_xlim(w1-0.05, w2+0.05)
p3.set_ylim(-20,120.)
p3.set_xlabel(r'$\lambda$ [micron]')
p3.set_ylabel(r'S/N')


# In[ ]:


sens_file = HOME_PATH + '/grizli/CONF/sens_0720_2020.fits'
sens_cat2 = Table.read(sens_file)
print(sens_cat2.colnames)

from scipy import interpolate

R1 = interpolate.interp1d(sens_cat1['WAVELENGTH'],sens_cat1['SENSITIVITY'])
w10,w11 = sens_cat1['WAVELENGTH'][0],sens_cat1['WAVELENGTH'][-1]

R2 = interpolate.interp1d(sens_cat2['WAVELENGTH'],sens_cat2['SENSITIVITY'])
w20,w21 = sens_cat2['WAVELENGTH'][0],sens_cat2['WAVELENGTH'][-1]

fig = plt.figure()
p1 = fig.add_subplot(111)

for i, key in enumerate(beams.keys()):
    print(key)
    print()
    w, f, e = beams[key].beam.optimal_extract(beams[key].model+beams[key].grism.data['SCI'], bin=0)
    #w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], 
    #                                          ivar=1./(beams[key].grism.data['ERR'])**2,bin=0)
    
    print(w)
    
    # Method 4
    
    i10 = np.argmin(np.abs(w-w10))
    i11 = np.argmin(np.abs(w-w11))
    
    i20 = np.argmin(np.abs(w-w10))
    i21 = np.argmin(np.abs(w-w11))
    
    w1 = w[i10:i11]
    f1 = f[i10:i11]
    e1 = e[i10:i11]
    
    w2 = w[i20:i21]
    f2 = f[i20:i21]
    e2 = e[i20:i21]
    
    response1 = [R1(w_el) for w_el in w1]
    
    response2 = [R2(w_el) for w_el in w2]
    
    flux1 = f1/response1
    noise1 = e1/response1
    
    flux2 = f2/response2
    noise2 = e2/response2
    
    
    # plot spectrum
    p1.plot(w1/1.e4, flux1, c="k") # linestyle='steps-mid')
    #p1.plot(w1/1.e4, noise1, c="r") # linestyle='steps-mid')
    p1.plot(w2/1.e4, flux2, c="g") # linestyle='steps-mid')
    #p1.plot(w2/1.e4, noise2, c="r") # linestyle='steps-mid')

#p.legend(fontsize=14)
p1.set_xlim(0.9, 2.05)
#p.set_ylim(0.0,0.025)
#p1.set_ylim(-2e-22,6e-22) # shows noise level
p1.set_xlabel(r'$\lambda$ [micron]')
p1.set_ylabel(r'F$_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')



# In[ ]:


#/Users/gwalth/Dropbox/Research/sims/Galacticus/new_spectra
#SPEC_PATH = "/Users/gwalth/Dropbox/Research/aXeSIM/Roman/aXeSIM_Roman/SIMDATA/new_spectra"
#SPEC_PATH = "/Users/gwalth/Dropbox/Research/sims/Galacticus/bkup_2022.0120/spectra"
SPEC_PATH = "/Users/gwalth/Dropbox/Research/sims/Galacticus/sed_galacticus_ForGreg"
# ATLAS number 2118, modspec=1752
#primer_number=2118 # old
primer_number=1571 

#specdat = SPEC_PATH + "/ATLAS_1deg_subsample_spec_%06d.dat" % (primer_number)
specdat = SPEC_PATH + "/ATLAS_1deg_spec_%06d.dat" % (primer_number)
specdata = Table.read(specdat,format="ascii",names=('wave','flux'))
print(specdata)

fig = plt.figure()
p1 = fig.add_subplot(111)
p1.plot(specdata["wave"],specdata["flux"],label="Input spectra")
# [Ang] [erg/s/cm^2/Ang]
p1.set_xlim(8750,20750)
p1.set_yscale("log")
#p1.set_ylabel("e-/s per erg/s/cm$^2$/Ang") # according to axe_manual
#p1.legend()


# In[ ]:


beam = 'beam__%05d.grism.A.fits' % (id)
#new_beam = '{0}_{1:05d}.beams.fits'.format(root, id)
new_beam = beam.replace(".fits","_GLW.fits")

hdu = pyfits.open(beam)
print(hdu[0].header)
hdu.info()


ext = 0
hdu[ext].header['EXPTIME'] = hdu['SCI'].header['EXPTIME']
hdu.writeto(new_beam,clobber=True)
hdu.info()


# In[ ]:


#mb = multifit.MultiBeam([old_beam], fcontam=0.2, group_name=root, psf=False, min_sens=0.05)
mb = multifit.MultiBeam([new_beam], fcontam=0.2, group_name=root, psf=False, min_sens=0.05)
mb.write_master_fits()

####################################################################################
# Limited set of red stellar templates
#tstar = grizli.utils.load_templates(fwhm=1200, line_complexes=True, 
#                                    fsps_templates=True, stars=True)

# Fit spectral types.  Just makes the plot now, outputs not saved 
#fig, result, tfit = mb.xfit_star(tstar=tstar, fit_background=False,
#                                 spline_correction=True, spline_args={'Rspline':5})
####################################################################################

#fwhm = 325 # km/s
fwhm = 650 # km/s

# Fitting templates

# First is set with combined emission line complexes for the redshift fit 
# (don't allow infinite freedom) of the line ratios / fluxes
t0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                     full_line_list=None,  continuum_list=None, 
                                     fsps_templates=True)

# Second set has individual line templates for fitting the line fluxes
t1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                     full_line_list=None, continuum_list=None, 
                                     fsps_templates=True)


#fit = mb.xfit_redshift(templates=t0, zr=[0.65, 1.6], dz=[0.004, 0.0002], fitter='nnls')

###################################################################################
#hdu, fig = mb.drizzle_grisms_and_PAs(fcontam=0.2, flambda=False, kernel='point', 
#                                     size=32, zfit=tfit, diff=True)


# In[ ]:


print(fit.keys())


# ## Velocity resolution
# [top](#Table-of-Contents)

# In[ ]:


# https://wfirst.ipac.caltech.edu/sims/Param_db.html#wfi_grism
# R = 461*wav # wav in microns (2pix)

# Euclid
# R = 380*wav # wav in microns

c = 3e5 # km/s

# R = wav/delta_wav
# u = delta_wav/wav*c
# u = c/R

# u = (z - z0)*c/(1+z0)
# u = dz*c/(1+z0)
# dz = u(1+z0)/c
# dz = delta_wav/wav * (1+z0)

 
R = 380 # Euclud
#R = 461 # Roman

u = lambda wav: c/(R*wav)

winc = 0.25
w = np.arange(1.0,2.0+winc,winc)
for w0 in w:
    print("v = %.2f km/s (%.2f micron)" % (u(w0),w0))

print()   
for w0 in w:
    print("R = %.2f (%.2f micron)" % (R*w0,w0))

print()
z = 1.
print("z = %.2f" % (z))
for w0 in w:
    print("dz = %.6f (%.2f micron)" % (u(w0)*(1+z)/c,w0))

    
# R = wav/delta_wav
# R = 461*wav # microns
# 1/delta_wav = 461
# delta_wav = 1/461. 
print()
print("delta_wav = %.7f microns" % (1/R))
print("delta_wav = %.3f Angstroms" % (1/R*10000.))




# In[ ]:


# dz = delta_wav/wav * (1+z0)

dz = 10/12400.* (1+0.88615)
print(dz)


# ## Fit redshift to source
# [top](#Table-of-Contents)

# In[ ]:


beam = 'beam__%05d.grism.A.fits' % (id)
#new_beam = '{0}_{1:05d}.beams.fits'.format(root, id)
new_beam = beam.replace(".fits","_GLW.fits")

hdu = pyfits.open(beam)
print(hdu[0].header)
hdu.info()


ext = 0
hdu[ext].header['EXPTIME'] = hdu['SCI'].header['EXPTIME']
hdu.writeto(new_beam,clobber=True)
hdu.info()


# In[ ]:


#mb = multifit.MultiBeam([old_beam], fcontam=0.2, group_name=root, psf=False, min_sens=0.05)
mb = multifit.MultiBeam([new_beam], fcontam=0.2, group_name=root, psf=False, min_sens=0.05)
mb.write_master_fits()

####################################################################################
# Limited set of red stellar templates
#tstar = grizli.utils.load_templates(fwhm=1200, line_complexes=True, 
#                                    fsps_templates=True, stars=True)

# Fit spectral types.  Just makes the plot now, outputs not saved 
#fig, result, tfit = mb.xfit_star(tstar=tstar, fit_background=False,
#                                 spline_correction=True, spline_args={'Rspline':5})
####################################################################################

#fwhm = 325 # km/s
fwhm = 650 # km/s

# Fitting templates

# First is set with combined emission line complexes for the redshift fit 
# (don't allow infinite freedom) of the line ratios / fluxes
t0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                     full_line_list=None,  continuum_list=None, 
                                     fsps_templates=True)

# Second set has individual line templates for fitting the line fluxes
t1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                     full_line_list=None, continuum_list=None, 
                                     fsps_templates=True)


#fit = mb.xfit_redshift(templates=t0, zr=[0.65, 1.6], dz=[0.004, 0.0002], fitter='nnls')

###################################################################################
#hdu, fig = mb.drizzle_grisms_and_PAs(fcontam=0.2, flambda=False, kernel='point', 
#                                     size=32, zfit=tfit, diff=True)


# In[ ]:


fitting.run_all(id, t0, t1, fit_only_beams=True, fwhm=fwhm, zr=[0.05, 3.0], 
                dz=[0.004, 0.0002], fitter=['nnls', 'bounded'], group_name=root)


# ## Compare to primer

# In[ ]:


print(root,id)

# Primer redshift
filt = primer["NUMBER"] == primer_number # ATLAS NUMBER
z_true = primer["REDSHIFT"][filt]
print("z_true = %.6f" % (z_true))


# Redshift fit to the spectrum
full_hdu = pyfits.open('{0}_{1:05d}.full.fits'.format(root, id))
head = full_hdu[0].header
z_fit = head['REDSHIFT']
print("z_fit  = %.6f" % (z_fit))

dz = z_true - z_fit
print("dz     = %.6f" % (dz))


# In[ ]:


#print(primer)
print(primer.keys())
modspec_number = 8127
#print(primer_number)
#print(primer[modspec_number-1])

filt = primer["NUMBER"] == primer_number # ATLAS NUMBER
print(primer[filt])
print(primer["REDSHIFT"][filt])

filt = primer["SPECTEMP"] == modspec_number # MODSPEC NUMBER
print(primer[filt])


# ## Coordinates check
# [top](#Table-of-Contents)

# In[ ]:


fig = plt.figure(figsize=[10,10])

mag_limit = 25

filt_pri = primer['MAG_F1600W'] < mag_limit
filt_sex = phot['MAG_AUTO'] < mag_limit

ax = fig.add_subplot(111)
ax.scatter(primer['RA'][filt_pri], primer['DEC'][filt_pri], s=20,
               edgecolor='green', facecolor='none', alpha=0.8, label="Primer")
ax.scatter(phot['X_WORLD'][filt_sex], phot['Y_WORLD'][filt_sex], s=100,
               edgecolor='orange', facecolor='none', alpha=0.7, label="SExtractor")

ax.set_xlabel("RA [deg]")
ax.set_ylabel("Dec [deg]")
ax.invert_xaxis()

ax.legend(loc=1)


# <span style="color:red">
# Different sims may have different WCS centers!
# </span>
# Use this test as a way to check the correct center was used.
# 
# 
# aXeSIM conf | hdf5
# :----------:|:---:
# Roman.G150.wcs_v1.RN0.conf | galacticus_ATLAS_1deg2_subsample.hdf5  
# Roman.G150.wcs_v2.RN0.conf | galacticus_ForGreg.hdf5

# ## SED check
# 
# [top](#Table-of-Contents)

# In[ ]:


#/Users/gwalth/Dropbox/Research/sims/Galacticus/new_spectra
#SPEC_PATH = "/Users/gwalth/Dropbox/Research/aXeSIM/Roman/aXeSIM_Roman/SIMDATA/new_spectra"
#SPEC_PATH = "/Users/gwalth/Dropbox/Research/sims/Galacticus/bkup_2022.0120/spectra"
SPEC_PATH = "/Users/gwalth/Dropbox/Research/sims/Galacticus/sed_galacticus_ForGreg"
# ATLAS number 2118, modspec=1752
#number = 2118 # old
number = 1571

#specdat = SPEC_PATH + "/ATLAS_1deg_subsample_spec_%06d.dat" % (number)
specdat = SPEC_PATH + "/ATLAS_1deg_spec_%06d.dat" % (number)
print(specdat)
specdata = Table.read(specdat,format="ascii",names=('wave','flux'))
print(specdata)

fig = plt.figure()
p1 = fig.add_subplot(111)
p1.plot(specdata["wave"],specdata["flux"],label="Input spectra")
# [Ang] [erg/s/cm^2/Ang]
p1.set_xlim(8750,20750)
p1.set_yscale("log")

p1.set_xlabel("Wavelength [Ang]")
#p1.set_ylabel("e-/s per erg/s/cm$^2$/Ang") # according to axe_manual
#p1.legend()


# In[ ]:


SPEC_PATH = "/Users/gwalth/Dropbox/Research/sims/Galacticus/sed_galacticus_ForGreg"
AXE_PATH = "/Users/gwalth/data/Roman/grizli/my_roman_sims/Prep/"

axe_image = "Roman_ATLAS_1deg_random2022_2022-01-27T06:04:27_RN0_v6_images.fits"
axe_spec = "Roman_ATLAS_1deg_random2022_2022-01-27T06:04:27_RN0_v6_spectra.fits"


# In[ ]:


f = AXE_PATH + axe_spec
pf = pyfits.open(f)
N = len(pf)

#number = 2118 # old
number = 1571

for i in np.arange(N)+1:
    
    head = pf[i].header
    
    specname = head['SPECNAME']
    extname = head['EXTNAME']
    mag_ab = head['MAG_AB']
    
    number_found = int(specname.split("_")[-1])
    if number == number_found:
        print(i,number_found)
    
    # SPECNAME= 'ATLAS_1deg_spec_029260' / Name of spectrum 

    #tab = pf[ext].data
    #flux = tab["flux"]
    #wav = tab["wav_nm"]


# In[ ]:


f = AXE_PATH + axe_spec
#f = axe_spec
#ext = 4468 # old
ext = 8127
pf = pyfits.open(f)
print(len(pf))

#print(key)

tab = pf[ext].data

flux = tab["flux"]
wav = tab["wav_nm"]


delta_wav = 1/461. # microns
delta_wav *= 1000. # nm


################################################
# convolve with the resolution of the instrument
################################################
# from iris_snr_sim (essentially from Tuan Do)

delt = 2.0*(delta_wav)/(wav[1]-wav[0])


stddev = delt/2*sqrt(2*log(2))
psf_func = models.Gaussian1D(amplitude=1.0, stddev=stddev)
x = np.arange(4*int(delt)+1)-2*int(delt)
psf = psf_func(x)
psf /= psf.sum() # normalize

new_flux = np.convolve(flux, psf,mode='same')


fig = plt.figure()
p1 = fig.add_subplot(111)
p1.plot(wav,flux,label="aXeSIM spectra",c="r",alpha=0.5)
p1.plot(wav,new_flux,label="Convolved spectra",c="g",alpha=0.8)
# [nm] [erg/s/cm^2/Ang]

#######################
#######################
#######################

# ATLAS number 2118, modspec=1752
#number = 2118 # old
number = 1571 # new

#specdat = SPEC_PATH + "/ATLAS_1deg_subsample_spec_%06d.dat" % (number)
specdat = SPEC_PATH + "/ATLAS_1deg_spec_%06d.dat" % (number)
print(specdat)
specdata = Table.read(specdat,format="ascii",names=('wave','flux'))
print(specdata)

p1.plot(specdata["wave"]/10.,specdata["flux"],label="Input spectra")
# [Ang] [erg/s/cm^2/Ang]

######################
sens_file = HOME_PATH + '/grizli/CONF/sens_0720_2020.fits'
sens_cat2 = Table.read(sens_file)
print(sens_cat2.colnames)

from scipy import interpolate

R = interpolate.interp1d(sens_cat2['WAVELENGTH'],sens_cat2['SENSITIVITY'])
w0,w1 = sens_cat2['WAVELENGTH'][0],sens_cat2['WAVELENGTH'][-1]

w, f, e = beams[key].beam.optimal_extract(beams[key].model+beams[key].grism.data['SCI'], bin=0)
#w, f, e = beams[key].beam.optimal_extract(beams[key].grism.data['SCI'], 
#                                          ivar=1./(beams[key].grism.data['ERR'])**2,bin=0)
    
#print(w)
    
# Method 4

i0 = np.argmin(np.abs(w-w0))
i1 = np.argmin(np.abs(w-w1))

w = w[i0:i1]
f = f[i0:i1]
e = e[i0:i1]
    
response = [R(w_el) for w_el in w]
    
flux2 = f/response
noise2 = e/response
    
# plot spectrum
p1.plot(w/10., flux2, c="orange", label="extracted spectra") # linestyle='steps-mid')
# factor of 10 off, could be WAV_AB being wrong (i.e. 1600 instead of 160)

p1.set_xlim(875.0,2075.0)
#p1.set_ylim(1e-15,1e-21)
p1.set_yscale("log")
p1.set_xlabel("Wavelength [nm]")
#p1.set_ylabel("e-/s per erg/s/cm$^2$/Ang") # according to axe_manual

p1.legend(loc=1)



pf[ext].header


# In[ ]:


f = AXE_PATH + axe_image

ext = 3
pf = pyfits.open(f)
print(len(pf))
pf[ext].header

#tab = pf[ext].data

#fig = plt.figure()
#p1 = fig.add_subplot(111)
#p1.plot(tab["wav_nm"],tab["flux"],label="Input spectra")
## [Ang] [erg/s/cm^2/Ang]
#p1.set_xlim(875.0,2075.0)
#p1.set_yscale("log")
##p1.set_ylabel("e-/s per erg/s/cm$^2$/Ang") # according to axe_manual
##p1.legend()


# In[ ]:


conf_path = "/Users/gwalth/data/Roman/grizli/grizli/templates/fsps/"

#L = glob.glob(conf_path + "*_v3_nolines_???.dat")
L = glob.glob(conf_path + "*_v3_???.dat")

print(L)

fig = plt.figure()
ax1 = fig.add_subplot(111)

for l in L:
    tbl = Table.read(l, format="ascii")
    #print(tbl)

    ax1.plot(tbl["wave"],tbl["flux"])
    
ax1.set_xlim(3000,12500)
ax1.set_ylim(1e-7,1e-1)
ax1.set_yscale("log")
ax1.set_xlabel("Wavelength [Ang]")
plt.show()


# In[ ]:




