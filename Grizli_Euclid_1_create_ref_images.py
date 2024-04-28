### Euclid Parameters and Requirements
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


#################################################################################
# usage:
#
# python Grizli_Euclid_1_create_ref_images.py
# python Grizli_Euclid_2_run_SExtractor.py ../sims/Euclid/TestPoints/config.yaml
# python Grizli_Euclid_3_prep.py ../sims/Euclid/TestPoints/config.yaml
# python Grizli_Euclid_4_model.py ../sims/Euclid/TestPoints/config.yaml
# python Grizli_Euclid_5_fit_redshifts.py DET11 ../sims/Euclid/TestPoints/config.yaml
#
#
#
#
# python Grizli_Euclid_2_run_SEP.py ../sims/Euclid/TestPoints/config.yaml
#
#
#
#################################################################################

import time
T_START = time.time()

import grizli_functions
#from grizli_functions import wcs_pixel_scale, check_sims, create_circular_mask
from grizli_functions import (
        add_noise, 
        wcs_pixel_scale, 
        check_sims, 
        display_grizli,
        fake_euclid_direct,
        fake_euclid_ref,
        total_image,
        euclid_det, 
        read_slitless_headers,
        write_individual_slitless, 
        plot_slitless, 
        map_src_to_det,
        calc_photflam,
)
print(grizli_functions.__file__)
# import jwst failure is ok!

import yaml
import glob, os, sys
#from collections import OrderedDict

import matplotlib as mpl    
import matplotlib.pyplot as plt
#from matplotlib.gridspec import GridSpec
#from matplotlib.ticker import MultipleLocator

#from IPython.display import Image

mpl.rcParams['figure.figsize'] = (10.0, 6.0)
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.dpi'] = 72

import numpy as np
#from scipy import integrate
#from math import cos, sin, atan2, pi
#from math import sqrt, log


import astropy
#import astropy.units as u
#from astropy.coordinates import SkyCoord

import astropy.io.fits as pyfits
#from astropy import wcs
from astropy.table import Table, unique, join, vstack
#from astropy.modeling import models

#import drizzlepac
#import photutils

#import grizli
#import grizli.model
#import grizli.multifit
#from grizli import utils, multifit, fitting
#import grizli.fake_image
#from grizli.pipeline import auto_script
#from grizli import prep

print('\n Python version: ', sys.version)
#print('\n Grizli version: ', grizli.__version__)
print('\n Astropy version: ', astropy.__version__)


####################################
# paramters
####################################
plot = 1
mag_zero = [25.6, 25.04, 25.26, 25.21, 26.0]
yaml_file = f"config.yaml"


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
####################################


# ## Install templates for redshift fitting
# [top](#Table-of-Contents)
# 
# Run only once for the install

#### WFC3 and ACS calibs
#grizli.utils.fetch_default_calibs()

#### WFC3 PSF and Pickles stars
#grizli.utils.fetch_config_files()

#### Templates used in fitting
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


#YAML_PATH = os.getcwd()

#os.chdir('../')
#os.chdir('/Users/gwalth/data/Roman/grizli/sims/')
os.chdir('/Users/gwalth/data/Roman/grizli/sims/Euclid')

#os.chdir('/local/RomanSims/grizli/sims/') # cygnusd
HOME_PATH = os.getcwd()
print('HOME_PATH = ', HOME_PATH)
#root = "SIM_10_18_22"
#root = "SIM_12_23_22"
#root = "TestPoints"
#root = "TestPoints_v2"
#root = "TestPoints_v3"

#root = "TUTESTGAL"
#root = "Env9"
#root = "Env9_V1_RN"
#root = "Env9_V1_TEST"
root = "FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11"



YAML_PATH = os.path.join(HOME_PATH, root)

os.chdir(os.path.join(HOME_PATH, root, 'Prep'))



#slitless_files = glob.glob('EUC_SIM_*.fits')
#slitless_files.sort()
#print(len(slitless_files))

# First test
#slitless_files = ['NISPS_TIPS_TestPoints_highSNR_mod1_14324_2023_05_26_frame1.fits']
#zodi_files = ['NISPS_TIPS_TestPoints_highSNR_JustZodi_14324_2023_05_26_frame1.fits']
#catalog_files = ['TestPoints_highSNR_mod1_14324.fits']

# TUTESTGAL
#slitless_files = ['NISPS_TIPS_TUTESTGAL_noSTARS_RN_frame1.fits']
#zodi_files = []
#catalog_files = ['TUGALCAT_TEST.fits']


# Env9
#slitless_files = ['NISPS_TIPS_Env9_mod3_14826_Hcut23_REFMAGeq1_noSTARS_RN_frame1.fits']
#zodi_files = []
#catalog_files = ['Env9_mod3_14826_Hcut23_REFMAGeq1.fits']

# Env9 V1
#slitless_files = ['Env9_V1_TestArray_RN_2023_12_15_frame1.fits']
#zodi_files = []
#catalog_files = ['Env9_mod3_15329_base_Felgt2em16_Visit1_TESTarray.fits']

# Env9 V1 TEST
#slitless_files = [
#    'NISP_SLITLESS_FRAME1_Env9_V1_TESTarray_noMW_RNDC_2024-02-05.fits',
#    'NISP_SLITLESS_FRAME2_Env9_V1_TESTarray_noMW_RNDC_2024-02-05.fits',
#    'NISP_SLITLESS_FRAME3_Env9_V1_TESTarray_noMW_RNDC_2024-02-05.fits',
#    'NISP_SLITLESS_FRAME4_Env9_V1_TESTarray_noMW_RNDC_2024-02-05.fits',
#]
#zodi_files = []
#catalog_files = ['Env9_mod3_15329_base_Felgt2em16_Visit1_TESTarrayNEW.fits']

# FSpatch_mod3
slitless_files = [
     'NISP_SLITLESS_FRAME1_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
     #'NISP_SLITLESS_FRAME2_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
     #'NISP_SLITLESS_FRAME3_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
     #'NISP_SLITLESS_FRAME4_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
]
#slitless_files = [
#     'NISP_SLITLESS_FRAME1_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
#     'NISP_SLITLESS_FRAME2_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
#     'NISP_SLITLESS_FRAME3_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
#     'NISP_SLITLESS_FRAME4_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits',
#]
zodi_files = []
catalog_files = ['FSpatch_mod3_16183_TESTarrayCALIB_V1.fits']


print(slitless_files)
print(catalog_files)


#slitless_files = glob.glob("EUC_SIM_NISR*.fits")
#catalog_files = glob.glob("CATALOG*.fits")
#print(slitless_files)
#print(catalog_files)


# ## Directory Structure
# 
# I was structing it similar to Grizli with Prep, RAW and Extraction directories.  If this were real mission data, the stage that we recieved from Anihita would have been drizzled images and spectra which would go into the Prep directories.
# 
# This is just showing that we have the right directories and we can find all of the files.

# Clean
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

os.system('rm *_direct.fits')
os.system('rm *_slitless.fits')
os.system('rm *_wcs.fits')
os.system('rm *_final*.fits')
os.system('rm *_final.cat')
os.system('rm Euclid_GrismFLT.pickle')


# detector header names
all_det = euclid_det()


## Write individual files for each extension of the slitless spectra
# Grizli is easier to manage when writing out all of the files. 
# At some point we'll want to read the data extensions directly into Grizli, 
# this is currently a kludge.
#all_slitless = ["Euclid_FRAME%i" % (i+1) + "_DET%s_slitless.fits" for i,sf in enumerate(slitless_files)]
all_slitless = [write_individual_slitless(sf, file_str = "Euclid_FRAME%i" % (i+1) + "_DET%s_slitless.fits", rot=1) for i,sf in enumerate(slitless_files)]
all_zodi = [write_individual_slitless(sf, file_str ="Zodi_DET%s_slitless.fits") for sf in zodi_files]

os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

if plot:
    file = "Euclid_FRAME1_DET11_slitless.fits"
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


# Plot the slitless extensions
if plot:
    plot_slitless(slitless_files[0], vmin=500, vmax=1700, verb=0)


## Read slitless headers and plot the image coordinates of the detectors relative to each other
heads = read_slitless_headers(slitless_files[0], verb=0, plot=0)
print(heads)


# ## Read the source catalog
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

print(catalog_files)
primer = Table.read(catalog_files[0]) 
print(primer.colnames)
print(len(primer))


filt = primer['RA'] < 200.
print(filt)
print(primer[filt])

index = np.arange(len(primer))
print(index[filt])
# for TestPoints
#primer.remove_row(index[filt][0])

#filt = primer['VIS'] == -99
#print(primer[filt])

#filt = primer['NIR_Y'] == -99
#print(primer[filt])

#filt = primer['NIR_J'] == -99
#print(primer[filt])

#filt = primer['NIR_H'] == -99
#print(primer[filt])


print([col for col in primer.colnames if "TU_" in col])
Euclid_bands = ['VIS','NIR_Y','NIR_J','NIR_H']
Euclid_bands_flux = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG'] 

# add magnitudes band into catalog
for bk,fk in zip(Euclid_bands,Euclid_bands_flux):
    fnu_Jy = primer[fk] # Jy    
    mab = -2.5*np.log10(fnu_Jy) + 8.90   # Jy --> AB mag            
    mab[np.isinf(mab)]=-99.
    primer[bk] = mab    


#primer[:10].show_in_notebook()
#primer[Euclid_bands][:10].show_in_notebook()
#primer[Euclid_bands_flux][:10].show_in_notebook()

# ## Plot the magnitude histogram of sources

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

if plot:
    fig = plt.figure()

    ax1 = fig.add_subplot(111)
    ax1.scatter(primer['RA'],primer['DEC'],s=0.1, alpha=1.0, edgecolors="none")
    ax1.set_aspect(1.)
    ax1.set_xlabel("RA [deg]")
    ax1.set_ylabel("Dec [deg]")

    plt.show()

os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))
#direct_thumbnail_files = glob.glob('NIS_catalog_file_??.thm.beamA.fits')
direct_thumbnail_files = glob.glob('catalog_??_??.1.thm.beamA.fits')
direct_thumbnail_files.sort()
print(direct_thumbnail_files)
print(len(direct_thumbnail_files))

#os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))
## Map the sources the sources from the catalog to the detectors 
#det_tbl, det_dict = map_src_to_det(plot=0)
#det_dict = map_src_to_det(primer, heads, plot=0)


# ## Print each of the detectors WCS header information

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

#primer["TU_FNU_H_NISP_MAG","NIR_H"][:10].show_in_notebook()

#    # effective wavelength (pivot wavelength)
#    # http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/
#    numerator = integrate.trapz(tran,wav)
#    denominator = integrate.trapz(tran/wav**2,wav)
#    effwav_sq = numerator/denominator
#    effwav = np.sqrt(effwav_sq)



# ## Create input table for the reference image

#thumb_temp = 'NIS_catalog_file_%s.thm.beamA.fits'
thumb_temp = 'catalog_%s_%s.1.thm.beamA.fits'
os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))

source_ids = []
thumbs = []

t0 = time.time()

k = 0
for i in range(0,4): 
    for j in range(0,4):
        print(i+1,j+1)
        num = '%i%i' % (i+1,j+1)
        #hdus = pyfits.open(thumb_temp % (num))
        hdus = pyfits.open(thumb_temp % (num,num))
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

#all_tbl[:10].show_in_notebook()

#print(primer)
print(primer.colnames)


# ## Create reference images

# Mauri et al. 2020
# (ropper et al. 2016
os.chdir(os.path.join(HOME_PATH, root, 'Prep', 'Input_Thumbnails'))

# noise characteristics of "direct image"
# Scaramella et al. 2022 - Euclid preparation I. The Euclid Wide Survey
# RGS000, RGS180, RGS000_rot, RGS180_rot
spec_exptime = 574 # seconds
#spec_gain = 2.0
spec_gain = 6.0


nexp = 4 

# VIS 
readnoise = 4.5  # e-
background = 1.2 # e-/s
#pixel_size = 12. # um
dir_gain = 1.0   # e-/photon

# VIS
flux_key = "TU_FNU_VIS_MAG"
mag_key = "VIS"
id_key = "SOURCE_ID"
dir_exptime = 570 # seconds
wav_cen = 7102.613 # calculated from pivot
wav_width = 9000. - 5600.
eff_tot = 0.70
#eff_tot = 0.74
output = "Euclid-VIS_ref.fits"
filt = "VIS"
instr = "VIS"
ZP = 25.6
photflam = calc_photflam(ZP, wav_cen) # photplam


t0 = time.time()

print("Creating VIS reference image")
test_fluxes_vis, test_mags_vis, test_offset = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, id_key=id_key, 
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output, 
                                         filt=filt, instr = instr, photflam=photflam)



t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))


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
id_key = "SOURCE_ID"
dir_exptime = 112. # seconds
wav_cen = 10809. # Ang
wav_width = 2627. # Ang
eff_tot = 0.772
output = "Euclid-NISP_Y_ref.fits"
filt = "NISP_Y"
instr = "NISP"
ZP = 25.04
photflam = calc_photflam(ZP, wav_cen) # photplam


print("Creating NISP_Y reference image")
test_fluxes_nisp_y, test_mags_nisp_y, _ = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, id_key=id_key,
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output,
                                         filt=filt, instr = instr, photflam=photflam)

t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))


t0 = time.time()
# J
flux_key = "TU_FNU_J_NISP_MAG"
mag_key = "NIR_J"
id_key = "SOURCE_ID"
dir_exptime = 112. # seconds
wav_cen = 13673. # Ang
wav_width = 3994. # Ang
eff_tot = 0.790
output = "Euclid-NISP_J_ref.fits"
filt = "NISP_J"
instr = "NISP"
ZP = 25.26
photflam = calc_photflam(ZP, wav_cen) # photplam


print("Creating NISP_J reference image")
test_fluxes_nisp_j, test_mags_nisp_j, _ = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, id_key=id_key,
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output,
                                         filt=filt, instr = instr, photflam=photflam)
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))
# In[ ]:



t0 = time.time()
# H
flux_key = "TU_FNU_H_NISP_MAG"
mag_key = "NIR_H"
id_key = "SOURCE_ID"
wav_cen = 17714. # Ang
wav_width = 4999. # Ang
dir_exptime = 112. # seconds
eff_tot = 0.782
output = "Euclid-NISP_H_ref.fits"
filt = "NISP_H"
instr = "NISP"
ZP = 25.21
photflam = calc_photflam(ZP, wav_cen) # photplam

print("Creating NISP_H reference image")
test_fluxes_nisp_h, test_mags_nisp_h, _ = fake_euclid_ref(final_tbl, ra_cen = ra_avg, dec_cen = dec_avg, 
                                         pixel_scale = 0.3, flux_key=flux_key, mag_key=mag_key, id_key=id_key,
                                         gain=dir_gain, background=background, exptime=dir_exptime, 
                                         nexp=nexp, readnoise=readnoise, wav_cen=wav_cen, 
                                         wav_width=wav_width, eff_tot=eff_tot, output=output,
                                         filt=filt, instr = instr, photflam=photflam)
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))

print(len(test_offset))




os.chdir(os.path.join(HOME_PATH, root, 'Prep'))


coord_offset_tbl = Table(rows=test_offset)
#print(coord_offset_tbl)
# in pixels
coord_offset_tbl.write("coord_offset.fits", overwrite=True)


t0 = time.time()
ref_files = ["Euclid-VIS_ref.fits", "Euclid-NISP_Y_ref.fits", "Euclid-NISP_J_ref.fits",
             "Euclid-NISP_H_ref.fits"]
total_image(ref_files, output="Euclid-total_ref.fits", img_ext='REF') 
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))

T_END = time.time()

print()
print("Finished in %.1f seconds" % (T_END-T_START))

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

yaml_dict = {
    'HOME_PATH': HOME_PATH,
    'root': root,
    'ref_files': ref_files,
    'mag_zero': mag_zero,
    'slitless_files': slitless_files,
    'zodi_files': zodi_files,
    'catalog_files': catalog_files,
    'all_slitless': all_slitless,
    'all_zodi': all_zodi,
    'spec_exptime': spec_exptime,
    'spec_gain': spec_gain,
    'nexp': nexp,         
}

os.chdir(YAML_PATH)
with open(yaml_file, 'w',) as f:
    yaml.dump(yaml_dict, f, sort_keys=False)
