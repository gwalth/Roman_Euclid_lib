#!/usr/bin/env python
#
import yaml
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
)
print(grizli_functions.__file__)
# import jwst failure is ok!

import glob, os, sys
#from collections import OrderedDict

import matplotlib as mpl    
import matplotlib.pyplot as plt
#from matplotlib.gridspec import GridSpec
#from matplotlib.ticker import MultipleLocator

#from IPython.display import Image

mpl.rcParams['figure.figsize'] = (10.0, 6.0)
mpl.rcParams['font.size'] = 14
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

import grizli
#import grizli.model
#import grizli.multifit
#from grizli import utils, multifit, fitting
#import grizli.fake_image
from grizli.pipeline import auto_script
#from grizli import prep

print('\n Python version: ', sys.version)
print('\n Grizli version: ', grizli.__version__)
print('\n Astropy version: ', astropy.__version__)


####################################
# paramters
####################################
plot = 1
#yaml_file = "config.yaml"
yaml_file = sys.argv[1]
det_img = "Euclid_total_ref.fits"
####################################



#YAML_PATH = os.getcwd()


with open(yaml_file, 'r') as f:
    yaml_dict = yaml.safe_load(f)
    print(yaml_dict)

HOME_PATH = yaml_dict["HOME_PATH"]
print("HOME_PATH =", HOME_PATH)
root = yaml_dict["root"]
print("root =", root)
YAML_PATH = os.path.join(HOME_PATH, root)


## Make SExtractor catalog
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

# Euclid (Schirmer et al. 2022, XVIII. The NISP photometric system)
#mag_zero = 25.04 # +- 0.05 mag   Y 
#mag_zero = 25.26 # +- 0.05 mag   J
#mag_zero = 25.21 # +- 0.05 mag   H
ref_files = yaml_dict["ref_files"]

           

#           VIS   Y      J      H      total
mag_zero = yaml_dict["mag_zero"]



#all_cat = []
#all_seg = []
#all_bkg = []

all_ref_files = ref_files + [det_img]
print(all_ref_files)



#from importlib import reload
#from grizli import model, multifit, grismconf







### Run SEP (~SExtractor clone) catalog on the "ir" combined image
### and generate a photometric catalog with aperture photometry in all available bands
#if not os.path.exists('{0}_phot.fits'.format(root)):
#    
#    multiband_catalog_args=kwargs['multiband_catalog_args']
#    tab = auto_script.multiband_catalog(field_root=root,
#                                    **multiband_catalog_args)
#
##    get_background=False # SExtractor background subtraction
#
##     tab = auto_script.multiband_catalog(field_root=root, threshold=1.8,
##                                         detection_background=get_background,
##                                         photometry_background=get_background) 
#    
#files = glob.glob('{0}-ir*'.format(root)) + glob.glob('*phot*fits')
#for file in files:
#    print(file)
#    
#phot = utils.GTable.gread('{0}_phot.fits'.format(root))
#print('{0}Metadata{0}'.format('\n'+'='*20+'\n'))
#for k in phot.meta:
#    print('{0}:\t{1}'.format(k, phot.meta[k]))


def split_euclid_fits_files(
        fits_file, 
        keys = ['REF', 'ERR'], 
        nkeys = ['sci', 'wht'], 
        suffix = "_ref",
    ):

    hdu = pyfits.open(fits_file)
    #head = hdu[0].header
        
    for k,nk in zip(keys, nkeys):    
        img = hdu[k].data
        head = hdu[k].header

        new_fits_file = fits_file.replace(suffix + ".fits","_drz_" + nk + ".fits")
        
        new_hdu = pyfits.PrimaryHDU(img, header=head)
        new_hdu.writeto(new_fits_file, overwrite=True, output_verify='fix')
        
        print("Writing",new_fits_file)

#os.symlink(src,dst)


all_ref_files = [
    "Euclid-VIS_ref.fits",
    "Euclid-NISP_Y_ref.fits",
    "Euclid-NISP_J_ref.fits",
    "Euclid-NISP_H_ref.fits",
    "Euclid-total_ref.fits",
]

# test case
#all_ref_files = [all_ref_files[-1]]
if 1:
    for i,ref in enumerate(all_ref_files):
        print(ref)
        split_euclid_fits_files(ref)
        #split_fits_files_v1(ndf)
        #split_fits_files_v2(ndf)

#sys.exit()

#print(grizli.prep.SEXTRACTOR_PHOT_APERTURES)
#print(grizli.prep.SEXTRACTOR_PHOT_APERTURES_ARCSEC)

#pixel_scale = 0.06
#SEXTRACTOR_PHOT_APERTURES_ARCSEC = [float(ap)*pixel_scale*u.arcsec for ap in SEXTRACTOR_PHOT_APERTURES.split(',')]
#print(SEXTRACTOR_PHOT_APERTURES_ARCSEC)

#pixel_scale = 0.3
#SEXTRACTOR_PHOT_APERTURES_ARCSEC = [float(ap)*pixel_scale*u.arcsec for ap in SEXTRACTOR_PHOT_APERTURES.split(',')]
#print(SEXTRACTOR_PHOT_APERTURES_ARCSEC)


prefix = "Euclid"
#utils.set_warnings()
#phot = auto_script.multiband_catalog(field_root=root, filters=['F158'], sci_image="hlss-000.0-f158_drz_sci.fits",
#                                     detection_filter='F158', get_all_filters=True)
#phot = auto_script.multiband_catalog(field_root=prefix, filters=['NISP_Y', 'NISP_J', 'NISP_H'],
#                                     detection_filter='total',) # get_all_filters=True)
#phot = auto_script.multiband_catalog(field_root=prefix, filters=['total'], sci_image="Euclid-total_drz_sci.fits",
#                                     detection_filter='total',) # get_all_filters=True)
#phot = auto_script.multiband_catalog(field_root=prefix, filters=['total'],
#                                     detection_filter='total',) # get_all_filters=True)

phot = auto_script.multiband_catalog(field_root=prefix, filters=['VIS', 'NISP_Y', 'NISP_J', 'NISP_H'],
#                                     detection_filter='NISP_H',) # get_all_filters=True)
                                     detection_filter='total',) # get_all_filters=True)


yaml_dict["all_ref_files"] = all_ref_files
yaml_dict["all_cat"] = glob.glob("*phot.fits")
yaml_dict["all_seg"] = glob.glob("*seg.fits")
yaml_dict["all_bkg"] = glob.glob("*bkg.fits")

yaml_dict["phot_mode"] = "SEP"

os.chdir(YAML_PATH)
with open(os.path.basename(yaml_file), 'w',) as f:
    yaml.dump(yaml_dict, f, sort_keys=False)

#sys.exit()


