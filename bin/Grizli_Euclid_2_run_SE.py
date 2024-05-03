#!/usr/bin/env python

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
#yaml_file = "config.yaml"
yaml_file = sys.argv[1]
det_img = "Euclid-total_ref.fits"
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


# remove segmentation FITS to redo "next" step
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
os.system('rm *_seg.fits')



## Make SExtractor catalog
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

# Euclid (Schirmer et al. 2022, XVIII. The NISP photometric system)
#mag_zero = 25.04 # +- 0.05 mag   Y 
#mag_zero = 25.26 # +- 0.05 mag   J
#mag_zero = 25.21 # +- 0.05 mag   H
ref_files = yaml_dict["ref_files"]
           

#           VIS   Y      J      H      total
mag_zero = yaml_dict["mag_zero"]



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
        
yaml_dict["all_ref_files"] = all_ref_files
yaml_dict["all_cat"] = all_cat
yaml_dict["all_seg"] = all_seg
yaml_dict["all_bkg"] = all_bkg

yaml_dict["phot_mode"] = "SExtractor"

os.chdir(YAML_PATH)
with open(os.path.basename(yaml_file), 'w',) as f:
    yaml.dump(yaml_dict, f, sort_keys=False)
