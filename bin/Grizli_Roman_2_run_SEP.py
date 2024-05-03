#!/usr/bin/env python

import grizli_functions

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


############################################################################################
# WORKS! but doesn't it never gets to filter_tab
############################################################################################
# cd test1
#root = "hlss"
#phot = auto_script.multiband_catalog(field_root=root, filters=['F158'], sci_image="hlss-000.0-f158_drz_sci.fits",
#                                     detection_filter='F158') #, get_all_filters=True)
############################################################################################
############################################################################################
## cd test2
#root = "hlss"
#phot = auto_script.multiband_catalog(field_root=root, filters=['F158'], sci_image="hlss-F158_drz_sci.fits",
#                                     detection_filter='F158') 
############################################################################################
# HST
############################################################################################
# cd test3
root = "j002836m3311"
phot = auto_script.multiband_catalog(field_root=root, filters=['f110w','f160w'],
                                     detection_filter='f160w') 
############################################################################################



sys.exit()

