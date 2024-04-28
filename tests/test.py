import numpy as np
from astropy.table import Table, unique, join, vstack

ZP = np.array([25.6, 25.04, 25.26, 25.21])


# wav_cen or pivot


photplam = np.array([
    7102.612795438403,      # Euclid VIS GLW 4/13/23
    10785.387735681754,  # Euclid NISP GLW 4/13/23
    13620.676219329825,  # Euclid NISP GLW 4/13/23
    17648.78006935212,   # Euclid NISP GLW 4/13/23
])


#photplam = np.array([ 7102.613 , 10809., 13673., 17714.]) # Ang


#  ZP = (-2.5*np.log10(model.photflam_list[fi]) - 21.10 -
#                  5*np.log10(model.photplam_list[fi]) + 18.6921)

# ZP = -2.5*np.log10(photflam) - 21.10 - 5*np.log10(photplam) + 18.6921

photflam = 10**(-0.4*(ZP + 21.10 - 18.6921 + 5*np.log10(photplam)))


print(photflam)

#import grizli.prep
import astropy.units as u
from grizli.prep import SEXTRACTOR_PHOT_APERTURES

print(SEXTRACTOR_PHOT_APERTURES)

pixel_scale = 0.06
SEXTRACTOR_PHOT_APERTURES_ARCSEC = [float(ap)*pixel_scale*u.arcsec for ap in SEXTRACTOR_PHOT_APERTURES.split(',')]
print(SEXTRACTOR_PHOT_APERTURES_ARCSEC)

pixel_scale = 0.3
SEXTRACTOR_PHOT_APERTURES_ARCSEC = [float(ap)*pixel_scale*u.arcsec for ap in SEXTRACTOR_PHOT_APERTURES.split(',')]
print(SEXTRACTOR_PHOT_APERTURES_ARCSEC)


ee_file = "/Users/gwalth/python/src/grizli_v1.6.0.dev42/grizli/data/hst_encircled_energy.fits"
ee = Table.read(ee_file, format="fits") # ref_cat in multimission
print(ee.colnames)
print(ee)




