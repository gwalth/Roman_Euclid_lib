#!/usr/bin/env python

# python Grizli_Euclid_5_single_fit_redshift.py DET11 214.6623540 57.8587647 ../sims/Euclid/TestPoints_v2/config.yaml


# python Grizli_Euclid_5_single_fit_redshift.py DET11 213.5306518 57.4122003 ../sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/config.yaml




import time
import yaml
import os, sys
import logging 

import numpy as np

import grizli
import grizli.utils
import grizli.multifit
from grizli import utils, multifit, fitting

from grizli_functions import check_sims2, euclid_det


#logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)
logging.basicConfig(filename='Grizli_Euclid_5_fit_redshifts.log', encoding='utf-8')

####################################
# paramters
####################################
det = sys.argv[1]
a0 = float(sys.argv[2])
d0 = float(sys.argv[3])
yaml_file = sys.argv[4]

print(sys.argv)

mag_limit = 30
#fwhm = 395 # km/s
fwhm = 400 # km/s
####################################

#YAML_PATH = os.getcwd()
#print(YAML_PATH)
with open(yaml_file, 'r') as f:
    yaml_dict = yaml.safe_load(f)
    print(yaml_dict)

HOME_PATH = yaml_dict["HOME_PATH"]
print("HOME_PATH =", HOME_PATH)
root = yaml_dict["root"]
print("root =", root)
YAML_PATH = os.path.join(HOME_PATH, root)

t0 = time.time()

import pickle

os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

with open('Euclid_%s_GrismFLT.pickle' % (det), 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    Euclid = pickle.load(f)
    
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))



# bright source near the top of DET11

print(Euclid)
print(Euclid.catalog)




dr = np.sqrt((Euclid.catalog['X_WORLD']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (Euclid.catalog['Y_WORLD']-d0)**2)*3600.
id = Euclid.catalog['NUMBER'][np.argmin(dr)]
ind = np.argmin(dr)
obj_mag = Euclid.catalog['MAG_AUTO'][np.argmin(dr)]
print('ID:%d, IND:%d, mag=%.2f, dr=%.2f"' % (id, ind, obj_mag, np.min(dr)))






#beams = OrderedDict()


ix = Euclid.catalog['id'] == id
x0, y0 = Euclid.catalog['x_flt'][ix][0], Euclid.catalog['y_flt'][ix][0]
print(Euclid.direct.instrument, x0, y0)
#print(Euclid.wcs.pscale)
#dim = 18*0.135/sim.flt_wcs.pscale 
#beam = grizli.model.BeamCutout(id=id, x=x0, y=y0, 
#                               cutout_dimensions=np.cast[int]((dim, dim)), 
#                               conf=sim.conf, GrismFLT=sim)

print(Euclid.object_dispersers[id])


#sys.exit()



T0 = time.time()

os.chdir(os.path.join(HOME_PATH, root, 'Extractions'))


group_name = root + "_" + det

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

is_cgs, spectrum_1d, b = Euclid.object_dispersers[id]
cutout = grizli.model.BeamCutout(Euclid, b['A'], min_sens=0,) # min_mask=0) 

cutout.beam.compute_model()  
cutout.contam = cutout.beam.cutout_from_full_image(Euclid.model)
if id in Euclid.object_dispersers:
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



T1 = time.time()

print()
print("Finished in %.1f seconds" % (T1-T0))
