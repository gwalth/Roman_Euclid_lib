import time
import yaml
import os, sys
import logging 

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
yaml_file = sys.argv[2]
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

failures = []

with open('Euclid_%s_GrismFLT.pickle' % (det), 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    Euclid = pickle.load(f)
    
t1 = time.time()

print()
print("Finished in %.1f seconds" % (t1-t0))


euclid_all,euclid_magcut,euclid_extract = check_sims2(Euclid, mag_limit)
Nl = [id for id in euclid_extract['NUMBER']]
print(det,"number_of_sources =",len(Nl))

Nl = Nl[:20]

T0 = time.time()

os.chdir(os.path.join(HOME_PATH, root, 'Extractions'))



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

print("Fitting redshifts for %s" % (det))

group_name = root + "_" + det

t0 = time.time()

if not os.path.exists(det):
    os.mkdir(det)    
os.chdir(det)

for j,id in enumerate(Nl):

    print("id =",id)
    print("%i of %i" % (j+1,len(Nl)))

    #beams = OrderedDict()

    try:

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

        logging.info("Successful fit id = %i" % (id))

    except:
        print("ERROR fitting redshift!!! id = %i" % (id))

        #logging.debug('This message should go to the log file')
        #logging.info('So should this')
        #logging.warning('And this, too')
        logging.exception("Problem fitting id = %i" % (id))

        failures.append(id)

        continue


os.chdir("..")

t1 = time.time()
print()
print("Finished %s in %.1f seconds" % (det,t1-t0))
print()
    
with open('failures.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(failures, f, pickle.HIGHEST_PROTOCOL)
    
T1 = time.time()

print()
print("Finished in %.1f seconds" % (T1-T0))
