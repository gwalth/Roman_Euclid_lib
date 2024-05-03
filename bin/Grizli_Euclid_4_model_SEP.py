#!/usr/bin/env python

import glob, os, sys
import yaml
import time
import pickle

from grizli import model, multifit, grismconf
from grizli.pipeline import auto_script
from grizli.pipeline.auto_script import get_yml_parameters

yaml_file = sys.argv[1]
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

## Read the default parameters that can be edited and passed en-masse to `auto_script.go`
kwargs = get_yml_parameters()
print(list(kwargs.keys()))

#pad=200 # pixels
pad = 800 # I think this may be optimal given the spectra size (only seems to add 10-20 spectra)
prefix = "Euclid"

#os.chdir(os.path.join(HOME_PATH + root, 'Prep'))
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
         
if not os.path.exists('../Extractions'):
    os.mkdir('../Extractions')
    
files = glob.glob('*GrismFLT.fits')

if len(files) == 0:
    ### Grism contamination model

    # Which filter to use as direct image?  Will try in order of the list until a match is found.
    grism_prep_args = kwargs['grism_prep_args']
   
    files = glob.glob("*_flt.fits")
    files.sort()
    files = [files[0]]
    print(files)
    #files = glob.glob("*_slitless_final.fits")
    grism_prep_args['files'] = files
    
    # For now, turn off refining contamination model with polynomial fits
    #grism_prep_args['refine_niter'] = 2
    grism_prep_args['refine_niter'] = 0
    
    # Flat-flambda spectra
    grism_prep_args['init_coeffs'] = [1.0]
    
    grism_prep_args['mask_mosaic_edges'] = False
    
    # Fairly bright for speedup, these can be adjusted based on how deep the spectra/visits are
    grism_prep_args['refine_mag_limits'] = [17,24]
    grism_prep_args['prelim_mag_limit'] = 25
    
    grism_prep_args['gris_ref_filters'] = {'RED':['total']}
    
    print(grism_prep_args)
    print("GLW GLW GLW GLW HERE 1!")
    
    grp = auto_script.grism_prep(field_root=prefix, pad=pad, **grism_prep_args)

    print("GLW GLW GLW GLW HERE 2!")
    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(prefix), 
                            cpu_count=-1, sci_extn=1, pad=pad)
    
else:
    os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(prefix), 
                            cpu_count=-1, sci_extn=1, pad=pad) 
