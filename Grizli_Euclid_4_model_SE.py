import sys
import os
import time 
import yaml
import grizli
import grizli.model
import pickle

from grizli_functions import euclid_det
####################################
# paramters
####################################
# allow simulation of objects at the edges
#pad=200 # pixels
pad = 800 # I think this may be optimal given the spectra size (only seems to add 10-20 spectra)
mag_limit = 30 
#yaml_file = "config.yaml"
yaml_file = sys.argv[1]
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

all_final_slitless = yaml_dict["all_final_slitless"]
all_ref_files = yaml_dict["all_ref_files"]
all_seg = yaml_dict["all_seg"]

T0 = time.time()

os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

with open('all_phot_matched_clean.pickle', 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    all_phot_matched_clean = pickle.load(f)

#all_euclid = []

all_det = euclid_det()
#print(all_det)
#sys.exit()

for i in range(len(all_final_slitless)):
    
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

    #all_euclid.append(Euclid)
    ## ~ 5GB file
    with open('Euclid_%s_GrismFLT.pickle' % (all_det[i]), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(Euclid, f, pickle.HIGHEST_PROTOCOL)
    
    t1 = time.time()

    del Euclid
    
    print("Detector finished in %.1f seconds" % (t1-t0))
    
    
T1 = time.time()

print()
print("Finished in %.1f seconds" % (T1-T0))


#print(all_euclid)



#t0 = time.time()

#print(all_euclid)



## ~ 5GB file
#with open('Euclid_GrismFLT.pickle', 'wb') as f:
#    # Pickle the 'data' dictionary using the highest protocol available.
#    pickle.dump(all_euclid, f, pickle.HIGHEST_PROTOCOL)
    
#t1 = time.time()

#print()
#print("Finished in %.1f seconds" % (t1-t0))


