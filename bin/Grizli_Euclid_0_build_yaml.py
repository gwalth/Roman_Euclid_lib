#!/usr/bin/env python

import argparse, os, yaml
import shutil


#################################################################################
# usage:
#
#  export EUCLID_LIB=/Users/gwalth/Euclid/grizli/Roman_Euclid_lib
#
#  export EUCLID_SIM=/Users/gwalth/Euclid/grizli/sims/FSpatch_mod3_16183
#
# cd $EUCLID_SIM
#
# python $EUCLID_LIB/Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/Euclid/grizli/sims
#          --dir_root FSpatch_mod3_16183 --slitless slitless_input --catalog catalog_input
# python $EUCLID_LIB/Grizli_Euclid_1_create_ref_images.py config.yaml
# python $EUCLID_LIB/Grizli_Euclid_2_run_SE.py config.yaml
# python $EUCLID_LIB/Grizli_Euclid_3_prep_SE.py config.yaml
# python $EUCLID_LIB/Grizli_Euclid_4_model_SE.py config.yaml
# python $EUCLID_LIB/Grizli_Euclid_5_fit_redshifts.py DET11 config.yaml
#
#################################################################################
#
#  ls ../Raw/EuclidSIMS/FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/*/*.fits > slitless.cat
#  ls ../Raw/EuclidSIMS/FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/FSpatch_mod3_16183_TAcalib_newGRID_V1.fits > catalog.cat
#  ls -d ../Raw/EuclidSIMS/FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/*/Input_Thumbnails > thumbnails.cat
#  Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/data/Roman/grizli/sims/Euclid
#      --dir_root FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17 --slitless slitless.cat --catalog catalog.cat
#      --thumbnails thumbnails.cat
#      --create
#
#

yaml_file = f"config.yaml"

parser = argparse.ArgumentParser(description='Create yaml config file for Euclid grizli')

parser.add_argument('--home_path', metavar='path',type=str, nargs='?', required=True,
                    help='', )

parser.add_argument('--dir_root', metavar='root',type=str, nargs='?', required=True,
                    help='')

parser.add_argument('--slitless', metavar='root',type=str, nargs='?', required=True,
                    help='')

parser.add_argument('--thumbnails', metavar='root',type=str, nargs='?', required=True,
                    help='')

parser.add_argument('--catalog', metavar='root',type=str, nargs='?', required=True,
                    help='')

parser.add_argument('--create', action='store_true', help='')
#parser.add_argument('--zodi', metavar='root',type=str, nargs='?', default="",
#                    help='')


args = parser.parse_args()
#print args

HOME_PATH = args.home_path
root = args.dir_root
slitless_input = args.slitless
#zodi_input = args.zodi
catalog_input = args.catalog
thumbnail_input = args.thumbnails
create = args.create

#print(slitless_input)



slitless_files_full = [l.strip() for l in open(slitless_input,'r').readlines() if l[0] != "#"]
catalog_files_full = [l.strip() for l in open(catalog_input,'r').readlines() if l[0] != "#"]
thumbnail_files_full = [l.strip() for l in open(thumbnail_input,'r').readlines() if l[0] != "#"]
#print(thumbnail_files_full)

slitless_files = [os.path.basename(sf) for sf in slitless_files_full]
catalog_files = [os.path.basename(cf) for cf in catalog_files_full]
#thumbnail_files = [os.path.basename(tf) for tf in thumbnail_files_full]
thumbnail_files = ["_".join([tf.split("/")[-1], tf.split("/")[-2]]) for tf in thumbnail_files_full]
#print(thumbnail_files)

#sys.exit()

#if zodi_input:
#    zodi_files = [l.strip() for l in open(zodi_input,'r').readlines() if l[0] != "#"]
#else: zodi_files = []

if create:

    print("Building Grizli directory structure")

    if not os.path.exists("Prep"):
        os.mkdir("Prep")
    
    if not os.path.exists("Extractions"):
        os.mkdir("Extractions")


    
    
    for src in slitless_files_full:
        dst = os.path.join("Prep", os.path.basename(src))
        #print(os.path.abspath(src))
        #print(os.path.relpath(src, start="Prep"))
        new_src = os.path.relpath(src, start="Prep")
        try: 
            os.symlink(new_src, dst)
        except FileExistsError:
            os.remove(dst)
            os.symlink(new_src, dst)
    
    for src in catalog_files_full:
        dst = os.path.join("Prep", os.path.basename(src))
        new_src = os.path.relpath(src, start="Prep")
        try: 
            os.symlink(new_src, dst)
        except FileExistsError:
            os.remove(dst)
            os.symlink(new_src, dst)

    for i, (src_num, src) in enumerate(zip(thumbnail_files, thumbnail_files_full)):
        #print(src_num)
        dst = os.path.join("Prep", src_num)
        new_src = os.path.relpath(src, start="Prep")
        try: 
            os.symlink(new_src, dst)
            if not i: os.symlink(new_src, "Prep/Input_Thumbnails")  # Until I combine all of the thumbnails
        except FileExistsError:
            os.remove(dst)
            os.symlink(new_src, dst)
            if not i: os.symlink(new_src, "Prep/Input_Thumbnails")  # Until I combine all of the thumbnails
    
    
    src = os.path.join(os.environ["ROMAN_EUCLID_LIB"], "SExtractor")
    destination = shutil.copytree(src, "Prep", dirs_exist_ok=True) 
    
    src = os.path.join(os.environ["ROMAN_EUCLID_LIB"], "utils")
    destination = shutil.copytree(src, "Prep", dirs_exist_ok=True) 



yaml_dict = {
    'HOME_PATH': HOME_PATH,
    'root': root,
    'slitless_files': slitless_files,
    'thumbnail_files': thumbnail_files,
#    'zodi_files': zodi_files,
    'catalog_files': catalog_files,
}

YAML_PATH = os.path.join(HOME_PATH, root)
os.chdir(YAML_PATH)
with open(yaml_file, 'w',) as f:
    yaml.dump(yaml_dict, f, sort_keys=False)
