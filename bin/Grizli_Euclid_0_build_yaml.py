#!/usr/bin/env python

import argparse, os, yaml


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

yaml_file = f"config.yaml"

parser = argparse.ArgumentParser(description='Create yaml config file for Euclid grizli')

parser.add_argument('--home_path', metavar='path',type=str, nargs='?', required=True,
                    help='', )

parser.add_argument('--dir_root', metavar='root',type=str, nargs='?', required=True,
                    help='')

parser.add_argument('--slitless', metavar='root',type=str, nargs='?', required=True,
                    help='')

parser.add_argument('--catalog', metavar='root',type=str, nargs='?', required=True,
                    help='')

parser.add_argument('--zodi', metavar='root',type=str, nargs='?', default="",
                    help='')


args = parser.parse_args()
#print args

HOME_PATH = args.home_path
root = args.dir_root
slitless_input = args.slitless
zodi_input = args.zodi
catalog_input = args.catalog


slitless_files = [l.strip() for l in open(slitless_input,'r').readlines() if l[0] != "#"]
catalog_files = [l.strip() for l in open(catalog_input,'r').readlines() if l[0] != "#"]

if zodi_input:
    zodi_files = [l.strip() for l in open(zodi_input,'r').readlines() if l[0] != "#"]
else: zodi_files = []


yaml_dict = {
    'HOME_PATH': HOME_PATH,
    'root': root,
    'slitless_files': slitless_files,
    'zodi_files': zodi_files,
    'catalog_files': catalog_files,
}

YAML_PATH = os.path.join(HOME_PATH, root)
os.chdir(YAML_PATH)
with open(yaml_file, 'w',) as f:
    yaml.dump(yaml_dict, f, sort_keys=False)
