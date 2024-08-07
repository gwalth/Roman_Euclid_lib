# Roman/Euclid Library

The Roman/Euclid library contains python packages wrapped around Grizli to help reduce and analyze
Roman/Euclid simulations and eventually data.  


## (1) Install Conda and Grizli

- Fetch the Grizli repo which includes the specialized Euclid branch

`git clone https://github.com/gwalth/grizli.git`

- Checkout the Euclid branch

`cd grizli`

`git checkout euclid_sims_v1_gwalth`

- Generate a conda environment named "grizli39"

`conda create -n grizli39 python=3.9 pip numpy cython`

- Activate the environment.

`conda activate grizli39`

- Compile and install the grizli module.

```
cd grizli
pip install --editable . -r requirements.txt
```

- Install eazy-py
```
git clone --recurse-submodules https://github.com/gbrammer/eazy-py.git
cd eazy-py
pip install . -r requirements.txt
```


## (2) Setup Grizli

Set GRIZLI_ROOT environmental variable to where you will put your Euclid directory. 
`export GRIZLI_ROOT="/Users/gwalth/Euclid"`


```
export GRIZLI="${GRIZLI_ROOT}/grizli"
export iref="${GRIZLI}/iref"
export jref="${GRIZLI}/jref"
```

```
mkdir ${GRIZLI_ROOT}
mkdir ${GRIZLI}
mkdir ${GRIZLI}/CONF
mkdir ${GRIZLI}/templates
mkdir ${GRIZLI}/iref
mkdir ${GRIZLI}/jref
```

Make the directories where your simulated files will live.  Prep and Extractions are hard coded in the notebook.
```
mkdir ${GRIZLI_ROOT}/my_euclid_sims
mkdir ${GRIZLI_ROOT}/my_euclid_sims/Prep
mkdir ${GRIZLI_ROOT}/my_euclid_sims/Extractions
```


## (3) Install Roman_Euclid_lib
cd ${GRIZLI_ROOT}
```
git clone https://github.com/gwalth/Roman_Euclid_lib.git
cd Roman_Euclid_lib
source libenv_setup.sh
```

## Download Files


## Copy the necessary files

 Copy configurations files for the Euclid grisms
```
cd ${GRIZLI_ROOT}/grizli/CONF
mkdir Euclid
cd Euclid
cp /local/SIRsim/SIM_10_18_22_singleframe/frame_1/CONF* .
```

```
cd ${GRIZLI_ROOT}/my_roman_sims/Prep
tar -xvf ~/Downloads/Euclid_prep_glw_v1.tar
```

Copy direct images, slitless spectra and primer catalog to your working directory
```
cd ${GRIZLI_ROOT}/my_roman_sims/Prep
cp -r /local/SIRsim/SIM_10_18_22_singleframe/SIM_intermediate_files/midfiles_frame_1/Input_Thumbnails .
cp /local/SIRsim/SIM_10_18_22_singleframe/data/EUC_SIM_NISRGS000-0-1_20220913T230154.141Z_TEST_SC8_NIS_S1.fits .
cp ~/Downloads/CATALOG_WP9_INTERMEDIATE_RUN_v2_NewCoord_mod.fits .
```

## Running the notebook

```
source ~/conda_setup.sh                 # conda
conda activate grizli39 
cd /Users/gwalth/data/Roman/grizli
source grizli_env.sh                    # grizli
cd Roman_Euclid_lib
source libenv_setup.sh                  # Roman_Euclid_lib
```

source ~/sandbox.sh  # sex and ds9


# modify Grizli_Euclid_1_create_ref_images.py to files and catalogs
```

export EUCLID_SIM=/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11
export EUCLID_SIM=/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v1
export EUCLID_SIM=/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v2
export EUCLID_SIM=/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17


cd $EUCLID_SIM
#Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/data/Roman/grizli/sims/Euclid --dir_root FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v1 --slitless slitless_input --catalog catalog_input
#Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/data/Roman/grizli/sims/Euclid --dir_root FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v2 --slitless slitless_input --catalog catalog_input
Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/data/Roman/grizli/sims/Euclid --dir_root FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17 --slitless slitless.cat --catalog catalog.cat --thumbnails thumbnails.cat --create
Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/data/Roman/grizli/sims/Euclid --dir_root FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17 --slitless slitless.cat --catalog catalog.cat --thumbnails thumbnails.cat 

Grizli_Euclid_1_create_ref_images.py config.yaml
Grizli_Euclid_2_run_SE.py config.yaml
Grizli_Euclid_3_prep_SE.py config.yaml
Grizli_Euclid_4_model_SE.py config.yaml

sh test_fit_redshifts.sh
sh test_fit_redshifts_v2.sh

cd $EUCLID_SIM/Extractions
ls *stack.fits | wc   # should be 16
ls *beams.fits | wc   # should be 16
python $ROMAN_EUCLID_LIB/tests/test_inspect2d.py "*stack.fits"
python $ROMAN_EUCLID_LIB/tests/test_inspect2d.py "*beams.fits"
```

```
python $ROMAN_EUCLID_LIB/tests/test_display1d.py FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11 DET11 160
```

