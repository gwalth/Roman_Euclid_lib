# Roman/Euclid Library

The Roman/Euclid library contains python packages wrapped around grizli to help reduce and analyze
Roman/Euclid simulations and eventually data.  


## Instructions for installing with Conda

- Fetch the grizli repo which includes the specialized Euclid branch

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


## Setup the directory structure for the Jupyter notebook

Set BASE environmental variable to where you will put your Euclid directory. 
`export BASE="/Users/gwalth/Euclid"`


```
export GRIZLI="${BASE}/grizli"
export iref="${GRIZLI}/iref"
export jref="${GRIZLI}/jref"
```

```
mkdir ${BASE}
mkdir ${GRIZLI}
mkdir ${GRIZLI}/CONF
mkdir ${GRIZLI}/templates
mkdir ${GRIZLI}/iref
mkdir ${GRIZLI}/jref
```

Make the directories where your simulated files will live.  Prep and Extractions are hard coded in the notebook.
```
mkdir ${BASE}/my_euclid_sims
mkdir ${BASE}/my_euclid_sims/Prep
mkdir ${BASE}/my_euclid_sims/Extractions
```

## Download Files
### Source Extractor files
https://drive.google.com/file/d/1jQAypoi_YpD0BWY9SnlaP1G51dhX6Aay/view?usp=share_link


## Copy the necessary files

 Copy configurations files for the Euclid grisms
```
cd ${BASE}/grizli/CONF
mkdir Euclid
cd Euclid
cp /local/SIRsim/SIM_10_18_22_singleframe/frame_1/CONF* .
```

Copy the SExtractor parameter files for Euclid
```
cd ${BASE}/my_roman_sims/Prep
tar -xvf ~/Downloads/Euclid_prep_glw_v1.tar
```

Copy direct images, slitless spectra and primer catalog to your working directory
```
cd ${BASE}/my_roman_sims/Prep
cp -r /local/SIRsim/SIM_10_18_22_singleframe/SIM_intermediate_files/midfiles_frame_1/Input_Thumbnails .
cp /local/SIRsim/SIM_10_18_22_singleframe/data/EUC_SIM_NISRGS000-0-1_20220913T230154.141Z_TEST_SC8_NIS_S1.fits .
cp ~/Downloads/CATALOG_WP9_INTERMEDIATE_RUN_v2_NewCoord_mod.fits .
```

## Running the notebook



cd /Users/gwalth/data/Roman/grizli
```
git clone https://github.com/gwalth/Roman_Euclid_lib.git
```

```
source ~/conda_setup.sh
source ~/sandbox.sh
conda activate grizli39
cd /Users/gwalth/data/Roman/grizli
source grizli_env.sh
cd Roman_Euclid_lib
```


# modify Grizli_Euclid_1_create_ref_images.py to files and catalogs
```

export EUCLID_SIM=/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v1
export EUCLID_SIM=/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v2


cd $EUCLID_SIM
python $BASE/Roman_Euclid_lib/Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/data/Roman/grizli/sims/Euclid --dir_root FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v1 --slitless slitless_input --catalog catalog_input
python $BASE/Roman_Euclid_lib/Grizli_Euclid_0_build_yaml.py --home_path /Users/gwalth/data/Roman/grizli/sims/Euclid --dir_root FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_v2 --slitless slitless_input --catalog catalog_input

python $BASE/Roman_Euclid_lib/Grizli_Euclid_1_create_ref_images.py config.yaml
python $BASE/Roman_Euclid_lib/Grizli_Euclid_2_run_SE.py config.yaml
python $BASE/Roman_Euclid_lib/Grizli_Euclid_3_prep_SE.py config.yaml
python $BASE/Roman_Euclid_lib/Grizli_Euclid_4_model_SE.py config.yaml

cd $BASE/Roman_Euclid_lib/
sh test_fit_redshifts.sh

cd $EUCLID_SIM/Extractions
ls *stack.fits | wc   # should be 16
ls *beams.fits | wc   # should be 16
python $BASE/Roman_Euclid_lib/tests/test_inspect2d.py "*stack.fits"
python $BASE/Roman_Euclid_lib/tests/test_inspect2d.py "*beams.fits"
```



