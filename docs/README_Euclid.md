# Euclid_notebooks


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



```
git clone https://github.com/gwalth/Euclid_notebooks.git
```

