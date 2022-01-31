# Roman_notebooks


## Instructions for installing with Conda

- Fetch the grizli repo which includes the specialized Roman branch

`git clone https://github.com/gwalth/grizli.git`

- Checkout the Roman branch

`cd grizli`

`git checkout roman_sims_v1_gwalth`

- Generate a conda environment named "grizli_1.3.2"

`conda env create -f environment.yml -n grizli_1.3.2`

## Fixes issue:
### ImportError: cannot import name 'MalformedPolygonError'

`conda install tweakwcs=0.7.1`

- Activate the environment.

`source activate grizli_1.3.2`

- Compile and install the grizli module.

`python setup.py install`

Which will install in your Python distribution site-packages directory

If you foresee modifying the code regularly, use the develop option instead.  This
enables easier modification of the code.

`python setup.py develop`

## Setup the directory structure for the Jupyter notebook

Set BASE environmental variable to where you will put your Roman directory. 
`export BASE="/Users/gwalth/Roman"`


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
mkdir ${BASE}/my_roman_sims
mkdir ${BASE}/my_roman_sims/Prep
mkdir ${BASE}/my_roman_sims/Extractions
```

## Download Files
### Direct image
https://drive.google.com/file/d/1hH9NYz5FQ51rLKbifFB8IGeTEqdq0rP7/view?usp=sharing
### Slitless spectra
https://drive.google.com/file/d/1vlNAcJaOd8S27QkNjW6kWWyd_7_DmIkm/view?usp=sharing
### Primer catalog
https://drive.google.com/file/d/10RkOCTSFrgGsxY0wrsimuiv53qMC8rTa/view?usp=sharing
### CONF files
https://drive.google.com/file/d/1R2UWGL1M9Q3HrdZXVsH1Rwx-k4PkaNXi/view?usp=sharing
### SE files
https://drive.google.com/file/d/1BHPsQtHg_r2glOF6_LitsVrOCidDaC4E/view?usp=sharing


## Copy the necessary files

 Copy configurations files for the Roman grism
```
cd ${BASE}
tar -xvf ~/Downloads/Roman_grizli_conf_glw_v1.tar
```
Copy the SExtractor parameter files for Roman
```
cd ${BASE}/my_roman_sims
tar -xvf ~/Downloads/Roman_prep_glw_v1.tar
```

Copy direct image, slitless spectra and primer catalog to your working directory
```
cd ${BASE}/my_roman_sims/Prep
cp ~/Downloads/Roman_ATLAS_1deg_*_direct.fits .
cp ~/Downloads/Roman_ATLAS_1deg_*_slitless.fits .
cp ~/Downloads/ATLAS_1deg_primer_*.cat .
```

## Running the notebook



```
git clone https://github.com/gwalth/Roman_notebooks.git
```


