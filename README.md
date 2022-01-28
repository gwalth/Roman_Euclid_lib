# Roman_notebooks


## Instructions for installing with Conda

- Fetch the grizli repo which includes the specialized Roman branch

`git clone https://github.com/gwalth/grizli.git`

- Checkout the Roman branch

`cd grizli`

`git checkout roman_sims_v1_gwalth`

- Generate a conda environment named "grizli_1.3.2"

`conda env create -f environment.yml -n grizli_1.3.2`

# Fixes issue:
# ImportError: cannot import name 'MalformedPolygonError'

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

```
export BASE="/Users/gwalth/Roman"
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

```
mkdir ${BASE}/my_roman_sims
mkdir ${BASE}/my_roman_sims/Prep
mkdir ${BASE}/my_roman_sims/Extractions
```

# Copy the necessary files

```
cd ${BASE}
```

```
cp ~/data/Roman/grizli/grizli/CONF/Roman.G150*conf           grizli/CONF
cp ~/data/Roman/grizli/grizli/CONF/sens_0720_2020.fits       grizli/CONF
cp ~/data/Roman/grizli/grizli/CONF/Roman.G150.v1.6.sens.fits grizli/CONF
```

# Direct image
https://drive.google.com/file/d/10kvh6s4lqI9AVtMZYPuXTG8O-RSMFO9s/view?usp=drive_web
# Slitless spectra
https://drive.google.com/file/d/1ZHXJGAX0KGPf9pJH0JJM_add1ZNmFPtq/view?usp=drive_web

```
cp ~/Downloads/Roman_ptg01_WFI_G150_random_*  my_roman_sims/Prep
cp ~/data/Roman/grizli/my_roman_sims/Prep/ATLAS_1deg_subsample_primer_random.lis my_roman_sims/Prep
```

```
cp ~/data/Roman/grizli/my_roman_sims/Prep/Roman.param .
cp ~/data/Roman/grizli/my_roman_sims/Prep/Roman.sex .
cp ~/data/Roman/grizli/my_roman_sims/Prep/gauss_5.0_9x9.conv .
cp ~/data/Roman/grizli/my_roman_sims/Prep/default.nnw .
```


## Running the notebook



```
git clone https://github.com/gwalth/Roman_notebooks.git
```


