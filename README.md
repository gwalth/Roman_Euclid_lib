# Roman_notebooks


## Instructions for installing with Conda

- Fetch the grizli repo which includes the specialized Roman branch

`git clone https://github.com/gwalth/grizli.git`

- Checkout the Roman branch

`cd grizli`

`git checkout roman_sims_v1_gwalth`

- Generate a conda environment named "grizli_1.3.2"

`conda env create -f environment.yml -n grizli_1.3.2`

- Activate the environment.

`source activate grizli_1.3.2`

- Compile and install the grizli module.

`python setup.py install`

Which will install in your Python distribution site-packages directory

If you foresee modifying the code regularly, use the develop option instead.  This
enables easier modification of the code.

`python setup.py develop`

