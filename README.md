# MODIS Science Code

---

A set of classes and functions to extract Google Earth Engine data for MOD16 and MOD17 and calibrate model parameters. Examples of how to run key functions are located in the 'Examples' directory. 

Requirements
=================
- git
- conda
- Python 3
    - Google Earth Engine Python API
    - pandas
    - spotpy
    
The following install instructions will install the modis_science_code package and all additional requirements 
(with the exception of git and conda).

Install
=================
I would recommend using a conda environment for code in this project so that it doesn't affect your base python install.

If you are using a Unix-like OS, you can paste the following in the command line to create a conda environment with
 necessary packages:

	git clone https://github.com/colinbrust/modis_science_code.git
	cd modis_science_code
	chmod 777 setup.sh
	./setup.sh gee # NOTE: you can change 'gee' to whatever you want the conda environment to be named.
	source activate gee # If you speficied a name other than 'gee' change that here.
	python setup.py install

If you are using a Windows OS, running the following in the command line will create a conda environment with necessary 
packages:
    
    git clone https://github.com/colinbrust/modis_science_code.git
	cd modis_science_code
	conda env create -f environment.yml
	conda activate gee
	python setup.py install

If you haven't used the GEE Python API, you will have to run `earthengine authenticate` in the command line to 
authenticate your EarthEngine account after following the above steps.
