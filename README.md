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

Data Formatting
=================
This package requires all input data to be in .csv format. Flux tower observation data should be formatted as the 'template_example.csv' files in data/MOD16 and data/MOD17 with columns named 'name', 'date', 'lat', 'lon', and 'target' correspoinding to flux tower names, date formatted as YYYY-MM-DD, tower latitude, tower longitude, and tower observation values (ET in the case of MOD16 and GPP in the case of MOD17), respectively.

An input file similar to 'template_example.csv' is the only thing necessary to run and calibrate the models. However, if you want to calibrate across different groups (e.g. plant functional types), a file formatted as 'group_template.csv' is required that assigns a 'group' column to each unique tower site.

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
