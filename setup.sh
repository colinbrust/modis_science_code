#!/bin/bash
conda create -n $1 python=3
source activate $1
conda install -c conda-forge earthengine-api
conda install -c anaconda pandas
conda install -c anaconda scikit-learn
pip install spotpy
