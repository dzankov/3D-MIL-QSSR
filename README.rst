# 3D-MIL-QSSR

Source code for building 3D multi-conformer models for predicting catalysts enantioselectivity

# Installation

### Install miniconda

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh

### Create new environment with poetry

    $ conda create -n exp -c conda-forge "poetry=1.3.2" "python=3.10" -y
    $ conda activate exp

### Install source code

    $ git clone https://github.com/dzankov/3D-MIL-QSAR

### Install required packages

    $ cd 3D-MIL-QSSR
    $ poetry install --with cpu
    $ conda install -c conda-forge openbabel

# Usage

* Prepare the configuration file (see [config.yaml](config.yaml))


* Run model building:


    $ miqssr_build_model --config config.yaml

# Graphical User Interface (GUI) 
https://chematlas.chimie.unistra.fr/Predictor/qscer.php





