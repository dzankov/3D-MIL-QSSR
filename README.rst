3D-MIL-QSSR
============
Source code for building 3D multi-conformer models for predicting catalysts enantioselectivity

Installation
============

Install miniconda

.. code-block:: bash

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh

Create new environment with poetry

.. code-block:: bash

    $ conda create -n exp -c conda-forge "poetry=1.3.2" "python=3.10" -y
    $ conda activate exp

Install source code

.. code-block:: bash

    $ git clone https://github.com/dzankov/3D-MIL-QSAR

Install required packages

.. code-block:: bash

    $ cd 3D-MIL-QSSR
    $ poetry install --with cpu
    $ conda install -c conda-forge openbabel

Usage
============

* prepare the configuration file (see `config.yaml <https://github.com/dzankov/3D-MIL-QSSR/blob/main/config.yaml>`_)


* run model building:

.. code-block:: bash

    $ miqssr_build_model --config config.yaml

Graphical User Interface (GUI) 
============
https://chematlas.chimie.unistra.fr/Predictor/qscer.php





