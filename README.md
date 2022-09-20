## VaMPy - Vascular Modeling Pypeline
[![CircleCI](https://circleci.com/gh/KVSlab/VaMPy/tree/master.svg?style=svg)](https://circleci.com/gh/hkjeldsberg/VaMPy/tree/master)
[![codecov](https://codecov.io/gh/KVSlab/VaMPy/branch/master/graph/badge.svg?token=M2NMX6HOSZ)](https://codecov.io/gh/KVSlab/VaMPy)
[![Documentation Status](https://readthedocs.org/projects/vampy/badge/?version=latest)](https://vampy.readthedocs.io/en/latest/?badge=latest)

<p align="center">
    <img src=test/artery_pipeline.png width="830 height="370" alt="Output pre processing"/>
</p>
<p align="center">
    Pre-processed and simulated artery model. From left to right: A variable density volumetric mesh, zoomed in view of an outlet showing the four boundary layers, corresponding inlet flow rate, outlet flow split, and probes for velocity and pressure sampling. From the simulation results, we have shown the velocity field represented by vectors, and the time averaged wall shear stress (TAWSS) as one of the many hemodynamic indices computed by the post-processing scripts.
</p>

Description
-----------
The Vascular Modeling Pypeline (VaMPy) is a collection of fully automated scripts used to prepare, run, and analyze cardiac and atrial morphologies.  This includes pre-processing scripts for meshing and probe sampling, two [Oasis](https://github.com/mikaem/Oasis) problem files for simulating flow in the [internal carotid artery](https://en.wikipedia.org/wiki/Internal_carotid_artery) and the [left atrium](https://en.wikipedia.org/wiki/Atrium_(heart)), and a variety of post-processing scripts for computing WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters in patient-specific geometries. 

Authors
-------
VaMPy has been developed by

* Aslak Wigdahl Bergersen
* Christophe Chnafa
* Henrik A. Kjeldsberg

Licence
-------
VaMPy is licensed under the GNU GPL, version 3 or (at your option) any
later version.

VaMPy is Copyright (2018-2022) by the authors.

Documentation
-------------
For detailed installation notes and an introduction to VaMPy, please refer to the [documentation](https://vampy.readthedocs.io/en/latest/).

Installation
------------
For reference, VaMPy requires the following dependencies: morphMan, FEniCS, and Paramiko. 
If you are on Windows, macOS or Linux you can install all the general dependencies through [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
First install Anaconda or Miniconda on your computer, and add the `conda-forge` channel with:

    conda config --add channels conda-forge
    conda config --set channel_priority strict

Then create a conda environment for `morphMan` and `FEniCS` by executing the following command in a terminal

    conda create -n your_environment morphman fenics paramiko 

Next, you need to install [`Oasis`](https://github.com/mikaem/Oasis). You can do so with the following commands:

    git clone https://github.com/mikaem/Oasis
    pip install cppimport
    pip install --editable Oasis

Finally, you are ready to clone and use the Vascular Modeling Pypeline:

    git clone https://github.com/KVSLab/VaMPy.git
    cd VaMPy

Issues
------
Please report bugs and other issues through the issue tracker at:

https://github.com/KVSlab/VaMPy/issues
