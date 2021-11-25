## VaMPy - Vascular Modeling Pypeline
[![Build Status](https://travis-ci.com/KVSlab/Aneurysm_workflow.svg?token=qbve9tcy6am6sUJksBcu&branch=master)](https://travis-ci.com/KVSlab/Aneurysm_workflow)
[![codecov](https://codecov.io/gh/KVSlab/Aneurysm_workflow/branch/master/graph/badge.svg?token=M2NMX6HOSZ)](https://codecov.io/gh/KVSlab/Aneurysm_workflow)

<p align="center">
    <img src=test/processed_model.png width="830 height="370" alt="Output pre processing"/>
</p>
<p align="center">
    Meshed and processed aneurysm model. Volumetric mesh (left), inlet flow rate, outlet flow split, and probes (middle), and temporal wall shear stress gradient (right).
</p>

Description
-----------
The Vascular Modeling Pypeline (VaMPy) is a collection of scripts used to prepare, run, and analyze cardiac and atrial problems.  This includes pre-processing scripts for meshing and probe sampling, a [Oasis](https://github.com/mikaem/Oasis) problem file for simulating flow in the [internal carotid artery](https://en.wikipedia.org/wiki/Internal_carotid_artery), and a variety of post-processing scripts for computing WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters in patient-specific geometries. 


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

VaMPy is Copyright (2018-2021) by the authors.

Documentation
-------------
For detailed installation notes and an introduction to VaMPy, please refer to the [documentation](https://vampy.readthedocs.io/en/latest/).

Installation
------------
For reference, VaMPy requires the following dependencies: VTK > 8.1, Numpy <= 1.13, SciPy > 1.0.0, VMTK 1.4, ITK, Paramiko, and FEniCS. 
If you are on Windows, macOS or Linux you can install all the general dependencies through anaconda.
First install Anaconda or Miniconda (preferably the Python 3.6 version).
Then create two environments, one for `vmtk/vtk` and one for `fenics` by executing the following in a terminal

    conda create -n vmtk -c vmtk python=3.6 itk vtk vmtk paramiko
    conda create -n fenics -c conda-forge fenics

You might run into a problem with VMTK 1.4 if using Python 3, or if you are a Linux user, and have therefore provided a set of temporary fixes for these known issues [here](https://vampy.readthedocs.io/en/latest/installation.html#known-issues).

Next, you need to install [`Oasis`](https://github.com/mikaem/Oasis). You can do so with the following commands:

    conda activate fenics
    git clone https://github.com/mikaem/Oasis
    cd Oasis
    pip install . && pip install cppimport # add "--user" if you are on a cluster, or "-e" if you are changing the Oasis source code

Finally, you are ready to clone and use the `VaMPy` repository:

    git clone https://github.com/KVSLab/VaMPy.git
    cd VaMPy

Issues
------
Please report bugs and other issues through the issue tracker at:

https://github.com/KVSlab/VaMPy/issues
