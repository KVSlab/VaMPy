## Aneurysm workflow
[![Build Status](https://travis-ci.com/KVSlab/Aneurysm_workflow.svg?token=qbve9tcy6am6sUJksBcu&branch=master)](https://travis-ci.com/KVSlab/Aneurysm_workflow)
[![codecov](https://codecov.io/gh/KVSlab/Aneurysm_workflow/branch/master/graph/badge.svg?token=M2NMX6HOSZ)](https://codecov.io/gh/KVSlab/Aneurysm_workflow)

<p align="center">
    <img src=test/processed_model.png width="640 height="280" alt="Output pre processing"/>
</p>
<p align="center">
    Meshed aneurysm model showing inlet flow rate, outlet flow split, and probes.
</p>

Description
-----------
Aneurysm workflow is a collection of scripts to prepare, run, and analyze cardiac and atrial problems. This includes scripts for a variety of post-processing; WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters. 

The goal of the workflow is to provide research groups, and other individuals, with a set of tools for both pre- and post-processing patient-specific geometries.
Additionally, the workflow provides problem files for a set of vascular problems, which are readily available for simulation using the open-source CFD solver [Oasis](https://github.com/mikaem/Oasis).

Authors
-------
Aneurysm workflow has been developed by

* Aslak Wigdahl Bergersen
* Christophe Chnafa
* Henrik A. Kjeldsberg

Licence
-------
Aneurysm workflow is licensed under the GNU GPL, version 3 or (at your option) any
later version.

Aneurysm workflow is Copyright (2018-2021) by the authors.

Documentation
-------------
For detailed installation notes and an introduction to Aneurysm workflow, please refer to the [documentation](https://aneurysmworkflow.readthedocs.io/en/latest/).


Installation
------------
For reference, Aneurysm workflow requires the following dependencies: VTK > 8.1, Numpy <= 1.13, SciPy > 1.0.0, VMTK 1.4, ITK, Paramiko, and FEniCS. 
If you are on Windows, macOS or Linux you can install all the general dependencies through anaconda.
First install Anaconda or Miniconda (preferably the Python 3.6 version).
Then create two environments, one for `vmtk/vtk` and one for `fenics` by executing the following in a terminal

    conda create -n vtk -c vmtk python=3.6 itk vtk vmtk paramiko
    conda create -n fenics -c conda-forge fenics

You might run into a problem with VMTK 1.4 if using Python 3, or if you are a Linux user. Therefore, we have provided a set of temporary fixes for these known issues  [here](https://aneurysmworkflow.readthedocs.io/en/latest/installation.html#known-issues).

Now, you need to install [`Oasis`](https://github.com/mikaem/Oasis) and [`fenicstools`](https://github.com/mikaem/fenicstools/). You can do so with the following commands:

    conda activate fenics
    git clone https://github.com/mikaem/Oasis
    cd Oasis
    pip install .  # add "--user" if you are on a cluster, or "-e" if you are changing the Oasis source code

Similarly, `fenicstools` can be installed as follows:

    git clone https://github.com/mikaem/fenicstools.git
    cd fenicstools
    pip install . 
    pip install cppimport
    
If you have a cleaner install instruction please edit the above. 

Finally, you are ready to clone the `Aneurysm_workflow` repository:

    git clone https://github.com/KVSLab/Aneurysm_workflow.git
    cd Aneurysm_workflow


## Usage
First, use the automatedPreProcessing to create a mesh, boundary conditions, and probes for sampling. 

```
conda deactivate && conda activate vtk
python automatedPreProcessing/automatedPreProcessing.py -m diameter -i test/Case_test_artery/artery.vtp --aneurysm False -c 1.3
```

Then run a CFD simulation for two cycles with 10 000 time steps per cycle and default parameters with Oasis:
```
conda deactivate && conda activate fenics
oasis NSfracStep problem=Artery mesh_path=test/Case_test_artery/artery.xml.gz
```

Finally, you can create the WSS from the CFD simulation:
```
python postprocessing/compute_wss.py --case path_to_results/data/[run_number]/Solutions
```

You can also compute flow related metrics using `compute_flow_metrics.py`, but you would need to adapt how the files are read in to match with `compute_wss.py`.
To visualize velocity and pressure at the probes created by `Artery.py`, you can run the `visualize_probes.py` script, which has an additional dependency to [`Matplotlib`](https://github.com/matplotlib/matplotlib).
