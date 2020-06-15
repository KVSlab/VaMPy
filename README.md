# Aneurysm_workflow

This is a collection of scripts to run an aneurysm problem with Oasis. There are also scripts for a variety of post-processing; WSS-based matrics, more advanced turbulence metrics, and a variety of morphological parameters. The latter is implemented through automated neck plane detection, but are not adapted to the `Aneurysm_workflow` pipeline and are here merely for convenience.

## Authors
These scripts was written by
- Aslak Wigdahl Bergersen
- Christophe Chnafa
- Henrik A. Kjeldsberg

## Installation
You can choose how to install the dependencies, but the fastest way to get started is to first install anaconda or minicoda on your computer. Then create two environments, one for `vmtk/vtk` and one for `fenics` by executing the following in a terminal:
```
conda create -n vtk -c vmtk python=3.6 itk vtk vmtk paramiko
conda create -n fenics -c conda-forge fenics
```

You might run into a problem with vmtk (1.4) if using python 3. To fix this, please [follow these instructions](https://morphman.readthedocs.io/en/latest/installation.html#basic-installation) for fixing the centerline problem. For fixing mesh writing change line 263 to:
```
file = open(self.OutputFileName, 'rb')
````
and line 267 to
```
gzfile = gzip.open(self.OutputFileName, 'wb')
```
Please note that these changes are fixed in the development version of vmtk, but a new version has not been released in a long time.

Now, you need to install `Oasis`and `fenicstools`. You can do so with the following commands:
```
conda activate fenics
cd [path_to_your_installation_folder]
git clone https://github.com/mikaem/Oasis
cd Oasis
pip install .  # add "--user" if you are on a cluster, or "-e" if you are changing the Oasis source code
cd ..
pip install cppimport
git clone https://github.com/mikaem/fenicstools.git
cd fenicstools
pip install . 
```

## Usage
First, use the automatedPreProcessing to create a mesh, boundary conditions, and probes for sampling. 

```
conda deactivate
conda activate vtk
git clone https://github.com/KVSLab/Aneurysm_workflow.git
cd Aneurysm_workflow/automatedPreProcessing
python automatedPreProcessing.py -m diameter -i ../test/Case_test_71.vtp --aneurysm False -c 1.3
cd ..
```

Then run a CFD simulation for two cycles with 10 000 time steps per cycle and default parameters with Oasis:
```
conda deactivate
conda activate fenics
oasis NSfracStep problem=Artery mesh_path=test/Case_test_71.xml.gz
```

Finally, you can create the WSS from the CFD simulation:
