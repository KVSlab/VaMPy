.. title:: Installation

============
Installation
============
The Vascular Modeling Pypeline (VaMPy) is a collection of scripts to run an aneurysm problem with `Oasis <https://github.com/mikaem/Oasis>`_. There are also scripts for a variety of post-processing; WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters. The project is accessible through
`GitHub <https://github.com/KVSlab/VaMPy>`_.


Dependencies
============
The general dependencies of VaMPy are

* VMTK 1.4.0
* VTK 8.1.0
* Numpy <= 1.13
* SciPy 1.1.0
* Python (2.7 or >=3.5)


Basic Installation
==================
You can choose how to install the dependencies, but the fastest way to get started is to install the dependencies through Anaconda.
First, install Anaconda or Miniconda (preferably the Python 3.6 version) on your computer.
Then create two environments, one for `VMTK <http://www.vmtk.org/>`_ and one for `FEniCS <https://fenicsproject.org/>`_ by executing the following in a terminal window::

    conda create -n vmtk -c vmtk python=3.6 itk vtk vmtk paramiko
    conda create -n fenics -c conda-forge fenics

You can then activate your environment by running ``source activate [ENVIRONMENT NAME]``.

The next step is to install `Oasis <https://github.com/mikaem/Oasis>`_ and `fenicstools <https://github.com/mikaem/fenicstools/>`_.
You can do so with the following commands::

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

Now, all that is left is to clone the `VaMPy` repository::

    git clone https://github.com/KVSLab/VaMPy.git
    cd VaMPy

Now you are all set, and can start using the Vascular Modeling Pypeline.

Known issues
============

.. WARNING:: The `VMTK` version 1.4, the one currently distributed with Anaconda, has a Python3 bug in `vmtkcenterlines`, `vmtksurfacecurvature`, and `vmtkmeshwriter`. As a workaround you have to change these files. To find out where they are located you can use the ``which`` command  while in the ``vmtk`` environment. For `vmtkcenterlines` you can use the following command::
  
    which vmtkcenterlines

  Now copy the path up until ``vmtk`` and add ``lib/python3.6/site-packages/vmtk/vmtkcenterlines.py``.
  Please change the path separation symbol to match your operating system and change ``python3.6`` to the python version you are using. If you are using Miniconda, replace `anaconda3` with `miniconda3`.
  Using this path you can run the two following lines::

    sed -i -e 's/len(self.SourcePoints)\/3/len\(self.SourcePoints\)\/\/3/g' /Users/[Name]/anaconda3/envs/[your_environment]/lib/python3.6/site-packages/vmtk/vmtkcenterlines.py
    sed -i -e 's/len(self.TargetPoints)\/3/len\(self.TargetPoints\)\/\/3/g' /Users/[Name]/anaconda3/envs/[your_environment]/lib/python3.6/site-packages/vmtk/vmtkcenterlines.py

  Similarly, for `vmtksurfacecurvature.py`, run the following command::

    sed -i -e 's/(len(values) - 1)\/2/\(len\(values\) - 1\)\/\/2/g' /Users/[Name]/anaconda3/envs/[your_environment]/lib/python3.6/site-packages/vmtk/vmtksurfacecurvature.py

  Finally, to fix the issues in `vmtkmeshwriter.py`, change line 263 to::

    file = open(self.OutputFileName, 'rb')

  and line 267 to::

    gzfile = gzip.open(self.OutputFileName, 'wb')

  Please note that these changes are fixed in the development version of `VMTK`, but a new version has not been released in a while.


.. WARNING:: Some Linux users may experience the following Python compatibility issue::

    ModuleNotFoundError: No module named 'vtkRenderingOpenGL2Python'

  To fix this issue, a temporary solution is the install the ``llvm`` library directly in the virtual environment, using the following commands::

    conda config --set restore_free_channel true
    conda install llvm=3.3

.. WARNING:: Some Linux users may experience the following issue::

    ERROR: In ../Rendering/OpenGL2/vtkOpenGLRenderWindow.cxx, line 797

  To fix this issue, a temporary solution is to install `VTK` version `8.1.0` directly from the Anaconda channel. Assuming you have already tried to install `VTK` and `VMTK` and have the ``vmtk`` channel active, proceed with the following instructions.

  Remove `VMTK` 1.4 and `VTK` 8.1 that is installed from the `VMTK` Anaconda channel::

    conda uninstall vmtk vtk

  Install `VTK` 8.1.0 from the official Anaconda channel::

    conda install -c anaconda vtk

  Finally, install `VMTK` again::

    conda install -c vmtk vmtk

