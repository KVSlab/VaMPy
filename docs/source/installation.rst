.. title:: Installation

============
Installation
============
The Vascular Modeling Pypeline (VaMPy) is a collection of scripts used to prepare, run, and analyze cardiac and atrial morphologies.  This includes pre-processing scripts for meshing and probe sampling, a `Oasis <https://github.com/mikaem/Oasis>`_ problem file for simulating flow in the `internal carotid artery <https://en.wikipedia.org/wiki/Internal_carotid_artery>`_, and a variety of post-processing scripts for computing WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters in patient-specific geometries. The project is accessible through
`GitHub <https://github.com/KVSlab/VaMPy>`_.


Dependencies
============
The general dependencies of VaMPy are

* VMTK 1.4
* morphMan
* FEniCS
* Paramiko

Basic Installation
==================
You can choose how to install the dependencies, but the fastest way to get started is to install the dependencies through Anaconda.
First, install Anaconda or Miniconda (preferably the Python 3.6 version) on your computer.
Then create two environments, one for `morphMan <https://github.com/KVSlab/morphMan>`_ and `VMTK <http://www.vmtk.org/>`_, and one for `FEniCS <https://fenicsproject.org/>`_ by executing the following in a terminal window::

    conda create -n morphman -c vmtk -c morphman morphman python=3.6 paramiko
    conda create -n fenics -c conda-forge fenics

You can then activate your environment by running ``source activate [ENVIRONMENT NAME]``.
Windows users may need to install FEniCS as described `here <https://fenicsproject.org/download/>`_.

The next step is to install `Oasis <https://github.com/mikaem/Oasis>`_.
You can do so with the following commands::

    conda activate fenics
    cd [path_to_your_installation_folder]
    git clone https://github.com/mikaem/Oasis
    cd Oasis
    pip install cppimport && pip install . 

Now, all that is left is to clone the `VaMPy` repository::

    git clone https://github.com/KVSLab/VaMPy.git
    cd VaMPy

Now you are all set, and can start using the Vascular Modeling Pypeline.

Known issues
============

.. WARNING:: The `morphMan` dependency `VMTK` (version 1.4) currently distributed with Anaconda has a Python3 bug in ``vmtkcenterlines.py``, ``vmtksurfacecurvature.py``, and ``vmtkmeshwriter.py``. As a workaround you have to change these files. We have provided the script ``apply_vmtk_hotfixes.py``, which will automatically edit the files for you, provided you enter your local username, Anaconda version and environment name. The script may be executed with the following command::

    python apply_vmtk_hotfixes.py

  Alternatively, you may edit the files manually. To find out where they are located you can use the ``which`` command  while in the ``morphman`` environment. For `vmtkcenterlines` you can use the following command::
  
    which vmtkcenterlines

  Now copy the path up until ``morphman`` and add ``lib/python3.6/site-packages/vmtk/vmtkcenterlines.py``.
  Please change the path separation symbol to match your operating system and change ``python3.6`` to the python version you are using. If you are using Miniconda, replace `anaconda3` with `miniconda3`.
  Using this path you can run the two following lines::

    sed -i -e 's/len(self.SourcePoints)\/3/len\(self.SourcePoints\)\/\/3/g' /Users/[Username]/anaconda3/envs/morphman/lib/python3.6/site-packages/vmtk/vmtkcenterlines.py
    sed -i -e 's/len(self.TargetPoints)\/3/len\(self.TargetPoints\)\/\/3/g' /Users/[Username]/anaconda3/envs/morphman/lib/python3.6/site-packages/vmtk/vmtkcenterlines.py

  Similarly, for `vmtksurfacecurvature.py`, run the following command::

    sed -i -e 's/(len(values) - 1)\/2/\(len\(values\) - 1\)\/\/2/g' /Users/[Username]/anaconda3/envs/morphman/lib/python3.6/site-packages/vmtk/vmtksurfacecurvature.py

  Finally, to fix the issue in `vmtkmeshwriter.py`, run the following command::

    sed -i -e -r "s/file = open\(self\.OutputFileName, ?\'r\'\)/file = open\(self\.OutputFileName, \'rb\'\)/g" /Users/[Username]/anaconda3/envs/morphman/lib/python3.6/site-packages/vmtk/vmtkmeshwriter.py

  Please note that these changes are fixed in the development version of `VMTK`, but a new version has not been released in a while.


.. WARNING:: Some Linux users may experience the following Python compatibility issue::

    ModuleNotFoundError: No module named 'vtkRenderingOpenGL2Python'

  To fix this issue, a temporary solution is the install the ``llvm`` library directly in the virtual environment, using the following commands::

    conda config --set restore_free_channel true
    conda install llvm=3.3

.. WARNING:: Some Linux users may experience the following issue::

    ERROR: In ../Rendering/OpenGL2/vtkOpenGLRenderWindow.cxx, line 797

  To fix this issue, a temporary solution is to install `VTK` version `8.1.0` directly from the Anaconda channel. Assuming you have already tried to install `VTK` and `VMTK` and have the ``morphman`` environment active, proceed with the following instructions.

  Remove `VMTK` 1.4 and `VTK` 8.1 that is installed from the `VMTK` Anaconda channel::

    conda uninstall vmtk vtk

  Install `VTK` 8.1.0 from the official Anaconda channel::

    conda install -c anaconda vtk=8.1.0

  Finally, install `VMTK` again::

    conda install -c vmtk vmtk

