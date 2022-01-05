.. title:: Installation

============
Installation
============
The Vascular Modeling Pypeline (VaMPy) is a collection of scripts used to prepare, run, and analyze cardiac and atrial morphologies.  This includes pre-processing scripts for meshing and probe sampling, a `Oasis <https://github.com/mikaem/Oasis>`_ problem file for simulating flow in the `internal carotid artery <https://en.wikipedia.org/wiki/Internal_carotid_artery>`_, and a variety of post-processing scripts for computing WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters in patient-specific geometries. The project is accessible through
`GitHub <https://github.com/KVSlab/VaMPy>`_.


Dependencies
============
The general dependencies of VaMPy are

* VTK 8.1.0
* VMTK 1.4
* morphMan
* FEniCS
* Paramiko

You can choose how to install the dependencies, but the fastest way to get started is to install the dependencies through Anaconda.
You can chose installing either Anaconda or Miniconda on your computer, as described `here <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.

Installing morphMan on Windows or macOS
=======================================
The `morphMan <https://github.com/KVSlab/morphMan>`_ software can be installed on Windows and macOS by executing the following in a terminal window::

    conda create -n morphman -c vmtk -c morphman morphman python=3.6 paramiko

You can then activate your environment by running ``source activate morphman``.

Installing morphMan on Linux
============================
Linux users may experience various compatibility issues by following the instructions above.
Therefore, to install the `morphMan <https://github.com/KVSlab/morphMan>`_ software on Linux, execute the following commands::

    conda config --set restore_free_channel true
    conda create -n morphman -c anaconda vtk=8.1.0
    conda activate morphman
    conda install -c vmtk -c morphman morphman paramiko

You can then activate your environment by running ``source activate morphman``.

Installing FEniCS and Oasis
===========================
The next step is to create a separate environment for `FEniCS <https://fenicsproject.org/>`_, and installing `Oasis <https://github.com/mikaem/Oasis>`_.
You can do so with the following commands::

    conda create -n fenics -c conda-forge fenics
    conda activate fenics
    cd [path_to_your_installation_folder]
    git clone https://github.com/mikaem/Oasis
    cd Oasis
    pip install cppimport && pip install . 

Note that Windows users may need to install FEniCS as described `here <https://fenicsproject.org/download/>`_.

Downloading VaMPy
=================
All that is left is to clone the `VaMPy` repository::

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




