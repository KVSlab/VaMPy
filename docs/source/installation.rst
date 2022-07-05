.. title:: Installation

============
Installation
============
.. highlight:: console

The Vascular Modeling Pypeline (VaMPy) is a collection of scripts used to prepare, run, and analyze cardiac and atrial morphologies.  This includes pre-processing scripts for meshing and probe sampling, a `Oasis <https://github.com/mikaem/Oasis>`_ problem file for simulating flow in the `internal carotid artery <https://en.wikipedia.org/wiki/Internal_carotid_artery>`_, and a variety of post-processing scripts for computing WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters in patient-specific geometries. The project is accessible through
`GitHub <https://github.com/KVSlab/VaMPy>`_.


Dependencies
============
The main dependencies of VaMPy are

* morphMan
* Oasis (FEniCS)
* Paramiko
* Scipy

You can choose how to install the dependencies, but the fastest way to get started is to install the dependencies through `conda-forge`.
Start by installing `conda`, by downloading and installing `Anaconda <https://www.anaconda.com/products/distribution>`_ or `Miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ on your computer.

Installing morphMan and FEniCS
==============================
The `morphMan <https://github.com/KVSlab/morphMan>`_ software can be installed on all operative systems through Anaconda or `conda-forge`.
Similarly, the `FEniCS <https://fenicsproject.org/>`_ software is also available through Anaconda.
Note that Windows users may need to install FEniCS as described `here <https://fenicsproject.org/download/>`_.
To create a conda environment with both software's installed, start by adding the `conda-forge` channel with::

    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict

Once the `conda-forge` channel has been added, you can create a `conda` environment with `morphMan` and `FEniCS` installed with::

    $ conda create -n your_environment morphman fenics paramiko

.. note::
    Replace ``your_environment`` with the environment name.

This will also install the minor dependency to `paramiko`.
You can now activate your environment by running::

    $ conda activate your_environment

or::

    $ source activate your_environment

.. WARNING:: Some users may experience errors regarding compatibility if Anaconda already has been configured with certain channels. To resolve this issue you can remove conflicting channels using::

    $ conda config --remove channels [CHANNEL NAME]

  Alternatively, you set your Anaconda channel priority to *flexible*, with the following command::

    $ conda config --set channel_priority flexible


Installing Oasis
================
The next step is to download and install `Oasis <https://github.com/mikaem/Oasis>`_ on your machine.
Inside your conda environment run the following commands to clone and install Oasis, and the minor dependency to `cppimport`::

    $ git clone https://github.com/mikaem/Oasis
    $ pip install cppimport
    $ pip install --editable Oasis

The ``--editable`` flag installs the project in editable mode meaning that any changes to the original package will be reflected directly in your environment.

Downloading VaMPy
=================
All that is left is to clone the `VaMPy` repository::

    $ git clone https://github.com/KVSLab/VaMPy.git
    $ cd VaMPy

Now you are all set, and can start using the Vascular Modeling Pypeline.
