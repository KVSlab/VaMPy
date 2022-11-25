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


Downloading VaMPy
=================
All that is left is to clone the `VaMPy` repository::

    $ git clone https://github.com/KVSLab/VaMPy.git
    $ cd VaMPy

You can choose how to install the dependencies, but the fastest way to get started is to install the dependencies through `conda-forge`.

Conda
=====

Start by installing ``conda``, by downloading and installing `Anaconda <https://www.anaconda.com/products/distribution>`_ or `Miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ on your computer.

Using the ``environment.yml`` file in the root of the repository, you can call::

    $ conda env update --file environment.yml --name your_environment

You can now activate your environment by running::

    $ conda activate your_environment

or::

    $ source activate your_environment

.. WARNING:: Some users may experience errors regarding compatibility if Anaconda already has been configured with certain channels. To resolve this issue you can remove conflicting channels using::

    $ conda config --remove channels [CHANNEL NAME]

  Alternatively, you set your Anaconda channel priority to *flexible*, with the following command::

    $ conda config --set channel_priority flexible

.. WARNING:: After installing `morphMan`, you may experience an error during pre-processing, when the ``.vtu`` mesh is converted and compressed into ``.xml.gz`` format using the ``vmtkMeshWriter`` method.
    As a temporary fix you will need to update the ``vmtkMeshWriter`` script manually to avoid this error, located at ``/Users/[USERNAME]/miniconda3/envs/your_environment/lib/python3.10/site-packages/vmtk/vmtkmeshwriter.py``.
    To apply the fix, open the ``vmtkmeshwriter.py`` file, navigate to `line 264`, and change::

        file = open(self.OutputFileName,'r')

    to::

        file = open(self.OutputFileName,'rb').

Editable installations
======================
If you want to make changes to any of the underlying packages, you should remove them from the `environment.yml` file,
and install then from source, as described for `Oasis` in the next section.

Installing Oasis
################

The next step is to download and install `Oasis <https://github.com/mikaem/Oasis>`_ on your machine.
Inside your conda environment run the following commands to clone and install Oasis, and the minor dependency to `cppimport`::

    $ git clone https://github.com/mikaem/Oasis
    $ pip install cppimport
    $ pip install --editable Oasis

The ``--editable`` flag installs the project in editable mode meaning that any changes to the original package will be reflected directly in your environment.

Docker
======
A `Dockerfile` is supplied at the root of the repository, which can build a docker-image with all dependencies installed.
The Docker-image can be built with the following command::

    $ docker build -t name_of_image .

A Docker-container can then be started with the following command::

    $ docker run -ti --network=host -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $(pwd):/root/shared -w /root/shared --rm --shm-size=512m name_of_image

To run the VaMPy GUI, you need to call::

    $ xhost +local:root


.. WARNING::
    
    on your system before running the scripts. Remember to call::
        
        xhost -local:root
        
    on the host system when you are done running the Docker container.
    
Now you are all set, and can start using the Vascular Modeling Pypeline.
