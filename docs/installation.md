# Installation

`VaMPy` is a pure Python package that combines the morphological processing tools of
[morphMan](https://github.com/KVSlab/morphMan) and [vmtk](http://www.vmtk.org/), the computational fluid dynamics
solver [Oasis](https://github.com/mikaem/Oasis), and the finite element computing
platform [FEniCS](https://fenicsproject.org/) into a powerful and user-friendly computational fluid dynamics pipeline.

## Dependencies

The dependencies of `VaMPy` are:

* Python >= 3.8
* FEniCS >= 2018.1
* morphMan >= 1.2
* vmtk >= 1.4.0
* Oasis >= 2018.1
* paramiko >= 3.0
* cppimport >= 22.8.2

## Installing with `conda`

To install `VaMPy` and all its dependencies to a *local environment*, we recommend using `conda`.
Instructions for installing `VaMPy`
with `conda` can be found [here](install:conda).

## Installing with Docker

To install `VaMPy` and all its dependencies to an *isolated environment*, we recommend using the dedicated Docker
container. Instructions for installing `VaMPy` with Docker can be found [here](install:docker).

