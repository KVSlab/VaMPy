# VaMPy - Vascular Modeling Pypeline

_________________
[![GPL-3.0](https://img.shields.io/github/license/hkjeldsberg/vampy)](LICENSE)
[![codecov](https://codecov.io/gh/KVSlab/VaMPy/branch/master/graph/badge.svg?token=M2NMX6HOSZ)](https://codecov.io/gh/KVSlab/VaMPy)
[![CI](https://github.com/kvslab/vampy/actions/workflows/check_and_test_package.yml/badge.svg)](https://github.com/kvslab/vampy/actions/workflows/check_and_test_package.yml)
[![GitHub pages](https://github.com/kvslab/vampy/actions/workflows/deploy_pages.yml/badge.svg)](https://github.com/kvslab/vampy/actions/workflows/deploy_pages.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7707855.svg)](https://doi.org/10.5281/zenodo.7707855)

_________________

<p align="center">
    <img src=docs/figures/artery_pipeline.png width="830 height="370" alt="Output pre processing"/>
</p>
<p align="center">
    Pre-processed and simulated artery model. From left to right: A variable density volumetric mesh, zoomed in view of an outlet showing the four boundary layers, corresponding inlet flow rate, outlet flow split, and probes for velocity and pressure sampling. From the simulation results, we have shown the velocity field represented by vectors, and the time averaged wall shear stress (TAWSS) as one of the many hemodynamic indices computed by the post-processing scripts.
</p>

Description
-----------
The Vascular Modeling Pypeline (VaMPy) is a collection of fully automated scripts used to prepare, run, and analyze
cardiac and atrial morphologies. This includes pre-processing scripts for meshing and probe sampling,
two [Oasis](https://github.com/mikaem/Oasis) problem files for simulating flow in
the [internal carotid artery](https://en.wikipedia.org/wiki/Internal_carotid_artery) and
the [left atrium](https://en.wikipedia.org/wiki/Atrium_(heart)), and a variety of post-processing scripts for computing
WSS-based metrics, more advanced turbulence metrics, and a variety of morphological parameters in patient-specific
geometries.

Installation
------------
VaMPy is a Python package for Python >= 3.8, with main dependencies to [morphMan](https://github.com/KVSlab/morphMan)
and [Oasis](https://github.com/mikaem/Oasis). VaMPy and its dependencies can be installed with `conda` on Linux and
macOS as explained [here](https://kvslab.github.io/VaMPy/conda.html). The package can also be installed and run through
its latest `Docker` image supported by Windows, Linux, and macOS, and
explained [here](https://kvslab.github.io/VaMPy/docker.html).

Documentation
-------------
VaMPy's documentation is hosted [here](https://kvslab.github.io/VaMPy). This includes
two [tutorials](https://kvslab.github.io/VaMPy/tutorials.html), meant to guide the user through the basic steps of
performing a computational fluid dynamic simulation in a vascular body.

Licence
-------
VaMPy is licensed under the GNU GPL, version 3 or (at your option) any later version.

VaMPy is Copyright (2018-2023) by the authors.

Authors
-------
VaMPy has been developed by

* [Aslak W. Bergersen](https://github.com/aslakbergersen)
* [Henrik A. Kjeldsberg](https://github.com/HKjeldsberg)

Issues
------
Please report bugs and other issues through the issue tracker at:

https://github.com/KVSlab/VaMPy/issues
