---
title: 'VaMPy: Automated pipeline for modelling vascular geometries'
tags:
  - Python 
  - pre-processing
  - fluid simulation
  - post-processing
  - vascular geometries
authors:
  - name: Henrik A. Kjeldsberg
  - orcid: 0000-0002-7764-4248
    affiliation: 1
  - name: Aslak W. Bergersen
    orcid: 0000-0001-5063-3680
    affiliation: 1
  - name: Kristian Valen-Sendstad
    orcid: 0000-0002-2907-0171
    affiliation: 1
affiliations:
  - name: Department of Computational Physiology, Simula Research Laboratory
    index: 1
date: 15 March 2023
bibliography: paper.bib
---

# Summary

Cardiovascular diseases are overwhelming the healthcare systems, and the
costs are anticipated to increase in the years to come [@Murray1997],
not to the mention the personal tragedy for those affected [@gage1996effect].
Systemic risk factors are well known to correlate with cardiovascular diseases in general,
but, for instance, arterial plaques and brain aneurysms are focalized, highlighting
the role of local hemodynamics. Furthermore, blood-flow induced wall shear stress (WSS) is
known to contribute to vessel wall adaption and remodeling [@Malek1999b; @morbiducci2016atherosclerosis],
but is challenging to measure *in-vivo*. On the other hand, medical images are routinely available and have
been extensively used in combination with computational fluid dynamics to
study the initiation, progression, and outcome of vascular pathologies [@taylor2010image].

![voronoi centerline](./figure1.png)\

**Figure 1:** 
   A visualization of the Voronoi diagram (left) and the centerline (right) of a surface.

- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.
- A clear statement of need that illustrates the purpose of the software.
- A description of how this software compares to other commonly-used packages in this research area.
- Mentions (if applicable) of any ongoing research projects using the software or recent scholarly publications enabled by it.
- A list of key references including a link to the software archive.

# Acknowledgements

We acknowledge XXX for testing VaMPy, and the open-source projects [vtk](https://www.vtk.org/), [vmtk](http://www.vmtk.org), and [FEniCS](https://fenicsproject.org).

# References