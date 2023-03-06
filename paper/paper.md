
---
title: 'VaMPy: An Automated and Objective Pipeline for Modeling Vascular Geometries'
tags:
- Python
- pre-processing
- computational fluid dynamics
- post-processing
- vascular modeling
- automated objective pipeline

authors:
- name: Henrik A. Kjeldsberg
  orcid: 0000-0002-7764-4248
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

Hemodynamic forces, such as wall shear stress (WSS), play a role in vessel wall adaption and remodeling. However,
directly measuring these forces is challenging. Medical researchers commonly use image-based computational fluid
dynamics (CFD) to study vascular pathology, but modeling choices and simulation results vary widely. To address this, we
aim to create an automated CFD pipeline for modeling cardiovascular flows that is objective and consistent, where
modeling choices are backed up by rigorous research. The Vascular Modeling Pypeline (`VaMPy`) is an entry-level
high-performance CFD pipeline with high-level Python interface that lets the user easily extend or modify the
functionality.

# Statement of Need

Simulation of the cardiovascular flows has shown to become an indispensable research tool, which can potentially reveal
fundamental properties of the cardiovascular system in both physiological and pathological states. More specifically,
medical image-based computational fluid dynamics (CFD) [@taylor2010image] has been extensively used in the investigation
of disease initiation of, e.g., coronary artery disease [@taylor2013computational], carotid bifurcation [@fisher],
arteriovenous fistula [@Lee2005a], and aneurysms [@steinman2003image]. Numerous scientific studies have been conducted
to scrutinize and assess different components of a conventional image-based modeling process, with the objective of
creating a genuinely "patient-specific" CFD model. As highlighted and examined in a review [@steinman2019patient],
particular emphasis has been placed on investigating the influence of medical imaging techniques, segmentation methods,
flow velocities, and the impact of non-Newtonian rheology. However, recent challenge studies withinin aneurysm research
have brought to light a significant inter-laboratory variability. When 26 research groups were provided with identical
segmented surfaces and boundary conditions, the results showed large variability stemming from the various CFD solution
strategies [@Steinman2013]. Futhermore, when provided with identical medical images and no guidelines – reflecting
current research practice – results from 28 research groups showed significant variaiblity in the predicted
WSS [@valen2018real]. These results might point to a broader reproducibility issue. While modeling and simulating
cardiovascular flow can provide valuable and additional insight to vascular remodeling, establishing local computational
pipelines for medical image-based CFD remains a time-intensive process that is error-prone and a significant source of
variability.

With this in mind, the objective was to devise a comprehensive and resilient open-source research code enabling
open and reproduciable science, with an emphasis on user-friendliness, geared towards students, educators, and researchers. By
automating the process, we reduce the need for manual labor, which enables mass production of CFD results, and of equal
importance, significantly reduces the variability. The latter is also ensured by making all aspects of the modeling
choices based on state-of-the-art research shown to be the current gold-standard choices in aneurysm CFD modeling.
Thus, `VaMPy`
enables non-CFD-experts to perform objective and automated out-of-the-box CFD simulations, and to produce results of
publication quality.

![
Illustration of the computational fluid dynamics pipeline be executing the `Artery.py` demo in `VaMPy`. From left to right: (1) volumetric mesh, (2) boundary layers, (3) boundary conditions (flow rate at the inlet and flow split at the outlets) and probe points along the computational domain for which the velocity and pressure is evaluated during the simulation, (4) instantaneous velocity filed visualised with `Paraview` using glyphs, and (5) the resulting time averaged wall shear stress (TAWSS). \label{fig:pipeline}](Figure1.png)

# Overview of features

The first feature of `VaMPy` is the pre-processing pipeline, which is built upon the methods introduced
in `morphMan` [@kjeldsberg2019morphman], a published framework for objective and automated manipulation of vascular
morphologies. The pre-processing pipeline includes volumetric mesh generation, automated identification of inlet and
outlets for boundary conditions, and generation of probe points for velocity and pressure measurements within the
domain. Prior to meshing, the user may also add adjustable features such as flow extensions, surface smoothing, local
refinement, and generation of boundary layers. The volumetric meshing may be set to uniform or variable mesh density,
and in the two leftmost panels of \autoref{fig:pipeline} we show a meshed artery model with boundary layers and flow
extensions. Following the mesh generation, flow split boundary conditions are generated, and probe points are stored,
visualized in the middle panel of \autoref{fig:pipeline}.

The second feature of `VaMPy` is the CFD simulation pipeline, based on the solver `Oasis` [@oasis], which has been
verified and validated against spectral element methods and experimental
measurements [@khan2019direct; @bergersen2019fda]. `Oasis` is an open-source, finite element-based segregated
high-performance computing implementation of a space/time centered incremental pressure correction scheme. `Oasis` is
formal second-order accurate in time that ensures a solution that preserves kinetic energy while minimizing numerical
dispersion and diffusion errors [@Karniadakis2005]. A Womersley profile is prescribed at the inlet, where the inflow
waveform was obtained from older adults [@Hoi2011]. We prescribe the flow rate according to the square law, which
results in an average internal carotid artery flow rate of 245 mL/min for average sized arteries [@valen2015estimation].
At the outlets, we use a reduced order method to split the flow [@chnafa2017improved].

The third feature of `VaMPy` is post-processing, where we have scripts that compute the flow and simulation metrics,
hemodynamic indices, probe point visualization, and velocity and pressure conversion. The flow metrics include
parameters such as the friction velocity (and associated $l^+$ and $t^+$ values) [@valen2011direct], which allows the
user to assess the relative resolution and simulation quality. The script also computes the Kolmogorov scales, kinetic
energy, and turbulent kinetic energy, based on phase averaging multiple cardiac cycles. The script for computing
hemodynamic indices includes the most commonly computed ones, including wall shear stress (WSS), oscillatory shear
index (OSI), and relative residence time (RRT), and to demonstrate we have shown the time averaged WSS (TAWSS) in the
rightmost panel of \autoref{fig:pipeline}. The probe point visualization script creates a figure of velocity and
pressure traces at pre-determined points inside the domain. Finally, the conversion script creates viewable versions of
the compact velocity and pressure solutions, and may be visualized in software such as `ParaView` [@ayachit2015paraview].

![Example of an extension of `VaMPy` to cardiovascular flow in the left atrium and the associated hemodynamic stresses. From top left to bottom right: the volumetric rendering of velocity, the pressure field, volumetric rendering of the Q-criterion, TAWSS, OSI, and RRT. \label{fig:atrium}](Figure2.png)

# Extension to cardiac flows

The pipeline is fully automated and has been demonstrated and tailored towards simulations of cerebrovascular flows. The
demonstration shown in \autoref{fig:pipeline} is configured to be run on a laptop within a reasonable time frame, but to
perform simulations with adequate resolutions we refer to [@valen2014high; @valen2014mind; @khan2015narrowing]. Beside
cerebrovascular flows, `VaMPy` can easily be extended to also allow for simulation of other vascular territories. In
\autoref{fig:atrium} we show the application to modeling of the left atrium. More specifically, from top left to bottom
right the figure shows the instantaneous velocity magnitude, instantaneous pressure, vortex cores (Q-criterion), and the
time averaged quantities wall shear stress (TAWSS), oscillatory shear index (OSI), and relative residence time (RRT),
all of which are computed with the hemodynamics post-processing script.

# Acknowledgements

This work was supported by the SimCardioTest project (Digital transformation in Health and Care SC1-DTH-06-2020) under
grant agreement No. 101016496 and ERACoSysMed PARIS project under grant agreements No. 643271. The simulations were
performed on the Saga cluster, with resources provided by UNINETT Sigma2 – the National Infrastructure for High
Performance Computing and Data Storage in Norway, grant number nn9249k. We wish to thank Dr. Jørgen Dokken for technical
assistance, and Dr. Christophe Chnafa for his contribution to the boundary condition methodology.

# References