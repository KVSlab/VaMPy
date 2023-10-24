# Tutorials

These tutorials are meant to guide the user through the basic steps of performing a computational fluid dynamic (CFD)
simulation in a vascular body. The tutorials are divided into three main problems where we will first consider a model
of the [internal carotid artery](https://en.wikipedia.org/wiki/Internal_carotid_artery) (ICA)
with an aneurysm, followed by a model of the [left atrium](https://en.wikipedia.org/wiki/Atrium_(heart)). For the ICA,
we will demonstrate different approaches for meshing and pre-processing of the model, followed by a CFD simulation using
the Oasis software, before post-processing the results. For the left atrium model, we will investigate local refinement
in a user-defined region, also followed by CFD simulation and post-processing of the results. In the third tutorial,
we'll explore a moving domain simulation for a simplified left atrium model, and examine blood residence time by
simultaneously solving a scalar transport equation during the simulation.

We assume that you have VaMPy installed, meaning you have a working version
of [morphMan](https://github.com/KVSlab/morphMan) and [Oasis](https://github.com/mikaem/Oasis) on your computer. For the
moving domain tutorial you will need to install [OasisMove](https://github.com/KVSlab/OasisMove), a modified version
of `Oasis` for moving domain simulations. It is also beneficial, but not necessary, that you are familiar with
the [ParaView](https://www.paraview.org/) visualization software, which we will frequently be using to visualize the
geometries and results.

## Link to tutorials

- [Internal carotid artery](tutorial:artery)
- [Left atrium](tutorial:atrium)
- [Moving domain simulation in the Left atrium](tutorial:movingatrium)
