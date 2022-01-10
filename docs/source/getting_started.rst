.. title:: Using VaMPy

.. _getting_started:

====================================
Using the Vascular Modeling Pypeline
====================================
The Vascular Modeling Pypeline is a collection of scripts to prepare, run, and analyze vascular morphologies. This includes pre-processing scripts for generating a volumetric mesh, defining physiological boundary conditions, and inserting probes for sampling velocity and pressure. For the computational fluid dynamics (CFD) simulation, we have included an artery problem file used for running the simulation with `Oasis <https://github.com/mikaem/Oasis>`_.
Finally, there are a variety of post-processing scripts, which computes wall shear stress-based metrics, more advanced turbulence metrics, and a variety of morphological parameters. In this walkthrough, we exemplify the usage by preparing, simulating, and post-processing the `internal carotid artery <https://en.wikipedia.org/wiki/Internal_carotid_artery>`_, although the software may be readily used for other tubular or vascular shapes.

Pre-processing: Meshing and boundary conditions
===============================================
The first step of using the Vascular Modeling Pypeline is pre-processing. The pre-processing scripts are located inside the ``automatedPreProcessing`` folder, and we will be executing the ``automatedPreProcessing.py`` script to generate a mesh, boundary conditions, and probes for velocity and pressure sampling. Here we will perform pre-processing for the artery case located in the ``test`` folder.
Start by entering the ``morphman`` conda environment::

    conda deactivate && conda activate morphman

Then, to perform meshing, execute the following command::

    python automatedPreProcessing/automatedPreProcessing.py -m diameter -i test/Case_test_artery/artery.vtp -c 1.3

When complete, the script will save the volumetric mesh as ``artery.vtu``, alongside a compressed DOLFIN mesh in ``artery.xml.gz``, used for the following simulations.
The pre-processing script will also produce an info file and a probe file, named ``artery_info.json`` and ``artery_probes``, respectively.

Computational fluid dynamics simulations in Oasis
=================================================
The next step of using the Vascular Modeling Pypeline is performing the CFD simulations with `Oasis`.
First, activate the ``fenics`` conda environment::

    conda deactivate && conda activate fenics && cd simulation

Then, to run a CFD simulation for two cycles with 10 000 time steps per cycle and default parameters with Oasis, execute the following command::

    oasis NSfracStep problem=Artery mesh_path=../test/Case_test_artery/artery.xml.gz T=9.61 save_solution_after_cycle=0 && cd ..

Running the simulations will create the result folder ``results_artery`` (specific to the artery problem), with the results and corresponding mesh saved compactly in HDF5 format.

Post-processing: Hemodynamic indices, flow and simulation metrics, and probes
=============================================================================
Following the CFD simulations, the last usage of the Vascular Modeling Pypeline is the post-processing part.
You can start by computing the wall shear stress, oscillatory shear index and other hemodynamic indices by executing the following command::

    python automatedPostProcessing/compute_hemodynamic_indices.py --case simulation/results_artery/artery/data/[RUN_NUMBER]/Solutions

To compute fluid dynamic quantities and simulation metrics, you may execute the following command::

    python automatedPostProcessing/compute_flow_and_simulation_metrics.py --case simulation/results_artery/artery/data/[RUN_NUMBER]/Solutions

Finally, to visualize velocity and pressure at the probes created by ``Artery.py``, you can run the ``visualize_probes.py`` script, by executing the following command::

    python automatedPostProcessing/visualize_probes.py --case simulation/results_artery/artery/data/[run_number]/Probes

Note that this has an additional dependency to `Matplotlib <https://github.com/matplotlib/matplotlib>`_, which can quickly be installed with either `conda` or `pip`.

Features and issues
===================
The existing methods provide many degrees of freedom, however, if you need a specific method or functionality, please do not hesitate to propose enhancements in the `issue tracker <https://github.com/KVSlab/VaMPy/issues>`_, or create a `pull request <https://github.com/KVSlab/VaMPy/pulls>`_ with new features.
Similarly, we highly appreciate that you report any bugs or other issues you may experience in the `issue tracker <https://github.com/KVSlab/VaMPy/issues>`_.

