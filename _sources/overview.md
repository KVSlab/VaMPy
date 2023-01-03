# Overview

In this overview, we present VaMPy's functionality by
preparing, simulating, and post-processing an [internal carotid
artery](https://en.wikipedia.org/wiki/Internal_carotid_artery) geometry, although
the software may be readily used for other tubular or vascular shapes.

## Pre-processing: Meshing and boundary conditions

The first step of using the Vascular Modeling Pypeline is
pre-processing. The pre-processing scripts are located inside the
`automatedPreProcessing` folder, and we will be executing the
`automatedPreProcessing.py` script to generate a mesh, boundary
conditions, and probes for velocity and pressure sampling. Here we will
perform pre-processing for the artery case located in the `test` folder.
Start by entering the conda environment you created:

``` console
$ conda activate your_environment
```

Then, to perform meshing, execute the following command:

``` console
$ python automatedPreProcessing/automatedPreProcessing.py -m diameter -i test/Case_test_artery/artery.vtp -c 1.3
```

When complete, the script will save the volumetric mesh as `artery.vtu`,
alongside a compressed DOLFIN mesh in `artery.xml.gz`, used for the
following simulations. The pre-processing script will also produce an
info file and a probe file, named `artery_info.json` and
`artery_probes`, respectively.

## Computational fluid dynamics simulations in Oasis

The next step of using the Vascular Modeling Pypeline is performing the
CFD simulations with `Oasis`. 
For convenience, change directory to the `simulation` folder:

``` console
$ cd simulation
```

Then, to run a CFD simulation for two cycles with 10 000 time steps per
cycle and default parameters with Oasis, execute the following command:

``` console
$ oasis NSfracStep problem=Artery mesh_path=../test/Case_test_artery/artery.xml.gz T=9.61 save_solution_after_cycle=0
```

Running the simulations will create the result folder `results_artery`
(specific to the artery problem), with the results and corresponding
mesh saved compactly in HDF5 format.

(overview:post)=
## Post-processing: Hemodynamic indices, flow and simulation metrics, velocity and pressure, and probes

Following the CFD simulations, the last step of the Vascular Modeling
Pypeline is post-processing of the results. For reference, the
quantities computed by the scripts `compute_hemodynamic_indices.py` and
`compute_flow_and_simulation_metrics.py` can be viewed
[here](post:quantities).

You can start by computing the wall shear stress, oscillatory shear
index and other hemodynamic indices by executing the following command,
executed from the root directory:

``` console
$ python automatedPostProcessing/compute_hemodynamic_indices.py --case simulation/results_artery/artery/data/[RUN_NUMBER]/Solutions
```

To compute fluid dynamic quantities and simulation metrics, you may
execute the following command:

``` console
$ python automatedPostProcessing/compute_flow_and_simulation_metrics.py --case simulation/results_artery/artery/data/[RUN_NUMBER]/Solutions
```

As default, `compute_flow_and_simulation_metrics.py` computes the
average values over full cycle, for a desired number of cycles,
determined by the `--start-cycle` flag. Setting `--start-cycle 2` would
correspond to computing the averaged values from the second cardiac
cycle and onward. Alternatively, the user may supply the specific times
$t$ during the cardiac cycle $T$ to the `--times-to-average` flag, to
compute phase averaged quantities, where $t \in [0,T)$. Thus, to compute
the metrics at $t=0.2$ s representing peak systole, and to skip the
first cycle, the user may run the following command:

``` console
$ python automatedPostProcessing/compute_flow_and_simulation_metrics.py --case simulation/results_artery/artery/data/[RUN_NUMBER]/Solutions --start-cycle 2 --times-to-average 200 --T 1000 --dt 0.1
```

Note that the specified time is in milliseconds. A comparison between
cycle averaged and phase averaged kinetic energy (KE) and turbulent
kinetic energy (TKE) for a representative artery model is shown in
{numref}`phase`. The leftmost panels displays the model, specifically case
C0097 from the Aneurisk database. The top middle panel displays the
cycle averaged KE over the last four cycles, whereas the bottom middle
panel displays the phase averaged KE at $t=0.2$ s. Similarly, the top
right panel displays the cycle averaged TKE over the last four cycles,
and bottom right panel displays the phase averaged TKE at $t=0.2$ s. For
this illustration, the simulation was run over 5 cycles, and the metrics
were computed over the last four cycles.

```{figure} figures/phase_averaged.png
---
name: phase
---
Cycle and phase average kinetic energy and turbulent kinetic
energy. The leftmost panel disaply the model, the middle panels display
kinetic energy, and the rightmost panel display turbulent kinetic
energy.
```

Because velocity and pressure results are currently saved in the
compressed `.h5` format, they are not viewable in a software such as
[ParaView](https://www.paraview.org/). If it is desired to convert the
compressed velocity and pressure results to viewable `.xdmf` format, you
may execute the following command:

``` console
$ python automatedPostProcessing/compute_velocity_and_pressure.py --case simulation/results_artery/artery/data/[RUN_NUMBER]/Solutions
```

which will store the velocity and pressure into `velocity.xdmf` and
`pressure.xdmf`, respectively.

Finally, to visualize velocity and pressure at the probes created by
`Artery.py`, you can run the `visualize_probes.py` script, by executing
the following command:

``` console
$ python automatedPostProcessing/visualize_probes.py --case simulation/results_artery/artery/data/[run_number]/Probes
```

Note that this has an additional dependency to
[Matplotlib](https://github.com/matplotlib/matplotlib), which can
quickly be installed with either `conda` or `pip`.

