# Computational fluid dynamics 
## Simulations in `Oasis`

Following pre-processing, the next step of using the Vascular Modeling Pypeline is performing the computational fluid dynamics (CFD) simulations with `oasis`. 
Assuming `oasis` has been installed, start by navigating to the `simulation` folder:

``` console
$ cd src/vampy/simulation
```

We can now perform a CFD simulation for two cycles with 10 000 time steps per cycle and default parameters by executing the following command:
 
``` console
$ oasis NSfracStep problem=Artery mesh_path=../../../tests/Case_test_artery/artery.xml.gz save_solution_after_cycle=0
```

Running the simulations will create the result folder `results_artery` (specific to the `Artery.py` problem) located inside `src/vampy/simulation`, with the results and corresponding mesh saved compactly in HDF5 format.

## Adjusting simulation parameters

The default parameters for CFD simulation have been chosen based on the authors' experience and clinically reported hemodynamic parameters.
However, changing any of the parameters is simply done by adding them as command line arguments.
To demonstrate, consider the following Python snippet showing an overview of the parameters used for the `Artery.py` problem:

``` Python
# Parameters are in mm and ms
cardiac_cycle = float(commandline_kwargs.get("cardiac_cycle", 951))
number_of_cycles = float(commandline_kwargs.get("number_of_cycles", 2))

NS_parameters.update(
    # Fluid parameters
    nu=3.3018e-3,  # Kinematic viscosity: 0.0035 Pa-s / 1060 kg/m^3 = 3.3018E-6 m^2/s = 3.3018-3 mm^2/ms
    # Geometry parameters
    id_in=[],  # Inlet boundary ID
    id_out=[],  # Outlet boundary IDs
    area_ratio=[],  # Area ratio for the flow outlets
    area_inlet=[],  # Area of inlet in [mm^2]
    # Simulation parameters
    cardiac_cycle=cardiac_cycle,  # Duration of cardiac cycle [ms]
    T=cardiac_cycle * number_of_cycles,  # Simulation end time [ms]
    dt=0.0951,  # Time step size [ms]
    dump_probe_frequency=100,  # Dump frequency for sampling velocity & pressure at probes along the centerline
    save_solution_frequency=5,  # Save frequency for velocity and pressure field
    save_solution_after_cycle=1,  # Store solution after 1 cardiac cycle
    # Oasis specific parameters
    checkpoint=500,  # Checkpoint frequency
    print_intermediate_info=100,  # Frequency for printing solver statistics
    folder="results_artery",  # Preferred results folder name
    mesh_path=commandline_kwargs["mesh_path"],  # Path to the mesh
    # Solver parameters
    velocity_degree=1,  # Polynomial order of finite element for velocity. Normally linear (1) or quadratic (2)
    pressure_degree=1,  # Polynomial order of finite element for pressure. Normally linear (1)
    use_krylov_solvers=True,
    krylov_solvers=dict(monitor_convergence=False)
)
```
To start a simulation that runs for five cardiac cycles, use a coarser time step, and saves the solution less frequent we can run the following command:

``` console
$ oasis NSfracStep problem=Artery mesh_path=../../../tests/Case_test_artery/artery.xml.gz number_of_cycles=5 dt=0.951 save_solution_frequency=20
```
