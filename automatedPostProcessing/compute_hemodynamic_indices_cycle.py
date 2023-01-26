from __future__ import print_function

from os import path
from time import time

from dolfin import *

from postprocessing_common import STRESS, read_command_line

try:
    parameters["reorder_dofs_serial"] = False
except NameError:
    pass


def get_dataset_names(data_file, num_files=3000000, step=1, start=1, print_info=True,
                      vector_filename="/velocity/vector_%d"):
    """
    Read velocity fields datasets and extract names of files

    Args:
        data_file (HDF5File): File object of velocity
        num_files (int): Number of files
        step (int): Step between each data dump
        start (int): Step to start on
        print_info (bool): Prints info about data if true
        vector_filename (str): Name of velocity files

    Returns:
        names (list): List of data file names
    """
    check = True

    # Find start file
    t0 = time()
    while check:
        if data_file.has_dataset(vector_filename % start):
            check = False
            start -= step

        start += step

    # Get names
    names = []
    for i in range(num_files):

        index = start + i * step
        if data_file.has_dataset(vector_filename % index):
            names.append(vector_filename % index)

    t1 = time()

    # Print info
    if MPI.rank(MPI.comm_world) == 0 and print_info:
        print()
        print("=" * 6 + " Timesteps to average over " + "=" * 6)
        print("Length on data set names:", len(names))
        print("Start index:", start)
        print("Wanted num files:", num_files)
        print("Step between files:", step)
        print("Time used:", t1 - t0, "s")
        print()

    return names


def compute_hemodynamic_indices(case_path, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle, step,
                                average_over_cycles):
    """
    Loads velocity fields from completed CFD simulation,
    and computes and saves the following hemodynamic quantities:
    (1) WSS - Wall shear stress
    (2) TAWSS - Time averaged wall shear stress
    (3) TWSSG - Temporal wall shear stress gradient
    (4) OSI - Oscillatory shear index
    (5) RRT - Relative residence time
    (6) ECAP - Endothelial cell activation potential

    The resulting wall shear stress will be in units Pascal [Pa], given that the provided
    density (rho) is in [kg/m^3], the time step (dt) is in [ms], and viscosity (nu) is in [mm^2/ms].

    Args:
        velocity_degree (int): Finite element degree of velocity
        case_path (str): Path to results from simulation
        nu (float): Kinematic viscosity
        rho (float): Fluid density
        dt (float): Time step of simulation
        T (float): One cardiac cycle, in [ms]
        save_frequency (int): Frequency that velocity has been stored
        start_cycle (int): Determines which cardiac cycle to start from for post-processing
        step (int): Step size determining number of times data is sampled
        average_over_cycles (bool): Averages over cardiac cycles if True
    """
    # File paths
    file_path_u = path.join(case_path, "u.h5")
    mesh_path = path.join(case_path, "mesh.h5")
    file_u = HDF5File(MPI.comm_world, file_path_u, "r")

    # Start post-processing from 2nd cycle using every 10th time step, or 2000 time steps per cycle
    start = int(T / dt / save_frequency * (start_cycle - 1))

    # Get names of data to extract
    if MPI.rank(MPI.comm_world) == 0:
        print("Reading dataset names")

    dataset = get_dataset_names(file_u, start=start, step=step)

    # Read mesh saved as HDF5 format
    mesh = Mesh()
    with HDF5File(MPI.comm_world, mesh_path.__str__(), "r") as mesh_file:
        mesh_file.read(mesh, "mesh", False)

    # Load mesh
    bm = BoundaryMesh(mesh, 'exterior')

    if MPI.rank(MPI.comm_world) == 0:
        print("Define function spaces")
    V_b1 = VectorFunctionSpace(bm, "CG", 1)
    U_b1 = FunctionSpace(bm, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)

    if MPI.rank(MPI.comm_world) == 0:
        print("Define functions")
    u = Function(V)

    # RRT
    RRT = Function(U_b1)
    RRT_avg = Function(U_b1)

    # OSI
    OSI = Function(U_b1)
    OSI_avg = Function(U_b1)

    # ECAP
    ECAP = Function(U_b1)
    ECAP_avg = Function(U_b1)

    # WSS_mean
    WSS_mean = Function(V_b1)
    WSS_mean_avg = Function(V_b1)

    # TAWSS
    TAWSS = Function(U_b1)
    TAWSS_avg = Function(U_b1)

    # TWSSG
    TWSSG = Function(U_b1)
    TWSSG_avg = Function(U_b1)
    twssg = Function(V_b1)
    tau_prev = Function(V_b1)

    stress = STRESS(u, 0.0, nu, mesh)

    # Get number of saved steps and cycles
    saved_time_steps_per_cycle = int(T / dt / save_frequency / step)
    n_cycles = int(len(dataset) / saved_time_steps_per_cycle)
    # Set number of cycles to average over
    cycles = list(range(1, n_cycles + 1)) if average_over_cycles else []
    counters_to_save = [saved_time_steps_per_cycle * cycle for cycle in cycles]
    cycle_names = [""] + ["_cycle_{}".format(cycle) for cycle in cycles]

    # Create XDMF files for saving indices
    fullname = file_path_u.replace("u.h5", "%s%s.xdmf")
    fullname = fullname.replace("Solutions", "Hemodynamics")
    index_names = ["RRT", "OSI", "ECAP", "TAWSS", "TWSSG"]
    index_variables = [RRT, OSI, ECAP, TAWSS, TWSSG]
    index_variables_avg = [RRT_avg, OSI_avg, ECAP_avg, TAWSS_avg, TWSSG_avg]

    index_dict = dict(zip(index_names, index_variables))
    index_dict_cycle = dict(zip(index_names, index_variables_avg))

    indices = {}
    for cycle_name in cycle_names:
        for index in index_names + ["WSS"]:
            indices[index + cycle_name] = XDMFFile(MPI.comm_world, fullname % (index, cycle_name))
            indices[index + cycle_name].parameters["rewrite_function_mesh"] = False
            indices[index + cycle_name].parameters["flush_output"] = True
            indices[index + cycle_name].parameters["functions_share_mesh"] = True

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Start post processing", "=" * 10)

    counter = start
    for data in dataset:
        # Update file_counter
        counter += step

        file_u.read(u, data)

        if MPI.rank(MPI.comm_world) == 0:
            timestamp = file_u.attributes(data)["timestamp"]
            print("=" * 10, "Timestep: {}".format(timestamp), "=" * 10)

        # Compute WSS
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS (mean)")
        tau = stress()
        tau.vector()[:] = tau.vector()[:] * rho
        WSS_mean_avg.vector().axpy(1, tau.vector())

        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS (absolute value)")
        tawss = project(inner(tau, tau) ** (1 / 2), U_b1)
        TAWSS_avg.vector().axpy(1, tawss.vector())

        # Compute TWSSG
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute TWSSG")
        twssg.vector().set_local((tau.vector().get_local() - tau_prev.vector().get_local()) / dt)
        twssg.vector().apply("insert")
        twssg_ = project(inner(twssg, twssg) ** (1 / 2), U_b1)
        TWSSG_avg.vector().axpy(1, twssg_.vector())

        # Update tau
        if MPI.rank(MPI.comm_world) == 0:
            print("Update WSS \n")
        tau_prev.vector().zero()
        tau_prev.vector().axpy(1, tau.vector())

        # Save instantaneous WSS
        tau.rename("WSS", "WSS")
        indices["WSS"].write(tau, dt * counter)

        if len(cycles) != 0 and counter == counters_to_save[0]:
            # Get cycle number
            cycle = int(counters_to_save[0] / saved_time_steps_per_cycle)
            if MPI.rank(MPI.comm_world) == 0:
                print("========== Storing cardiac cycle {} ==========".format(cycle))

            # Get average over sampled time steps
            for index in [TWSSG_avg, TAWSS_avg, WSS_mean_avg]:
                index.vector()[:] = index.vector()[:] / saved_time_steps_per_cycle

            # Compute OSI, RRT and ECAP
            wss_mean = project(inner(WSS_mean_avg, WSS_mean_avg) ** (1 / 2), U_b1)
            wss_mean_vec = wss_mean.vector().get_local()
            tawss_vec = TAWSS_avg.vector().get_local()

            # Compute RRT, OSI, and ECAP based on mean and absolute WSS
            RRT_avg.vector().set_local(1 / wss_mean_vec)
            OSI_avg.vector().set_local(0.5 * (1 - wss_mean_vec / tawss_vec))
            ECAP_avg.vector().set_local(OSI.vector().get_local() / tawss_vec)

            for index in [RRT_avg, OSI_avg, ECAP_avg]:
                index.vector().apply("insert")

            # Rename displayed variable names
            for var, name in zip(index_variables_avg, index_names):
                var.rename(name, name)

            # Store solution
            for name, index in index_dict_cycle.items():
                indices[name + "_cycle_{}".format(cycle)].write(index)

            # Append solution to total solution
            for index, index_avg in zip(index_dict.values(), index_dict_cycle.values()):
                index_avg.vector().apply("insert")
                index.vector().axpy(1, index_avg.vector())

            WSS_mean_avg.vector().apply("insert")
            WSS_mean.vector().axpy(1, WSS_mean_avg.vector())

            # Reset tmp solutions
            for index_avg in index_dict_cycle.values():
                index_avg.vector().zero()

            WSS_mean_avg.vector().zero()

            counters_to_save.pop(0)

    print("=" * 10, "Saving hemodynamic indices", "=" * 10)
    n = n_cycles if average_over_cycles else (counter - start) // step
    TWSSG.vector()[:] = TWSSG.vector()[:] / n
    TAWSS.vector()[:] = TAWSS.vector()[:] / n
    WSS_mean.vector()[:] = WSS_mean.vector()[:] / n

    wss_mean = project(inner(WSS_mean, WSS_mean) ** (1 / 2), U_b1)
    wss_mean_vec = wss_mean.vector().get_local()
    tawss_vec = TAWSS.vector().get_local()

    # Compute RRT, OSI, and ECAP based on mean and absolute WSS
    RRT.vector().set_local(1 / wss_mean_vec)
    OSI.vector().set_local(0.5 * (1 - wss_mean_vec / tawss_vec))
    ECAP.vector().set_local(OSI.vector().get_local() / tawss_vec)

    for index in [RRT, OSI, ECAP]:
        index.vector().apply("insert")

    # Rename displayed variable names
    for name, var in index_dict.items():
        var.rename(name, name)

    # Write indices to file
    for name, index in index_dict.items():
        indices[name].write(index)

    print("========== Post processing finished ==========")
    print("Results saved to: {}".format(case_path))


if __name__ == '__main__':
    folder, nu, rho, dt, velocity_degree, _, _, T, save_frequency, _, start_cycle, step \
        = read_command_line()

    compute_hemodynamic_indices(folder, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle,
                                step, average_over_cycles=True)