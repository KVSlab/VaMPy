from __future__ import print_function

from os import path

from dolfin import *

from postprocessing_common import STRESS, read_command_line, get_dataset_names

try:
    parameters["reorder_dofs_serial"] = False
except NameError:
    pass


def compute_hemodynamic_indices(case_path, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle, step):
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

    # OSI
    OSI = Function(U_b1)

    # ECAP
    ECAP = Function(U_b1)

    # WSS_mean
    WSS_mean = Function(V_b1)
    wss_mean = Function(U_b1)

    # TAWSS
    TAWSS = Function(U_b1)
    tawss = Function(U_b1)

    # TWSSG
    TWSSG = Function(U_b1)
    twssg_ = Function(U_b1)
    twssg = Function(V_b1)
    tau_prev = Function(V_b1)

    stress = STRESS(u, 0.0, nu, mesh)

    # Create XDMF files for saving indices
    fullname = file_path_u.replace("u.h5", "%s.xdmf")
    fullname = fullname.replace("Solutions", "Hemodynamics")
    index_names = ["WSS", "RRT", "OSI", "ECAP", "TAWSS", "TWSSG"]
    index_variables = [RRT, OSI, ECAP, TAWSS, TWSSG]

    indices = {}
    for index in index_names:
        indices[index] = XDMFFile(MPI.comm_world, fullname % index)
        indices[index].parameters["rewrite_function_mesh"] = False
        indices[index].parameters["flush_output"] = True
        indices[index].parameters["functions_share_mesh"] = True

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Start post processing", "=" * 10)

    counter = start
    for data in dataset:
        file_u.read(u, data)

        if MPI.rank(MPI.comm_world) == 0:
            timestamp = file_u.attributes(data)["timestamp"]
            print("=" * 10, "Timestep: {}".format(timestamp), "=" * 10)

        # Compute WSS
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS (mean)")
        tau = stress()
        tau.vector()[:] = tau.vector()[:] * rho
        WSS_mean.vector().axpy(1, tau.vector())

        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS (absolute value)")
        tawss = project(inner(tau, tau) ** (1 / 2), U_b1)
        TAWSS.vector().axpy(1, tawss.vector())

        # Compute TWSSG
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute TWSSG")
        twssg.vector().set_local((tau.vector().get_local() - tau_prev.vector().get_local()) / dt)
        twssg.vector().apply("insert")
        twssg_ = project(inner(twssg, twssg) ** (1 / 2), U_b1)
        TWSSG.vector().axpy(1, twssg_.vector())

        # Update tau
        if MPI.rank(MPI.comm_world) == 0:
            print("Update WSS \n")
        tau_prev.vector().zero()
        tau_prev.vector().axpy(1, tau.vector())

        # Save instantaneous WSS
        tau.rename("WSS", "WSS")
        indices["WSS"].write(tau, dt * counter)

        # Update file_counter
        counter += step

    print("=" * 10, "Saving hemodynamic indices", "=" * 10)
    n = (counter - start) // step
    TWSSG.vector()[:] = TWSSG.vector()[:] / n
    TAWSS.vector()[:] = TAWSS.vector()[:] / n
    WSS_mean.vector()[:] = WSS_mean.vector()[:] / n

    TAWSS.rename("TAWSS", "TAWSS")
    TWSSG.rename("TWSSG", "TWSSG")

    wss_mean = project(inner(WSS_mean, WSS_mean) ** (1 / 2), U_b1)
    wss_mean_vec = wss_mean.vector().get_local()
    tawss_vec = TAWSS.vector().get_local()

    # Compute RRT and OSI based on mean and absolute WSS
    RRT.vector().set_local(1 / wss_mean_vec)
    RRT.vector().apply("insert")
    RRT.rename("RRT", "RRT")

    OSI.vector().set_local(0.5 * (1 - wss_mean_vec / tawss_vec))
    OSI.vector().apply("insert")
    OSI.rename("OSI", "OSI")

    # Compute ECAP based on OSI and TAWSS
    ECAP.vector().set_local(OSI.vector().get_local() / tawss_vec)
    ECAP.vector().apply("insert")
    ECAP.rename("ECAP", "ECAP")

    for writer, value in zip(list(indices.values())[1:], index_variables):
        writer.write(value)

    print("========== Post processing finished ==========")
    print("Results saved to: {}".format(case_path))


if __name__ == '__main__':
    folder, nu, rho, dt, velocity_degree, _, _, T, save_frequency, _, start_cycle, step = read_command_line()
    compute_hemodynamic_indices(folder, nu, rho, T, dt, velocity_degree, save_frequency, start_cycle, step)
