from os import path
from time import time

import numpy as np
from dolfin import *

from postprocessing_common import read_command_line, epsilon

set_log_active(False)


def get_dataset_names(data_file, num_files=3000000, step=1, start=1, print_info=True, vector_filename="/velocity/vector_%d"):
    """
    Read velocity fields datasets and extract names of files

    Args:
        data_file (HDF5File): File object of velocity
        num_files (int): Number of files
        step (int): Step between each data dump
        start (int): Step to start on
        print_info (bool): Prints info about data if true
        velocity_filename (str): Name of velocity files

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
        step = 1
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


def rate_of_strain(ssv, u, v, mesh, h):
    """
    Computes rate of strain

    Args:
        ssv (Function): Function to save rate of strain to
        u (Function): Function for velocity field
        v (Function): Test function for velocity
        mesh: Mesh to compute strain rate on
        h (float): Cell diameter of mesh

    Returns:
        ssv (Function): Rate of strain
    """
    eps = epsilon(u)
    f = sqrt(inner(eps, eps))
    x = assemble(inner(f, v) / h * dx(mesh))
    ssv.vector().set_local(x.get_local())
    ssv.vector().apply("insert")

    return ssv


def rate_of_dissipation(ssv, u, v, mesh, h, nu):
    """
    Computes rate of dissipation

    Args:
        ssv (Function): Function to save rate of dissipation to
        u (Function): Function for velocity field
        v (Function): Test function for velocity
        mesh: Mesh to compute dissipation rate on
        h (float): Cell diameter of mesh
        nu (float): Viscosity

    Returns:
        ssv (Function): Rate of dissipation
    """
    eps = epsilon(u)
    f = 2 * nu * inner(eps, eps)
    x = assemble(inner(f, v) / h * dx(mesh))
    ssv.vector().set_local(x.get_local())
    ssv.vector().apply("insert")

    return ssv


def compute_flow_metrics(folder, nu, dt):
    """
    Computes several flow field characteristics
    for velocity field stored at 'folder' location
    for flow_metrics given viscosity and time step

    Args:
        folder (str): Path to simulation results
        nu (float): Viscosity
        dt (float): Time step

    """
     # File paths

    file_path_u = path.join(folder, "u.h5")
    mesh_path = path.join(folder, "mesh.h5")

    f = HDF5File(MPI.comm_world, file_path_u, "r")

    # Get names of data to extract
    start = 0
    if MPI.rank(MPI.comm_world) == 0:
        print("The post processing starts from", start)

    dataset_names = get_dataset_names(f, start=start)

    # Get mesh information
    mesh = Mesh()
    with HDF5File(MPI.comm_world, mesh_path, "r") as mesh_file:
        mesh_file.read(mesh, "mesh", False)

    # Function space
    velocity_degree = 1
    DG = FunctionSpace(mesh, 'DG', 0)
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)
    Vv = FunctionSpace(mesh, "CG", velocity_degree)

    # Mesh info
    h = CellDiameter(mesh)
    dl = project(h, DG)

    # Functions for storing values
    v = TestFunction(DG)
    u = Function(V)
    u_mean = Function(V)
    u_prime = Function(V)

    # pluss-values
    l_pluss_avg = Function(DG)
    l_pluss = Function(DG)
    t_pluss_avg = Function(DG)
    t_pluss = Function(DG)

    # Kolmogorov scales
    length_scale = Function(DG)
    length_scale_avg = Function(DG)
    time_scale = Function(DG)
    time_scale_avg = Function(DG)
    velocity_scale = Function(DG)
    velocity_scale_avg = Function(DG)

    # Inner grad(u), grad(u)
    turbulent_dissipation = Function(DG)
    turbulent_dissipation_avg = Function(DG)
    ssv = Function(DG)
    ssv_avg = Function(DG)
    dissipation = Function(DG)
    dissipation_avg = Function(DG)

    # Energy
    kinetic_energy = Function(Vv)
    kinetic_energy_avg = Function(Vv)
    turbulent_kinetic_engergy = Function(Vv)
    turbulent_kinetic_engergy_avg = Function(Vv)

    u0 = Function(Vv)
    u1 = Function(Vv)
    u2 = Function(Vv)
    u0_prime = Function(Vv)
    u1_prime = Function(Vv)
    u2_prime = Function(Vv)

    # CFL
    CFL = Function(DG)
    CFL_avg = Function(DG)

    # Create xdmf files
    fullname = file_path_u.replace("u.h5", "%s.xdmf")
    fullname = fullname.replace("Solutions", "mesh_info_testTKE")
    var_name = ["u_mean", "l_pluss", "t_pluss", "CFL", "ssv", "length_scale", "time_scale", "velocity_scale", "u_mag",
                "dl", "dissipation", "kinetic_energy", "turbulent_kinetic_engergy", "turbulent_dissipation", "u_prime",
                "u_viz"]
    # TODO: Add WSS, RRT, OSI, TWSSG, WSSG
    metrices = {}
    for vn in var_name:
        if MPI.rank(MPI.comm_world) == 0:
            print(fullname % vn)
        metrices[vn] = XDMFFile(MPI.comm_world, fullname % vn)
        metrices[vn].parameters["rewrite_function_mesh"] = False
        metrices[vn].parameters["flush_output"] = True

    metrices["dl"].write(dl)

    # Get u mean
    u_mean_file_path = file_path_u.replace("u.h5", "u_mean.h5")
    tmp_file = HDF5File(MPI.comm_world, u_mean_file_path, "r")
    tmp_file.read(u, "u_mean/vector_0")
    tmp_file.close()
    assign(u_mean, u)
    
    for data in dataset_names:
        if MPI.rank(MPI.comm_world) == 0:
            print(data)

        # Time step and velocity
        f.read(u, data)
      
        # Compute CFL
        u_mag = project(sqrt(inner(u, u)), DG)
        CFL.vector().set_local(u_mag.vector().get_local() / dl.vector().get_local() * dt)
        CFL.vector().apply("insert")
        CFL_avg.vector().axpy(1, CFL.vector())

        # Compute rate-of-strain
        rate_of_strain(ssv, u, v, mesh, h)
        ssv_avg.vector().axpy(1, ssv.vector())

        # Compute l+ 
        u_star = np.sqrt(ssv.vector().get_local() * nu)
        l_pluss.vector().set_local(u_star * dl.vector().get_local() / nu)
        l_pluss.vector().apply("insert")
        l_pluss_avg.vector().axpy(1, l_pluss.vector())

        # Compute t+
        t_pluss.vector().set_local(nu / u_star ** 2)
        t_pluss.vector().apply("insert")
        t_pluss_avg.vector().axpy(1, t_pluss.vector())

        # Compute Kolmogorov
        rate_of_dissipation(dissipation, u, v, mesh, h, nu)
        ssv_ = dissipation.vector().get_local()
        dissipation_avg.vector().axpy(1, dissipation.vector())

        # Compute u_prime 
        u_prime.vector().set_local(u.vector().get_local() - u_mean.vector().get_local())
        u_prime.vector().apply("insert")

        # Compute Turbulent dissipation
        rate_of_dissipation(turbulent_dissipation, u_prime, v, mesh, h, nu)
        turbulent_dissipation_avg.vector().axpy(1, turbulent_dissipation.vector())

        # Compute length scale
        length_scale.vector().set_local((nu ** 3 / ssv_) ** (1. / 4))
        length_scale.vector().apply("insert")
        length_scale_avg.vector().axpy(1, length_scale.vector())

        # Compute time scale
        time_scale.vector().set_local((nu / ssv_) ** 0.5)
        time_scale.vector().apply("insert")
        time_scale_avg.vector().axpy(1, time_scale.vector())

        # Compute velocity scale
        velocity_scale.vector().set_local((ssv_ * nu) ** (1. / 4))
        velocity_scale.vector().apply("insert")
        velocity_scale_avg.vector().axpy(1, velocity_scale.vector())

        # Compute both kinetic energy and turbulent kinetic energy
        assign(u0, u.sub(0))
        assign(u1, u.sub(1))
        assign(u2, u.sub(2))
        
        kinetic_energy.vector().set_local(
            0.5 * (u0.vector().get_local() ** 2 + u1.vector().get_local() ** 2 + u2.vector().get_local() ** 2))
        kinetic_energy.vector().apply("insert")
        kinetic_energy_avg.vector().axpy(1, kinetic_energy.vector())

        assign(u0_prime, u_prime.sub(0))
        assign(u1_prime, u_prime.sub(1))
        assign(u2_prime, u_prime.sub(2))

        turbulent_kinetic_engergy.vector().set_local(
            0.5 * (u0_prime.vector().get_local() ** 2 + u1_prime.vector().get_local() ** 2 + u2_prime.vector().get_local() ** 2))
        turbulent_kinetic_engergy.vector().apply("insert")
        turbulent_kinetic_engergy_avg.vector().axpy(1, turbulent_kinetic_engergy.vector())

    # Get avg
    N = len(dataset_names)
    l_pluss_avg.vector()[:] = l_pluss_avg.vector()[:] / N
    t_pluss_avg.vector()[:] = t_pluss_avg.vector()[:] / N
    length_scale_avg.vector()[:] = length_scale_avg.vector()[:] / N
    time_scale_avg.vector()[:] = time_scale_avg.vector()[:] / N
    velocity_scale_avg.vector()[:] = velocity_scale_avg.vector()[:] / N
    CFL_avg.vector()[:] = CFL_avg.vector()[:] / N
    dissipation_avg.vector()[:] = dissipation_avg.vector()[:] / N
    kinetic_energy_avg.vector()[:] = kinetic_energy_avg.vector()[:] / N
    turbulent_kinetic_engergy_avg.vector()[:] = turbulent_kinetic_engergy_avg.vector()[:] / N
    turbulent_dissipation_avg.vector()[:] = turbulent_dissipation_avg.vector()[:] / N

    # Store average data
    metrices["CFL"].write(CFL_avg)
    metrices["l_pluss"].write(l_pluss_avg)
    metrices["t_pluss"].write(t_pluss_avg)
    metrices["length_scale"].write(length_scale_avg)
    metrices["time_scale"].write(time_scale_avg)
    metrices["velocity_scale"].write(velocity_scale_avg)
    metrices["dissipation"].write(dissipation_avg)
    metrices["kinetic_energy"].write(kinetic_energy_avg)
    metrices["turbulent_kinetic_engergy"].write(turbulent_kinetic_engergy_avg)
    metrices["turbulent_dissipation"].write(turbulent_dissipation_avg)
    metrices["u_mean"].write(u_mean)

    # Print info
    flow_metrics = [("dx", dl), ("l+", l_pluss_avg), ("t+", t_pluss_avg), ("Length scale", length_scale_avg),
         ("Time scale", time_scale_avg), ("Velocity scale", velocity_scale), ("CFL", CFL_avg), ("SSV", ssv_avg),
         ("dissipation", dissipation), ("turbulent dissipation", turbulent_dissipation),
         ("turbulent_kinetic_engergy", turbulent_kinetic_engergy), ("kinetic_energy", kinetic_energy)]

    for metric in flow_metrics:
        sum_ = MPI.sum(MPI.comm_world, np.sum(metric[1].vector().get_local()))
        num = MPI.sum(MPI.comm_world, metric[1].vector().get_local().shape[0])
        mean = sum_ / num
        max_ = MPI.max(MPI.comm_world, metric[1].vector().get_local().max())
        min_ = MPI.min(MPI.comm_world, metric[1].vector().get_local().min())

        if MPI.rank(MPI.comm_world) == 0:
            print(metric[0], "mean:", mean)
            print(metric[0], "max:", max_)
            print(metric[0], "min:", min_)


if __name__ == '__main__':
    folder, nu, dt, velocity_degree = read_command_line()
    compute_flow_metrics(folder, nu, dt)
