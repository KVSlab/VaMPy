from argparse import ArgumentParser
from os import path
from time import time

import numpy as np
from dolfin import *

set_log_active(False)


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--case', type=str, default="/results_folder/1/VTK", help="Path to simulation results",
                        metavar="PATH")
    parser.add_argument('--nu', type=float, default=1e-6, help="Kinematic viscosity used in simulation")
    parser.add_argument('--dt', type=float, default=0.0951, help="Time step of simulation")

    args = parser.parse_args()

    return args.case, args.nu, args.dt


def get_dataset_names(f, num_files=3000000, step=1, start=1, print_info=True, velocity_filename="velocity%s"):
    check = True

    # Find start file
    t0 = time()
    while check:
        if f.has_dataset(velocity_filename % start):
            check = False
            start -= step
        start += step

    # Get names
    names = []
    for i in range(num_files):
        step = 1
        index = start + i * step
        if f.has_dataset(velocity_filename % index):
            names.append(velocity_filename % index)

    t1 = time()

    # Print info
    if MPI.rank(MPI.comm_world) == 0 and print_info:
        print("")
        print("=" * 6 + " Timesteps to average over " + "=" * 6)
        print("Length on data set names:", len(names))
        print("Start index:", start)
        print("Wanted num files:", num_files)
        print("Step between files:", step)
        print("Time used:", t1 - t0, "s")

    return names


def rate_of_strain(ssv, u, v, mesh, h, nu):
    epsilon = 0.5 * (grad(u) + grad(u).T)
    f = sqrt(inner(epsilon, epsilon))
    x = assemble(inner(f, v) / h * dx(mesh))
    ssv.vector().set_local(x.array())
    ssv.vector().apply("insert")

    return ssv


def rate_of_dissipation(ssv, u, v, mesh, h, nu):
    epsilon = 0.5 * (grad(u) + grad(u).T)
    f = 2 * nu * inner(epsilon, epsilon)
    x = assemble(inner(f, v) / h * dx(mesh))
    ssv.vector().set_local(x.array())
    ssv.vector().apply("insert")

    return ssv


def save_h5(filename, field, timestep):
    f_xml = File(filename + timestep + ".xml.gz")
    f_xml << field
    del f_xml

    f_pvd = File(filename + ".pvd")
    f_pvd << field
    del f_pvd


def compute_info(folder, nu, dt):
    file_path_x = path.join(folder, "u0.h5")
    file_path_y = path.join(folder, "u1.h5")
    file_path_z = path.join(folder, "u2.h5")

    f_0 = HDF5File(MPI.comm_world, file_path_x, "r")
    f_1 = HDF5File(MPI.comm_world, file_path_y, "r")
    f_2 = HDF5File(MPI.comm_world, file_path_z, "r")

    # Get names of data to extract
    start = 0
    if MPI.rank(MPI.comm_world()) == 0:
        print
        "The post processing starts from", start
    dataset_names = get_dataset_names(f_0, start=start)
    # dataset_names = get_dataset_names(f_,  start=start) #beforeV

    # Get mesh information
    case = folder.split("_")[0]
    mesh = Mesh()
    f_0.read(mesh, "Mesh", False)

    # Function space
    u_order = 1
    # u_order = int(folder.split("_")[1].split("P")[1])
    DG = FunctionSpace(mesh, 'DG', 0)
    V = VectorFunctionSpace(mesh, "CG", u_order)
    Vv = FunctionSpace(mesh, "CG", u_order)

    # Mesh info
    h = CellDiameter(mesh)
    dl = project(h, DG)

    # Functions for storing values
    v = TestFunction(DG)
    u = Function(V)
    u_mean = Function(V)
    u_prime = Function(V)
    ssv = Function(DG)
    ssv_avg = Function(DG)
    dissipation = Function(DG)
    dissipation_avg = Function(DG)

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
    fullname = file_path_x.replace("u0.h5", "%s.xdmf")
    fullname = fullname.replace("VTK", "mesh_info_testTKE")
    var_name = ["u_mean", "l_pluss", "t_pluss", "CFL", "ssv", "length_scale", "time_scale", "velocity_scale", "u_mag",
                "dl", "dissipation", "kinetic_energy", "turbulent_kinetic_engergy", "turbulent_dissipation", "u_prime",
                "u_viz"]
    # TODO: add WSS, RRT, OSI, TWSSG, WSSG
    f = {}
    for vn in var_name:
        if MPI.rank(mpi_comm_world()) == 0:
            print
            fullname % vn
        f[vn] = XDMFFile(mpi_comm_world(), fullname % vn)
        f[vn].parameters["rewrite_function_mesh"] = False
        f[vn].parameters["flush_output"] = True

    f["dl"] << dl
    del f["dl"]

    # Get u mean
    u_mean_file_path = file_path_x.replace("u0.h5", "u%d_mean.h5")
    for i in range(3):
        tmp_file = HDF5File(mpi_comm_world(), u_mean_file_path % i, "r")
        tmp_file.read(u0, "u_mean/avg")
        tmp_file.close()
        assign(u_mean.sub(i), u0)

    for data in dataset_names:
        if MPI.rank(mpi_comm_world()) == 0:
            print
            data

        # Timestep and velocity
        timestep = data.split("velocity")[-1]
        f_0.read(u0, data)
        f_1.read(u1, data)
        f_2.read(u2, data)
        assign(u.sub(0), u0)
        assign(u.sub(1), u1)
        assign(u.sub(2), u2)
        # f["u_viz"].write(u)

        # Compute CFL
        u_mag = project(sqrt(inner(u, u)), DG)
        # f["u_mag"] << u_mag

        CFL.vector().set_local(u_mag.vector().array() / dl.vector().array() * dt)
        CFL.vector().apply("insert")
        CFL_avg.vector().axpy(1, CFL.vector())
        # f["CFL"] << CFL

        # Compute rate-of-strain
        rate_of_strain(ssv, u, v, mesh, h, nu)
        ssv_avg.vector().axpy(1, ssv.vector())
        # f["ssv"] << ssv

        # Compute l+ and t+
        u_star = np.sqrt(ssv.vector().array() * nu)
        l_pluss.vector().set_local(u_star * dl.vector().array() / nu)
        l_pluss.vector().apply("insert")
        l_pluss_avg.vector().axpy(1, l_pluss.vector())
        # f["l_pluss"] << l_pluss

        t_pluss.vector().set_local(nu / u_star ** 2)
        t_pluss.vector().apply("insert")
        t_pluss_avg.vector().axpy(1, t_pluss.vector())
        # f["t_pluss"] << t_pluss

        # Compute Kolmogorov
        rate_of_dissipation(dissipation, u, v, mesh, h, nu)
        ssv_ = dissipation.vector().array()
        # f["dissipation"] << dissipation
        dissipation_avg.vector().axpy(1, dissipation.vector())

        u_prime.vector().set_local(u.vector().array() - u_mean.vector().array())
        u_prime.vector().apply("insert")
        # f["u_prime"] << u_prime

        rate_of_dissipation(turbulent_dissipation, u_prime, v, mesh, h, nu)
        # f["turbulent_dissipation"] << turbulent_dissipation
        turbulent_dissipation_avg.vector().axpy(1, turbulent_dissipation.vector())

        length_scale.vector().set_local((nu ** 3 / ssv_) ** (1. / 4))
        length_scale.vector().apply("insert")
        length_scale_avg.vector().axpy(1, length_scale.vector())
        # f["length_scale"] << length_scale

        time_scale.vector().set_local((nu / ssv_) ** 0.5)
        time_scale.vector().apply("insert")
        time_scale_avg.vector().axpy(1, time_scale.vector())
        # f["time_scale"] << time_scale

        velocity_scale.vector().set_local((ssv_ * nu) ** (1. / 4))
        velocity_scale.vector().apply("insert")
        velocity_scale_avg.vector().axpy(1, velocity_scale.vector())
        # f["velocity_scale"] << velocity_scale

        # Energy
        assign(u0_prime, u_prime.sub(0))
        assign(u1_prime, u_prime.sub(1))
        assign(u2_prime, u_prime.sub(2))

        kinetic_energy.vector().set_local(
            0.5 * (u0.vector().array() ** 2 + u1.vector().array() ** 2 + u2.vector().array() ** 2))
        kinetic_energy.vector().apply("insert")
        kinetic_energy_avg.vector().axpy(1, kinetic_energy.vector())
        # f["kinetic_energy"] << kinetic_energy

        turbulent_kinetic_engergy.vector().set_local(
            0.5 * (u0_prime.vector().array() ** 2 + u1_prime.vector().array() ** 2 + u2_prime.vector().array() ** 2))
        turbulent_kinetic_engergy.vector().apply("insert")
        turbulent_kinetic_engergy_avg.vector().axpy(1, turbulent_kinetic_engergy.vector())
        # f["turbulent_kinetic_engergy"] << turbulent_kinetic_engergy

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
    f["CFL"] << CFL_avg
    f["l_pluss"] << l_pluss_avg
    f["t_pluss"] << t_pluss_avg
    f["length_scale"] << length_scale_avg
    f["time_scale"] << time_scale_avg
    f["velocity_scale"] << velocity_scale_avg
    f["dissipation"] << dissipation_avg
    f["kinetic_energy"] << kinetic_energy_avg
    f["turbulent_kinetic_engergy"] << turbulent_kinetic_engergy_avg
    f["turbulent_dissipation"] << turbulent_dissipation_avg
    f["u_mean"] << u_mean

    # print info
    a = [("dx", dl), ("l+", l_pluss_avg), ("t+", t_pluss_avg), ("Length scale", length_scale_avg),
         ("Time scale", time_scale_avg), ("Velocity scale", velocity_scale), ("CFL", CFL_avg),
         ("SSV", ssv_avg), ("dissipation", dissipation),
         ("turbulent dissipation", turbulent_dissipation),
         ("turbulent_kinetic_engergy", turbulent_kinetic_engergy),
         ("kinetic_energy", kinetic_energy)]

    for item in a:
        sum_ = MPI.sum(mpi_comm_world(), np.sum(item[1].vector().array()))
        num = MPI.sum(mpi_comm_world(), item[1].vector().array().shape[0])
        mean = sum_ / num
        max_ = MPI.max(mpi_comm_world(), item[1].vector().array().max())
        min_ = MPI.min(mpi_comm_world(), item[1].vector().array().min())

        if MPI.rank(mpi_comm_world()) == 0:
            print
            item[0], "mean:", mean
            print
            item[0], "max:", max_
            print
            item[0], "min:", min_


if __name__ == '__main__':
    folder, nu, dt = read_command_line()
    compute_info(folder, nu, dt)
