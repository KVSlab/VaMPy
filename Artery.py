import pickle
from os import path, makedirs, getcwd

import numpy as np
import dolfin
from fenicstools import Probes
from oasis.problems.NSfracStep import *

from Womersley import make_womersley_bcs, compute_boundary_geometry_acrn

set_log_level(50)


def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, 'params.dat'), 'r')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
    else:
        # Override some problem specific parameters
        # parameters are in mm and ms
        NS_parameters.update(
            nu = 3.3018e-3,
            T  = 951*2,
            dt = 0.0951, # 10 000 steps per cycle
            velocity_degree = 1,
            folder = "results_folder",
            save_step = 10000,
            checkpoint = 500,
            no_of_cycles = 2,
            mesh_path=commandline_kwargs["mesh_path"],
            id_in = [],
            id_out = [],
            area_ratio = [],
            dump_stats = 500,
            store_data = 5,
            compute_flux = 5,
            plot_interval = 10e10,
            print_intermediate_info = 100,
            use_krylov_solvers = True,
            krylov_solvers = dict(monitor_convergence=False)
        )

    caseName = NS_parameters["mesh_path"].split(".")[0]
    NS_parameters["folder"] = path.join(NS_parameters["folder"], caseName)


def mesh(mesh_path, **NS_namespace):
    # Read mesh
    mesh_folder = path.join(path.dirname(path.abspath(__file__)), mesh_path)
    return Mesh(mesh_folder)


def create_bcs(u_, t, NS_expressions, V, Q, area_ratio, mesh, folder, mesh_path, nu,
               id_in, id_out, velocity_degree, pressure_degree, no_of_cycles,
               T, **NS_namespace):
    # Mesh function
    fd = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())

    # Extract flow split ratios
    info = open(path.join(path.dirname(path.abspath(__file__)), mesh_path.split(".")[0] \
                           + ".txt"), "r").readlines()
    for line in info:
        if "inlet_area" in line:
            inlet_area = float(line.split(":")[-1])
        elif "idFileLine" in line:
            _, _, id_in_, id_out_, Q_mean = line.split()
            id_in.append(int(id_in_))
            id_out[:] = [int(p) for p in id_out_.split(",")]
            Q_mean = float(Q_mean)
        elif "areaRatioLine" in line:
            area_ratio[:] = [float(p) for p in line.split()[-1].split(",")]

    # Womersley boundary condition at inlet
    t_values = np.linspace(0, T)
    t_values, Q_ = np.load(path.join(path.dirname(path.abspath(__file__)), "ICA_values"))
    Q_values = Q_mean * Q_
    t_values *= 1000
    tmp_a, tmp_c, tmp_r, tmp_n = compute_boundary_geometry_acrn(mesh, id_in[0], fd)
    inlet = make_womersley_bcs(t_values, Q_values, mesh, nu, tmp_a,
                                  tmp_c, tmp_r, tmp_n, V.ufl_element())
    NS_expressions["inlet"] = inlet

    # Set start time equal to t_0
    for uc in inlet:
        uc.set_t(t)

    # Create pressure boundary condition
    area_out = []
    for i, ind in enumerate(id_out):
        dsi = ds(ind, domain=mesh, subdomain_data=fd)
        area_out.append(assemble(Constant(1.0, name="one")*dsi))

    bc_p = []
    print("Initial pressure:")
    for i, ID in enumerate(id_out):
        p_initial = area_out[i] / sum(area_out)
        outflow = Expression("p", p=p_initial, degree=pressure_degree)
        bc = DirichletBC(Q, outflow, fd, ID)
        bc_p.append(bc)
        NS_expressions[ID] = outflow
        print(ID,  p_initial)

    # Noslip condition at wall
    wall = Constant(0.0)

    # Create Boundary conditions for the velocity
    bc_wall = DirichletBC(V, wall, fd, 0)
    bc_inlet = [DirichletBC(V, inlet[i], fd, id_in[0]) for i in range(3)]

    # Return boundary conditions in dictionary
    return dict(u0=[bc_inlet[0], bc_wall],
                u1=[bc_inlet[1], bc_wall],
                u2=[bc_inlet[2], bc_wall],
                p=bc_p)


def get_file_paths(folder):
    if MPI.rank(MPI.comm_world) == 0:
        counter = 1
        to_check = path.join(folder, "data", "%d")
        while path.isdir(to_check % counter):
            counter += 1

        if counter > 1:
            counter -= 1
        if not path.exists(path.join(to_check % counter, "VTK")):
            makedirs(path.join(to_check % counter, "VTK"))
    else:
        counter = 0

    counter = int(MPI.max(MPI.comm_world, counter))

    common_path = path.join(folder, "data", str(counter), "VTK")
    #common_path1 = path.join(folder, "data", str(counter), "hdf5")
    file_u = [path.join(common_path, "u%d.xdmf" % i) for i in range(3)]
    file_p = path.join(common_path, "p.xdmf")
    file_nu = path.join(common_path, "nut.xdmf")
    file_mesh = path.join(common_path, "mesh.xdmf")
    file_u_mean = [path.join(common_path, "u%d_mean.xdmf" % i) for i in range(3)]
    file_u_mean_vec = path.join(common_path, "u_mean.xdmf")
    files = {"u": file_u, "p": file_p, "u_mean": file_u_mean, "nut": file_nu,
             "u_mean_vec": file_u_mean_vec, "mesh": file_mesh}

    #common_path1 = path.join(folder, "data", str(counter), "hdf5")
    #mesh_path0 = path.join(common_path1, "mesh.h5")
    #p_path = path.join(common_path1, "p.h5")
    #u0_path = path.join(common_path1, "u0.h5")
    #u1_path = path.join(common_path1, "u1.h5")
    #u2_path = path.join(common_path1, "u2.h5")
    #ile2 = [mesh_path0,p_path,u0_path,u1_path,u2_path]

    return files


def pre_solve_hook(mesh, V, Q, newfolder, folder, u_, mesh_path,
                   restart_folder, velocity_degree, **NS_namespace):

    # Create point for evaluation
    fd = MeshFunction("size_t", mesh, 2, mesh.domains())
    n = FacetNormal(mesh)
    eval_dict = {}
    rel_path = path.join(path.dirname(path.abspath(__file__)), mesh_path.split(".")[0] + "_probe_point")
    probe_points = np.load(rel_path, encoding='latin1', fix_imports=True, allow_pickle=True) #Original Code
    #probe_points = np.loadtxt(rel_path, delimiter=',')
    #probe_points = np.load(rel_path, allow_pickle=False)

     #Store points file in checkpoint
    if MPI.rank(MPI.comm_world) == 0:
        probe_points.dump(path.join(newfolder, "Checkpoint", "points"))

    eval_dict["centerline_u_x_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_y_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_z_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_p_probes"] = Probes(probe_points.flatten(), Q)

    if restart_folder is None:
        # Get files to store results
        files = get_file_paths(folder)
        NS_parameters.update(dict(files=files))
    else:
        files = NS_namespace["files"]

    with XDMFFile(MPI.comm_world, files["mesh"]) as mesh_file:
        mesh_file.write(mesh)

    writer = {}
    for key, value in files.items():
        if isinstance(value, list):
            values = value
            for i, value in enumerate(values):
                writer[key + str(i)] = dolfin.XDMFFile(MPI.comm_world, value)
                writer[key + str(i)].parameters["flush_output"] = True
                writer[key + str(i)].parameters["functions_share_mesh"] = True
                writer[key + str(i)].parameters["rewrite_function_mesh"] = False
        else:
            writer[key] = dolfin.XDMFFile(MPI.comm_world, value)
            writer[key].parameters["flush_output"] = True
            writer[key].parameters["functions_share_mesh"] = True
            writer[key].parameters["rewrite_function_mesh"] = False

    return dict(eval_dict=eval_dict, fd=fd, n=n, writer=writer)


def beta(err, p):
    if p < 0:
        if err >= 0.1:
            return 0.5
        else:
            return 1 -5*err**2
    else:
        if err >= 0.1:
            return 1.5
        else:
            return 1  + 5*err**2


def w(P):
    return 1 / ( 1 + 20*abs(P))


def temporal_hook(u_, p_, p, Q, mesh, tstep, compute_flux, dump_stats, eval_dict,
                  newfolder, id_in, id_out, fd, n, store_data, NS_expressions,
                  area_ratio, t, writer, **NS_namespace):

    # Update boundary condition
    for uc in NS_expressions["inlet"]:
        uc.set_t(t)

    # Compute flux and update pressure condition
    if tstep > 2 and tstep % 1 == 0:

        Q_in = abs(assemble(dot(u_, n)*ds(id_in[0], domain=mesh, subdomain_data=fd)))
        Q_outs =  []
        Q_ideals = []
        for i, out_id in enumerate(id_out):
            Q_out = abs(assemble(dot(u_, n)*ds(out_id, domain=mesh, subdomain_data=fd)))
            Q_outs.append(Q_out)

            Q_ideal = area_ratio[i]*Q_in
            Q_ideals.append(Q_ideal)

            p_old = NS_expressions[out_id].p

            # Gin and Steinman et al., A Dual-Pressure Boundary Condition
            # for use in Simulations of Bifurcating Conduits
            R_optimal = area_ratio[i]
            R_actual = Q_out / Q_in

            M_err = abs(R_optimal / R_actual)
            R_err = abs(R_optimal - R_actual)

            if p_old < 0:
                E = 1 + R_err / R_optimal
            else:
                E = -1 * ( 1 + R_err / R_optimal )


            # 1) Linear update to converge first 100 tsteps of first cycle
            delta = (R_optimal - R_actual) / R_optimal
            if tstep < 100:
                h = 0.1
                if p_old > 1 and delta < 0:
                    NS_expressions[out_id].p  = p_old
                else:
                    NS_expressions[out_id].p  = p_old * ( 1 - delta*h)

            # 2) Dual pressure BC
            else:
                if p_old > 2 and delta < 0:
                    NS_expressions[out_id].p  = p_old
                else:
                    NS_expressions[out_id].p  = p_old * beta(R_err, p_old) * M_err ** E

    if MPI.rank(MPI.comm_world) == 0 and tstep % 10 == 0:
        print("="*10, tstep, "="*10)
        print("Sum of Q_out = {:0.4f} Q_in = {:0.4f}".format(sum(Q_outs), Q_in))
        for i, out_id in enumerate(id_out):
            print(("({:d}) New pressure {:0.4f}").format(out_id, NS_expressions[out_id].p))
        for i, out_id in enumerate(id_out):
            print(("({:d}) area ratio {:0.4f}, ideal: {:0.4f} actual:" + \
                  " {:0.4f}").format(out_id, area_ratio[i], Q_ideals[i], Q_outs[i]))
        print()

    # Sample velocity in points
    eval_dict["centerline_u_x_probes"](u_[0])
    eval_dict["centerline_u_y_probes"](u_[1])
    eval_dict["centerline_u_z_probes"](u_[2])
    eval_dict["centerline_p_probes"](p_)

    # Store sampled velocity
    if tstep % dump_stats == 0:
        filepath = path.join(newfolder, "Stats")
        if MPI.rank(MPI.comm_world) == 0:
            if not path.exists(filepath):
                makedirs(filepath)

        arr_u_x = eval_dict["centerline_u_x_probes"].array()
        arr_u_y = eval_dict["centerline_u_y_probes"].array()
        arr_u_z = eval_dict["centerline_u_z_probes"].array()
        arr_p = eval_dict["centerline_p_probes"].array()

        # Dump stats
        if MPI.rank(MPI.comm_world) == 0:
            num = eval_dict["centerline_u_x_probes"].number_of_evaluations()
            pp = (path.join(filepath, "u_x_%s.probes" % str(tstep)))
            arr_u_x.dump(path.join(filepath, "u_x_%s.probes" % str(tstep)))
            arr_u_y.dump(path.join(filepath, "u_y_%s.probes" % str(tstep)))
            arr_u_z.dump(path.join(filepath, "u_z_%s.probes" % str(tstep)))
            arr_p.dump(path.join(filepath, "p_%s.probes" % str(tstep)))

        # Clear stats
        MPI.barrier(MPI.comm_world)
        eval_dict["centerline_u_x_probes"].clear()
        eval_dict["centerline_u_y_probes"].clear()
        eval_dict["centerline_u_z_probes"].clear()
        eval_dict["centerline_p_probes"].clear()

    # Save velocity and pressure
    #common_path1 = path.join(folder, "data", str(counter), "hdf5")
    #mesh_path = path.join(common_path1, "mesh.h5")
    #p_path = path.join(common_path1, "p.h5")
    #u0_path = path.join(common_path1, "u0.h5")
    #u1_path = path.join(common_path1, "u1.h5")
    #u2_path = path.join(common_path1, "u2.h5")
    # Save velocity and pressure
    if tstep % store_data == 0:
        # Name functions
        u_[0].rename("u0", "velocity-x")
        u_[1].rename("u1", "velocity-y")
        u_[2].rename("u2", "velocity-z")
        p_.rename("p", "pressure")

        folder = "results_folder" #folder is a local variable

        if MPI.rank(MPI.comm_world) == 0:
            counter = 1
            to_check = path.join(folder, "data", "%d")
            while path.isdir(to_check % counter):
                counter += 1

            if counter > 1:
                counter -= 1
            if not path.exists(path.join(to_check % counter, "hdf5")):
                makedirs(path.join(to_check % counter, "hdf5"))
        else:
            counter = 0

        counter = int(MPI.max(MPI.comm_world, counter))

        common_path1 = path.join(folder, "data", str(counter), "hdf5")
        mesh_path0 = path.join(common_path1, "mesh2.h5")
        p_path = path.join(common_path1, "p.h5")
        u0_path = path.join(common_path1, "u0.h5")
        u1_path = path.join(common_path1, "u1.h5")
        u2_path = path.join(common_path1, "u2.h5")

        if tstep == 10000:
            viz_mesh = HDF5File(MPI.comm_world, mesh_path0, "w")
            viz_mesh.write(mesh,"mesh")
            viz_mesh.close()
            viz_p = HDF5File(MPI.comm_world, p_path, "w")
            viz_p.write(p_,"/pressure",tstep)
            viz_p.close()
            viz_u0 = HDF5File(MPI.comm_world, u0_path, "w")
            viz_u0.write(u_[0],"/velocity",tstep)
            viz_u0.close()
            viz_u1 = HDF5File(MPI.comm_world, u1_path, "w")
            viz_u1.write(u_[1],"/velocity", tstep)
            viz_u1.close()
            viz_u2 = HDF5File(MPI.comm_world, u2_path, "w")
            viz_u2.write(u_[2],"/velocity", tstep)
            viz_u2.close()
        elif tstep>10000 and tstep % store_data == 0:
            viz_p = HDF5File(MPI.comm_world, p_path, "a")
            viz_p.write(p_,"/pressure",tstep)
            viz_p.close()
            viz_u0 = HDF5File(MPI.comm_world, u0_path, "a")
            viz_u0.write(u_[0],"/velocity", tstep)
            viz_u0.close()
            viz_u1 = HDF5File(MPI.comm_world, u1_path, "a")
            viz_u1.write(u_[1],"/velocity", tstep)
            viz_u1.close()
            viz_u2 = HDF5File(MPI.comm_world, u2_path, "a")
            viz_u2.write(u_[2],"/velocity", tstep)
            viz_u2.close()

#    save_step=100 #treating save_step as a local variable
#    if tstep>=10000 and tstep % save_step == 0:
        # Name functions
#        u_[0].rename("u0", "velocity-x")
#        u_[1].rename("u1", "velocity-y")
#        u_[2].rename("u2", "velocity-z")
#        p_.rename("p", "pressure")

        # Store files
#        components = {"u0": u_[0], "u1": u_[1], "u2": u_[2], "p": p_}
#        for key in components.keys():
#            if tstep == save_step:
#                writer[key].write_checkpoint(components[key], components[key].name(), tstep,
#                                             XDMFFile.Encoding.HDF5, append=False)
#            else:
#                writer[key].write_checkpoint(components[key], components[key].name(), tstep,
#                                             XDMFFile.Encoding.HDF5, append=True)
