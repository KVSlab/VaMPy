import json
import pickle
from os import path, makedirs

import numpy as np
from Womersley import make_womersley_bcs, compute_boundary_geometry_acrn
from fenicstools import Probes

from oasis.problems.NSfracStep import *

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
            # Fluid parameters
            nu=3.3018e-3,  # Viscosity
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            area_ratio=[],
            # Simulation parameters
            T=951 * 2,  # Run simulation for 2 cardiac cycles
            dt=0.0951,  # 10 000 steps per cycle
            no_of_cycles=2,
            dump_stats=100,
            store_data=5,
            store_data_tstep=10000,  # Start storing data at 2nd cycle
            save_step=200,
            checkpoint=500,
            print_intermediate_info=100,
            folder="results_artery",
            mesh_path=commandline_kwargs["mesh_path"],
            # Solver parameters
            velocity_degree=1,
            use_krylov_solvers=True,
            krylov_solvers=dict(monitor_convergence=False)
        )

    mesh_file = NS_parameters["mesh_path"].split("/")[-1]
    case_name = mesh_file.split(".")[0]
    NS_parameters["folder"] = path.join(NS_parameters["folder"], case_name)


def mesh(mesh_path, **NS_namespace):
    # Read mesh
    return Mesh(mesh_path)


def create_bcs(t, NS_expressions, V, Q, area_ratio, mesh, mesh_path, nu, id_in, id_out, pressure_degree, **NS_namespace):
    # Mesh function
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())

    # Read case parameters
    info_path = mesh_path.split(".")[0] + ".json"
    with open(info_path) as f:
        info = json.load(f)

    # Extract flow split ratios and inlet/outlet IDs
    id_info = info['idFileLine'].split()
    id_in.append(int(id_info[1]))
    id_out[:] = [int(p) for p in id_info[2].split(",")]
    Q_mean = float(id_info[3])
    area_ratio[:] = [float(p) for p in info['areaRatioLine'].split()[-1].split(",")]

    # Womersley boundary condition at inlet
    t_values, Q_ = np.load(path.join(path.dirname(path.abspath(__file__)), "ICA_values"))
    Q_values = Q_mean * Q_
    t_values *= 1000
    tmp_a, tmp_c, tmp_r, tmp_n = compute_boundary_geometry_acrn(mesh, id_in[0], boundary)
    inlet = make_womersley_bcs(t_values, Q_values, mesh, nu, tmp_a, tmp_c, tmp_r, tmp_n, V.ufl_element())
    NS_expressions["inlet"] = inlet

    # Set start time equal to t_0
    for uc in inlet:
        uc.set_t(t)

    # Create pressure boundary condition
    area_out = []
    for i, ind in enumerate(id_out):
        dsi = ds(ind, domain=mesh, subdomain_data=boundary)
        area_out.append(assemble(Constant(1.0, name="one") * dsi))

    bc_p = []
    print("Initial pressure:")
    for i, ID in enumerate(id_out):
        p_initial = area_out[i] / sum(area_out)
        outflow = Expression("p", p=p_initial, degree=pressure_degree)
        bc = DirichletBC(Q, outflow, boundary, ID)
        bc_p.append(bc)
        NS_expressions[ID] = outflow
        print(ID, p_initial)

    # No slip condition at wall
    wall = Constant(0.0)

    # Create Boundary conditions for the velocity
    bc_wall = DirichletBC(V, wall, boundary, 0)
    bc_inlet = [DirichletBC(V, inlet[i], boundary, id_in[0]) for i in range(3)]

    # Return boundary conditions in dictionary
    return dict(u0=[bc_inlet[0], bc_wall],
                u1=[bc_inlet[1], bc_wall],
                u2=[bc_inlet[2], bc_wall],
                p=bc_p)


def get_file_paths(folder):
    # Create folder where data and solutions (velocity, mesh, pressure) is stored
    common_path = path.join(folder, "VTK")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path):
            makedirs(common_path)

    file_p = path.join(common_path, "p.h5")
    file_u = [path.join(common_path, "u{}.h5".format(i)) for i in range(3)]
    file_mesh = path.join(common_path, "mesh.h5")
    files = {"u": file_u, "p": file_p, "mesh": file_mesh}

    return files


def pre_solve_hook(mesh, V, Q, newfolder, mesh_path, restart_folder, **NS_namespace):
    # Create point for evaluation
    boundary = MeshFunction("size_t", mesh, 2, mesh.domains())
    n = FacetNormal(mesh)
    eval_dict = {}
    rel_path = mesh_path.split(".")[0] + "_probe_point"
    probe_points = np.load(rel_path, encoding='latin1', fix_imports=True, allow_pickle=True)

    # Store points file in checkpoint
    if MPI.rank(MPI.comm_world) == 0:
        probe_points.dump(path.join(newfolder, "Checkpoint", "points"))

    eval_dict["centerline_u_x_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_y_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_z_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_p_probes"] = Probes(probe_points.flatten(), Q)

    if restart_folder is None:
        # Get files to store results
        files = get_file_paths(newfolder)
        NS_parameters.update(dict(files=files))
    else:
        files = NS_namespace["files"]

    # Save mesh as HDF5 file
    with HDF5File(MPI.comm_world, files["mesh"], "w") as mesh_file:
        mesh_file.write(mesh, "mesh")

    return dict(eval_dict=eval_dict, boundary=boundary, n=n)


def temporal_hook(u_, p_, mesh, tstep, dump_stats, eval_dict, newfolder, id_in, id_out, boundary, n, store_data,
                  NS_parameters, NS_expressions, area_ratio, t, store_data_tstep, **NS_namespace):
    # Update boundary condition
    for uc in NS_expressions["inlet"]:
        uc.set_t(t)

    # Compute flux and update pressure condition
    if tstep > 2 and tstep % 1 == 0:
        Q_ideals, Q_in, Q_outs = update_pressure_condition(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n,
                                                           tstep, u_)

    if MPI.rank(MPI.comm_world) == 0 and tstep % 10 == 0:
        print("=" * 10, tstep, "=" * 10)
        print("Sum of Q_out = {:0.4f} Q_in = {:0.4f}".format(sum(Q_outs), Q_in))
        for i, out_id in enumerate(id_out):
            print(("({:d}) New pressure {:0.4f}").format(out_id, NS_expressions[out_id].p))
        for i, out_id in enumerate(id_out):
            print(("({:d}) area ratio {:0.4f}, ideal: {:0.4f} actual:" +
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
    if tstep % store_data == 0 and tstep >= store_data_tstep:
        # Name functions
        u_[0].rename("u0", "velocity-x")
        u_[1].rename("u1", "velocity-y")
        u_[2].rename("u2", "velocity-z")
        p_.rename("p", "pressure")

        # Get save paths
        files = NS_parameters['files']
        file_mode = "w" if tstep == store_data_tstep else "a"
        p_path = files['p']

        # Save
        viz_p = HDF5File(MPI.comm_world, p_path, file_mode=file_mode)
        viz_p.write(p_, "/pressure", tstep)
        viz_p.close()

        for i in range(3):
            u_path = files['u'][i]
            viz_u = HDF5File(MPI.comm_world, u_path, file_mode=file_mode)
            viz_u.write(u_[i], "/velocity", tstep)
            viz_u.close()


def beta(err, p):
    """
    Adjusted choice of beta from
    Gin and Steinman et al., A Dual-Pressure Boundary Condition doi:10.1115/1.1504446 
    Ramped up to desired value if flow rate error (err) increases
    Args:
        err (float): Flow split error
        p (float): Pressure value
    Returns:
        beta (float): Variable factor in flow split method
    """
    if p < 0:
        if err >= 0.1:
            return 0.5
        else:
            return 1 - 5 * err ** 2
    else:
        if err >= 0.1:
            return 1.5
        else:
            return 1 + 5 * err ** 2


def update_pressure_condition(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n, tstep, u_):
    """
    Use Gin and Steinman et al., A Dual-Pressure Boundary Condition
    for use in Simulations of Bifurcating Conduits
    as pressure condition
    """
    Q_in = abs(assemble(dot(u_, n) * ds(id_in[0], domain=mesh, subdomain_data=boundary)))
    Q_outs = []
    Q_ideals = []
    for i, out_id in enumerate(id_out):
        Q_out = abs(assemble(dot(u_, n) * ds(out_id, domain=mesh, subdomain_data=boundary)))
        Q_outs.append(Q_out)

        Q_ideal = area_ratio[i] * Q_in
        Q_ideals.append(Q_ideal)

        p_old = NS_expressions[out_id].p

        R_optimal = area_ratio[i]
        R_actual = Q_out / Q_in

        M_err = abs(R_optimal / R_actual)
        R_err = abs(R_optimal - R_actual)

        if p_old < 0:
            E = 1 + R_err / R_optimal
        else:
            E = -1 * (1 + R_err / R_optimal)

        # 1) Linear update to converge first 100 tsteps of first cycle
        delta = (R_optimal - R_actual) / R_optimal
        if tstep < 100:
            h = 0.1
            if p_old > 1 and delta < 0:
                NS_expressions[out_id].p = p_old
            else:
                NS_expressions[out_id].p = p_old * (1 - delta * h)

        # 2) Dual pressure BC
        else:
            if p_old > 2 and delta < 0:
                NS_expressions[out_id].p = p_old
            else:
                NS_expressions[out_id].p = p_old * beta(R_err, p_old) * M_err ** E

    return Q_ideals, Q_in, Q_outs