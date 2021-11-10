import json
import pickle
from os import path, makedirs
from pprint import pprint

import numpy as np
from fenicstools import Probes
from oasis.problems.NSfracStep import *

from Womersley import make_womersley_bcs, compute_boundary_geometry_acrn

"""
Problem file for running CFD simulation in arterial models consisting of one inlet, and two or more outlets.
A Womersley velocity profile is applied at the inlet, and a flow split pressure condition is applied at the outlets,
following [1]. Flow rate for the inlet condition, and flow split values for the outlets are computed from the 
pre-processing script automatedPreProcessing.py. The simulation is run for two cycles (adjustable), but only the 
results/solutions from the second cycle are stored to avoid non-physiological effects from the first cycle.
One cardiac cycle is set to 0.951 s from [2], and scaled by a factor of 1000, hence all parameters are in [mm] or [ms].  

[1] Gin, Ron, Anthony G. Straatman, and David A. Steinman. "A dual-pressure boundary condition for use in simulations 
    of bifurcating conduits." J. Biomech. Eng. 124.5 (2002): 617-619. 
[2] Hoi, Yiemeng, et al. "Characterization of volumetric flow rate waveforms at the carotid bifurcations of older 
    adults." Physiological measurement 31.3 (2010): 291.
"""

set_log_level(50)


def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
    else:
        # Parameters are in mm and ms
        cardiac_cycle = 951
        number_of_cycles = 2

        NS_parameters.update(
            # Fluid parameters
            nu=3.3018e-3,  # Kinematic viscosity
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            area_ratio=[],
            area_inlet=[],
            # Simulation parameters
            cardiac_cycle=cardiac_cycle,  # Cardiac cycle [ms]
            T=cardiac_cycle * number_of_cycles,  # Simulation end time [ms]
            dt=0.0951,  # Time step size [ms]
            save_probe_frequency=100,  # Save frequency for sampling velocity & pressure at probes along the centerline
            save_solution_frequency=5,  # Save frequency for post processing
            save_solution_after_cycle=1,  # Store solution after 1 cardiac cycle
            # Oasis specific parameters
            checkpoint=500,  # Checkpoint frequency
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

    if MPI.rank(MPI.comm_world) == 0:
        print("=== Starting simulation for case: {} ===".format(case_name))
        print("Running with the following parameters:")
        pprint(NS_parameters)


def mesh(mesh_path, **NS_namespace):
    # Read mesh and print mesh information
    mesh = Mesh(mesh_path)
    print_mesh_information(mesh)

    return mesh


def print_mesh_information(mesh):
    comm = MPI.comm_world
    local_xmin = mesh.coordinates()[:, 0].min()
    local_xmax = mesh.coordinates()[:, 0].max()
    local_ymin = mesh.coordinates()[:, 1].min()
    local_ymax = mesh.coordinates()[:, 1].max()
    local_zmin = mesh.coordinates()[:, 2].min()
    local_zmax = mesh.coordinates()[:, 2].max()
    xmin = comm.gather(local_xmin, 0)
    xmax = comm.gather(local_xmax, 0)
    ymin = comm.gather(local_ymin, 0)
    ymax = comm.gather(local_ymax, 0)
    zmin = comm.gather(local_zmin, 0)
    zmax = comm.gather(local_zmax, 0)

    local_num_cells = mesh.num_cells()
    local_num_edges = mesh.num_edges()
    local_num_faces = mesh.num_faces()
    local_num_facets = mesh.num_facets()
    local_num_vertices = mesh.num_vertices()
    num_cells = comm.gather(local_num_cells, 0)
    num_edges = comm.gather(local_num_edges, 0)
    num_faces = comm.gather(local_num_faces, 0)
    num_facets = comm.gather(local_num_facets, 0)
    num_vertices = comm.gather(local_num_vertices, 0)
    volume = assemble(Constant(1) * dx(mesh))

    if MPI.rank(MPI.comm_world) == 0:
        print("=== Mesh information ===")
        print("X range: {} to {} (delta: {:.4f})".format(min(xmin), max(xmax), max(xmax) - min(xmin)))
        print("Y range: {} to {} (delta: {:.4f})".format(min(ymin), max(ymax), max(ymax) - min(ymin)))
        print("Z range: {} to {} (delta: {:.4f})".format(min(zmin), max(zmax), max(zmax) - min(zmin)))
        print("Number of cells: {}".format(sum(num_cells)))
        print("Number of cells per processor: {}".format(int(np.mean(num_cells))))
        print("Number of edges: {}".format(sum(num_edges)))
        print("Number of faces: {}".format(sum(num_faces)))
        print("Number of facets: {}".format(sum(num_facets)))
        print("Number of vertices: {}".format(sum(num_vertices)))
        print("Volume: {:.4f}".format(volume))
        print("Number of cells per volume: {:.4f}".format(sum(num_cells) / volume))


def create_bcs(t, NS_expressions, V, Q, area_ratio, area_inlet, mesh, mesh_path, nu, id_in, id_out, pressure_degree,
               **NS_namespace):
    # Mesh function
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())

    # Read case parameters
    parameters_file_path = mesh_path.split(".")[0] + ".json"
    with open(parameters_file_path) as f:
        parameters = json.load(f)

    # Extract flow split ratios and inlet/outlet IDs
    id_info = parameters['idFileLine'].split()
    id_in.append(int(id_info[1]))
    id_out[:] = [int(p) for p in id_info[2].split(",")]
    Q_mean = float(id_info[3])
    area_ratio[:] = [float(p) for p in parameters['areaRatioLine'].split()[-1].split(",")]
    area_inlet.append(float(parameters['inlet_area']))

    # Load normalized time and flow rate values
    t_values, Q_ = np.loadtxt(path.join(path.dirname(path.abspath(__file__)), "ICA_values")).T
    Q_values = Q_mean * Q_  # Specific flow rate * Flow wave form
    t_values *= 1000  # Scale to [ms]
    tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, id_in[0], boundary)

    # Create Womersley boundary condition at inlet
    inlet = make_womersley_bcs(t_values, Q_values, mesh, nu, tmp_area, tmp_center, tmp_radius, tmp_normal,
                               V.ufl_element())
    NS_expressions["inlet"] = inlet

    # Initialize inlet expressions with initial time
    for uc in inlet:
        uc.set_t(t)

    # Create pressure boundary condition
    area_out = []
    for i, ind in enumerate(id_out):
        dsi = ds(ind, domain=mesh, subdomain_data=boundary)
        area_out.append(assemble(Constant(1.0, name="one") * dsi))

    bc_p = []
    if MPI.rank(MPI.comm_world) == 0:
        print("=== Initial pressure and area fraction ===")
    for i, ID in enumerate(id_out):
        p_initial = area_out[i] / sum(area_out)
        outflow = Expression("p", p=p_initial, degree=pressure_degree)
        bc = DirichletBC(Q, outflow, boundary, ID)
        bc_p.append(bc)
        NS_expressions[ID] = outflow
        if MPI.rank(MPI.comm_world) == 0:
            print(("Boundary ID={:d}, pressure: {:0.6f}, area fraction: {:0.4f}".format(ID, p_initial, area_ratio[i])))

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
    common_path = path.join(folder, "Solutions")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path):
            makedirs(common_path)

    file_p = path.join(common_path, "p.h5")
    file_u = path.join(common_path, "u.h5")
    file_u_mean = path.join(common_path, "u_mean.h5")
    file_mesh = path.join(common_path, "mesh.h5")
    files = {"u": file_u, "u_mean": file_u_mean, "p": file_p, "mesh": file_mesh}

    return files


def pre_solve_hook(mesh, V, Q, newfolder, mesh_path, restart_folder, velocity_degree, cardiac_cycle,
                   save_solution_after_cycle, dt, **NS_namespace):
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

    # Save mesh as HDF5 file for post processing
    with HDF5File(MPI.comm_world, files["mesh"], "w") as mesh_file:
        mesh_file.write(mesh, "mesh")

    # Create vector function for storing velocity
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    U = Function(Vv)
    u_mean = Function(Vv)
    u_mean0 = Function(V)
    u_mean1 = Function(V)
    u_mean2 = Function(V)

    # Tstep when solutions for post processing should start being saved
    save_solution_at_tstep = int(cardiac_cycle * save_solution_after_cycle / dt)

    return dict(eval_dict=eval_dict, boundary=boundary, n=n, U=U, u_mean=u_mean, u_mean0=u_mean0, u_mean1=u_mean1, u_mean2=u_mean2, save_solution_at_tstep=save_solution_at_tstep)


def temporal_hook(u_, p_, mesh, tstep, save_probe_frequency, eval_dict, newfolder, id_in, id_out, boundary, n,
                  save_solution_frequency, NS_parameters, NS_expressions, area_ratio, t, save_solution_at_tstep,
                  U, area_inlet, nu, u_mean0, u_mean1, u_mean2, **NS_namespace):
    # Update boundary condition to current time
    for uc in NS_expressions["inlet"]:
        uc.set_t(t)

    # Compute flux and update pressure condition
    if tstep > 2:
        Q_ideals, Q_in, Q_outs = update_pressure_condition(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n,
                                                           tstep, u_)

    # Compute flow rates and updated pressure at outlets, and mean velocity and Reynolds number at inlet
    if MPI.rank(MPI.comm_world) == 0 and tstep % 10 == 0:
        U_mean = Q_in / area_inlet[0]
        diam_inlet = np.sqrt(4 * area_inlet[0] / np.pi)
        Re = U_mean * diam_inlet / nu
        print("=" * 10, "Time step " + str(tstep), "=" * 10)
        print("Sum of Q_out = {:0.4f}, Q_in = {:0.4f}, mean velocity (inlet): {:0.4f}, Reynolds number (inlet): {:0.4f}"
              .format(sum(Q_outs), Q_in, U_mean, Re))
        for i, out_id in enumerate(id_out):
            print(("For outlet with boundary ID={:d}: target flow rate: {:0.4f} mL/s, " +
                   "computed flow rate: {:0.4f} mL/s, pressure updated to: {:0.4f}")
                  .format(out_id, Q_ideals[i], Q_outs[i], NS_expressions[out_id].p))
        print()

    # Sample velocity and pressure in points/probes
    eval_dict["centerline_u_x_probes"](u_[0])
    eval_dict["centerline_u_y_probes"](u_[1])
    eval_dict["centerline_u_z_probes"](u_[2])
    eval_dict["centerline_p_probes"](p_)

    # Store sampled velocity and pressure
    if tstep % save_probe_frequency == 0:
        # Save variables along the centerline for CFD simulation
        # diagnostics and light-weight post processing
        filepath = path.join(newfolder, "Probes")
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

    # Save velocity and pressure for post processing
    if tstep % save_solution_frequency == 0 and tstep >= save_solution_at_tstep:
        # Assign velocity components to vector solution
        assign(U.sub(0), u_[0])
        assign(U.sub(1), u_[1])
        assign(U.sub(2), u_[2])

        # Get save paths
        files = NS_parameters['files']
        file_mode = "w" if tstep == save_solution_at_tstep else "a"
        p_path = files['p']
        u_path = files['u']

        # Save pressure
        viz_p = HDF5File(MPI.comm_world, p_path, file_mode=file_mode)
        viz_p.write(p_, "/pressure", tstep)
        viz_p.close()

        # Save velocity
        viz_u = HDF5File(MPI.comm_world, u_path, file_mode=file_mode)
        viz_u.write(U, "/velocity", tstep)
        viz_u.close()
    
        # Accumulate velocity
        u_mean0.vector().axpy(1, u_[0].vector())
        u_mean1.vector().axpy(1, u_[1].vector())
        u_mean2.vector().axpy(1, u_[2].vector())

def theend_hook(u_mean, u_mean0, u_mean1, u_mean2, T, dt, save_solution_at_tstep, save_solution_frequency, **NS_namespace):

    # get the file path
    files = NS_parameters['files']
    u_mean_path = files["u_mean"]

    # divide the accumlated veloicty by the number of steps
    NumTStepForAverage = (T/dt - save_solution_at_tstep) / save_solution_frequency + 1
    u_mean0.vector()[:] = u_mean0.vector()[:] /  NumTStepForAverage 
    u_mean1.vector()[:] = u_mean1.vector()[:] /  NumTStepForAverage 
    u_mean2.vector()[:] = u_mean2.vector()[:] /  NumTStepForAverage 

    assign(u_mean.sub(0), u_mean0)
    assign(u_mean.sub(1), u_mean1)
    assign(u_mean.sub(2), u_mean2)

    # Save u_mean
    with HDF5File(MPI.comm_world, u_mean_path, "w") as u_mean_file:
        u_mean_file.write(u_mean, "u_mean")


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
