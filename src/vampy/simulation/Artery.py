import json
import pickle
from os import makedirs
from pprint import pprint

import numpy as np
from dolfin import set_log_level
from oasis.problems.NSfracStep import *

from vampy.simulation.Probe import Probes  # type: ignore
from vampy.simulation.simulation_common import (
    get_file_paths,
    print_mesh_information,
    store_u_mean,
)
from vampy.simulation.Womersley import (
    compute_boundary_geometry_acrn,
    make_womersley_bcs,
)

# FEniCS specific command to control the desired level of logging, here set to critical errors
set_log_level(50)


def problem_parameters(
    commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace
):
    """
    Problem file for running CFD simulation in arterial models consisting of one inlet, and two or more outlets.
    A Womersley velocity profile is applied at the inlet, and a flow split pressure condition is applied at the outlets,
    following [1]. Flow rate for the inlet condition, and flow split values for the outlets are computed from the
    pre-processing script automatedPreProcessing.py. The simulation is run for two cycles (adjustable), but only the
    results/solutions from the second cycle are stored to avoid non-physiological effects from the first cycle.
    One cardiac cycle is set to 0.951 s from [2], and scaled by a factor of 1000, hence all parameters are in [mm] or
    [ms].

    [1] Gin, Ron, Anthony G. Straatman, and David A. Steinman. "A dual-pressure boundary condition for use in
        simulations of bifurcating conduits." J. Biomech. Eng. 124.5 (2002): 617-619.
    [2] Hoi, Yiemeng, et al. "Characterization of volumetric flow rate waveforms at the carotid bifurcations of older
        adults." Physiological measurement 31.3 (2010): 291.
    """
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, "params.dat"), "rb")
        NS_parameters.update(pickle.load(f))
        NS_parameters["restart_folder"] = restart_folder
    else:
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
            krylov_solvers=dict(monitor_convergence=False),
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


def create_bcs(
    t,
    NS_expressions,
    V,
    Q,
    area_ratio,
    area_inlet,
    mesh,
    mesh_path,
    nu,
    id_in,
    id_out,
    pressure_degree,
    **NS_namespace,
):
    # Mesh function
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())

    # Read case parameters
    info = mesh_path.split(".xml")[0] + "_info.json"
    with open(info) as f:
        info = json.load(f)

    id_in[:] = info["inlet_id"]
    id_out[:] = info["outlet_ids"]
    id_wall = min(id_in + id_out) - 1

    Q_mean = info["mean_flow_rate"]

    area_ratio[:] = info["area_ratio"]
    area_inlet.append(info["inlet_area"])

    # Load normalized time and flow rate values
    t_values, Q_ = np.loadtxt(
        path.join(path.dirname(path.abspath(__file__)), "ICA_values")
    ).T
    Q_values = (
        Q_mean * Q_
    )  # Specific flow rate = Normalized flow wave form * Prescribed flow rate
    t_values *= 1000  # Scale time in normalised flow wave form to [ms]
    tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(
        mesh, id_in[0], boundary
    )

    # Create Womersley boundary condition at inlet
    inlet = make_womersley_bcs(
        t_values,
        Q_values,
        mesh,
        nu,
        tmp_area,
        tmp_center,
        tmp_radius,
        tmp_normal,
        V.ufl_element(),
    )
    NS_expressions["inlet"] = inlet

    # Initialize inlet expressions with initial time
    for uc in inlet:
        uc.set_t(t)

    # Create pressure boundary condition
    area_out = []
    ds = Measure("ds", domain=mesh, subdomain_data=boundary)
    for i, ind in enumerate(id_out):
        dsi = ds(ind)
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
            print(
                (
                    "Boundary ID={:d}, pressure: {:0.6f}, area fraction: {:0.4f}".format(
                        ID, p_initial, area_ratio[i]
                    )
                )
            )

    # No slip condition at wall
    wall = Constant(0.0)

    # Create Boundary conditions for the velocity
    bc_wall = DirichletBC(V, wall, boundary, id_wall)
    bc_inlet = [DirichletBC(V, inlet[i], boundary, id_in[0]) for i in range(3)]

    # Return boundary conditions in dictionary
    return dict(
        u0=[bc_inlet[0], bc_wall],
        u1=[bc_inlet[1], bc_wall],
        u2=[bc_inlet[2], bc_wall],
        p=bc_p,
    )


# Oasis hook called before simulation start
def pre_solve_hook(
    mesh,
    V,
    Q,
    newfolder,
    mesh_path,
    restart_folder,
    velocity_degree,
    cardiac_cycle,
    save_solution_after_cycle,
    dt,
    **NS_namespace,
):
    # Mesh function
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())

    # Create point for evaluation
    n = FacetNormal(mesh)
    eval_dict = {}
    rel_path = mesh_path.split(".xml")[0] + "_probe_point"
    probe_points = np.load(
        rel_path, encoding="latin1", fix_imports=True, allow_pickle=True
    )

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

    # Save mesh as HDF5 file for post-processing
    with HDF5File(MPI.comm_world, files["mesh"], "w") as mesh_file:
        mesh_file.write(mesh, "mesh")

    # Create vector function for storing velocity
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    U = Function(Vv)
    u_mean = Function(Vv)
    u_mean0 = Function(V)
    u_mean1 = Function(V)
    u_mean2 = Function(V)

    # Time step when solutions for post-processing should start being saved
    save_solution_at_tstep = int(cardiac_cycle * save_solution_after_cycle / dt)

    return dict(
        eval_dict=eval_dict,
        boundary=boundary,
        n=n,
        U=U,
        u_mean=u_mean,
        u_mean0=u_mean0,
        u_mean1=u_mean1,
        u_mean2=u_mean2,
        save_solution_at_tstep=save_solution_at_tstep,
    )


# Oasis hook called after each time step
def temporal_hook(
    u_,
    p_,
    mesh,
    tstep,
    dump_probe_frequency,
    eval_dict,
    newfolder,
    id_in,
    id_out,
    boundary,
    n,
    save_solution_frequency,
    NS_parameters,
    NS_expressions,
    area_ratio,
    t,
    save_solution_at_tstep,
    U,
    area_inlet,
    nu,
    u_mean0,
    u_mean1,
    u_mean2,
    **NS_namespace,
):
    # Update boundary condition to current time
    for uc in NS_expressions["inlet"]:
        uc.set_t(t)

    # Compute flux and update pressure condition
    if tstep > 2:
        Q_ideals, Q_in, Q_outs = update_pressure_condition(
            NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n, tstep, u_
        )

    # Compute flow rates and updated pressure at outlets, and mean velocity and Reynolds number at inlet
    if MPI.rank(MPI.comm_world) == 0 and tstep % 10 == 0:
        U_mean = Q_in / area_inlet[0]
        diam_inlet = np.sqrt(4 * area_inlet[0] / np.pi)
        Re = U_mean * diam_inlet / nu
        print("=" * 10, "Time step " + str(tstep), "=" * 10)
        print(
            f"Sum of Q_out = {sum(Q_outs):0.4f}, Q_in = {Q_in:0.4f},"
            + f" mean velocity (inlet): {U_mean:0.4f},"
            + f" Reynolds number (inlet): {Re:0.4f}"
        )
        for i, out_id in enumerate(id_out):
            print(
                f"For outlet with boundary ID={out_id:d}: target flow rate: {Q_ideals[i]:0.4f} mL/s, "
                + f"computed flow rate: {Q_outs[i]:0.4f} mL/s, pressure updated to: {NS_expressions[out_id].p:0.4f}"
            )
        print()

    # Sample velocity and pressure in points/probes
    eval_dict["centerline_u_x_probes"](u_[0])
    eval_dict["centerline_u_y_probes"](u_[1])
    eval_dict["centerline_u_z_probes"](u_[2])
    eval_dict["centerline_p_probes"](p_)

    # Store sampled velocity and pressure
    if tstep % dump_probe_frequency == 0:
        # Save variables along the centerline for CFD simulation
        # diagnostics and light-weight post-processing
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

    # Save velocity and pressure for post-processing
    if tstep % save_solution_frequency == 0 and tstep >= save_solution_at_tstep:
        # Assign velocity components to vector solution
        assign(U.sub(0), u_[0])
        assign(U.sub(1), u_[1])
        assign(U.sub(2), u_[2])

        # Get save paths
        files = NS_parameters["files"]
        p_path = files["p"]
        u_path = files["u"]
        file_mode = "w" if not path.exists(p_path) else "a"

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


# Oasis hook called after the simulation has finished
def theend_hook(
    u_mean,
    u_mean0,
    u_mean1,
    u_mean2,
    T,
    dt,
    save_solution_at_tstep,
    save_solution_frequency,
    **NS_namespace,
):
    store_u_mean(
        T,
        dt,
        save_solution_at_tstep,
        save_solution_frequency,
        u_mean,
        u_mean0,
        u_mean1,
        u_mean2,
        NS_parameters,
    )


def beta(err, p):
    """
    Adjusted choice of beta for the dual-pressure boundary condition.
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
            return 1 - 5 * err**2
    else:
        if err >= 0.1:
            return 1.5
        else:
            return 1 + 5 * err**2


def update_pressure_condition(
    NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n, tstep, u_
):
    """
    Use a dual-pressure boundary condition as pressure condition at outlet.
    """
    ds = Measure("ds", domain=mesh, subdomain_data=boundary)
    Q_in = abs(assemble(dot(u_, n) * ds(id_in[0], subdomain_data=boundary)))
    Q_outs = []
    Q_ideals = []
    for i, out_id in enumerate(id_out):
        Q_out = abs(assemble(dot(u_, n) * ds(out_id)))
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

        # 1) Linear update to converge first 100 time steps of first cycle
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
                NS_expressions[out_id].p = p_old * beta(R_err, p_old) * M_err**E

    return Q_ideals, Q_in, Q_outs
