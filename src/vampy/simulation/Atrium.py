import importlib.util
import json
import pickle
from pprint import pprint

import numpy as np
from dolfin import set_log_level, MPI

from vampy.simulation.Probe import Probes  # type: ignore
from vampy.simulation.Womersley import make_womersley_bcs, compute_boundary_geometry_acrn
from vampy.simulation.simulation_common import store_u_mean, get_file_paths, print_mesh_information, \
    store_velocity_and_pressure_h5, dump_probes

# Check for oasis and oasismove
package_name_oasis = 'oasis'
package_name_oasismove = 'oasismove'
oasis_exists = importlib.util.find_spec(package_name_oasis)
oasismove_exists = importlib.util.find_spec(package_name_oasismove)
if oasismove_exists:
    from oasismove.problems.NSfracStep import *
elif oasis_exists:
    from oasis.problems.NSfracStep import *
else:
    print("Neither oasis nor oasismove is installed. Exiting simulation..")

# FEniCS specific command to control the desired level of logging, here set to critical errors
set_log_level(50)
comm = MPI.comm_world


def problem_parameters(commandline_kwargs, NS_parameters, scalar_components, Schmidt, NS_expressions, **NS_namespace):
    """
    Problem file for running CFD simulation in left atrial models consisting of arbitrary number of pulmonary veins (PV)
    (normally 3 to 7), and one outlet being the mitral valve. A Womersley velocity profile is applied at the inlets,
    where the total flow rate is split between the area ratio of the PVs. The mitral valve is considered open with a
    constant pressure of p=0. Flow rate and flow split values for the inlet condition are computed from the
    pre-processing script automatedPreProcessing.py. The simulation is run for two cycles (adjustable), but only the
    results/solutions from the second cycle are stored to avoid non-physiological effects from the first cycle.
    One cardiac cycle is set to 0.951 s from [1], and scaled by a factor of 1000, hence all parameters are in
    [mm] or [ms].

    [1] Hoi, Yiemeng, et al. "Characterization of volumetric flow rate waveforms at the carotid bifurcations of older
        adults." Physiological measurement 31.3 (2010): 291.
    """
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
    else:
        # Override some problem specific parameters
        # Parameters are in mm and ms
        cardiac_cycle = float(commandline_kwargs.get("cardiac_cycle", 951))
        number_of_cycles = float(commandline_kwargs.get("number_of_cycles", 2))

        NS_parameters.update(
            # Fluid parameters
            nu=3.3018868e-3,  # Viscosity [nu: 0.0035 Pa-s / 1060 kg/m^3 = 3.3018868E-6 m^2/s == 3.3018868E-3 mm^2/ms]
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            # Simulation parameters
            cardiac_cycle=cardiac_cycle,  # Run simulation for 1 cardiac cycles [ms]
            T=number_of_cycles * cardiac_cycle,  # Number of cycles
            dt=0.0951,  # # Time step size [ms]
            # Frequencies to save data
            dump_probe_frequency=100,  # Dump frequency for sampling velocity & pressure at probes along the centerline
            save_solution_frequency=5,  # Save frequency for velocity and pressure field
            save_solution_after_cycle=0,  # Store solution after 1 cardiac cycle
            # Oasis specific parameters
            checkpoint=100,  # Overwrite solution in Checkpoint folder each checkpoint
            print_intermediate_info=100,
            folder="results_atrium",
            mesh_path=commandline_kwargs["mesh_path"],
            # Solver parameters
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True,
            krylov_solvers=dict(monitor_convergence=False)
        )

    mesh_file = NS_parameters["mesh_path"].split("/")[-1]
    case_name = mesh_file.split(".")[0]
    NS_parameters["folder"] = path.join(NS_parameters["folder"], case_name)

    if MPI.rank(comm) == 0:
        print("=== Starting simulation for Atrium.py ===")
        print("Running with the following parameters:")
        pprint(NS_parameters)


def mesh(mesh_path, **NS_namespace):
    # Read mesh and print mesh information
    atrium_mesh = Mesh(mesh_path)
    print_mesh_information(atrium_mesh)

    return atrium_mesh


def create_bcs(NS_expressions, mesh, mesh_path, nu, t, V, Q, id_in, id_out, **NS_namespace):
    # Variables needed during the simulation
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())

    # Get IDs for inlet(s) and outlet(s)
    info_path = mesh_path.split(".xml")[0] + "_info.json"
    with open(info_path) as f:
        info = json.load(f)

    id_in[:] = info['inlet_ids']
    id_out[:] = info['outlet_id']
    id_wall = min(id_in + id_out) - 1

    Q_mean = info['mean_flow_rate']

    # Find corresponding areas
    ds = Measure("ds", domain=mesh, subdomain_data=boundary)
    area_total = 0
    for ID in id_in:
        area_total += assemble(Constant(1.0) * ds(ID))

    # Load normalized time and flow rate values
    t_values, Q_ = np.loadtxt(path.join(path.dirname(path.abspath(__file__)), "PV_values")).T
    Q_values = Q_mean * Q_  # Specific flow rate = Normalized flow wave form * Prescribed flow rate
    t_values *= 1000  # Scale time in normalised flow wave form to [ms]

    for ID in id_in:
        tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, ID, boundary)
        Q_scaled = Q_values * tmp_area / area_total

        # Create Womersley boundary condition at inlet
        inlet = make_womersley_bcs(t_values, Q_scaled, nu, tmp_center, tmp_radius, tmp_normal, V.ufl_element())
        NS_expressions[f"inlet_{ID}"] = inlet

    # Initial condition
    for ID in id_in:
        for i in [0, 1, 2]:
            NS_expressions[f"inlet_{ID}"][i].set_t(t)

    # Create inlet boundary conditions
    bc_inlets = {}
    for ID in id_in:
        bc_inlet = [DirichletBC(V, NS_expressions[f"inlet_{ID}"][i], boundary, ID) for i in range(3)]
        bc_inlets[ID] = bc_inlet

    # Set outlet boundary conditions, assuming one outlet (Mitral Valve)
    bc_p = [DirichletBC(Q, Constant(0), boundary, ID) for ID in id_out]

    # No slip on walls
    bc_wall = [DirichletBC(V, Constant(0.0), boundary, id_wall)]

    # Create lists with all velocity boundary conditions
    bc_u0 = [bc_inlets[ID][0] for ID in id_in] + bc_wall
    bc_u1 = [bc_inlets[ID][1] for ID in id_in] + bc_wall
    bc_u2 = [bc_inlets[ID][2] for ID in id_in] + bc_wall

    return dict(u0=bc_u0, u1=bc_u1, u2=bc_u2, p=bc_p)


def pre_solve_hook(V, Q, cardiac_cycle, dt, save_solution_after_cycle, mesh_path, mesh, newfolder, velocity_degree,
                   restart_folder, **NS_namespace):
    # Extract diameter at mitral valve
    info_path = mesh_path.split(".xml")[0] + "_info.json"
    with open(info_path) as f:
        info = json.load(f)

    outlet_area = info['outlet_area']
    D_mitral = np.sqrt(4 * outlet_area / np.pi)

    # Create point for evaluation
    n = FacetNormal(mesh)
    eval_dict = {}
    rel_path = mesh_path.split(".xml")[0] + "_probe_point.json"
    with open(rel_path, 'r') as infile:
        probe_points = np.array(json.load(infile))

    # Store points file in checkpoint
    if MPI.rank(comm) == 0:
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
    with HDF5File(comm, files["mesh"], "w") as mesh_file:
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

    return dict(D_mitral=D_mitral, n=n, eval_dict=eval_dict, U=U, u_mean=u_mean, u_mean0=u_mean0, u_mean1=u_mean1,
                u_mean2=u_mean2, save_solution_at_tstep=save_solution_at_tstep)


def temporal_hook(mesh, dt, t, save_solution_frequency, u_, NS_expressions, id_in, tstep, newfolder,
                  eval_dict, dump_probe_frequency, p_, save_solution_at_tstep, nu, D_mitral, U, u_mean0, u_mean1,
                  u_mean2, **NS_namespace):
    # Update inlet condition
    for ID in id_in:
        for i in [0, 1, 2]:
            NS_expressions["inlet_{}".format(ID)][i].set_t(t)

    # Compute flow rates and updated pressure at outlets, and mean velocity and Reynolds number at inlet
    if tstep % 10 == 0:
        # Compute and printCFL number
        DG = FunctionSpace(mesh, "DG", 0)
        U_vector = project(sqrt(inner(u_, u_)), DG).vector().get_local()
        h = mesh.hmin()

        # Collect velocities in parallel
        local_U_mean = U_vector.mean()
        local_U_max = U_vector.max()

        U_mean = comm.gather(local_U_mean, 0)
        U_max = comm.gather(local_U_max, 0)

        if MPI.rank(comm) == 0:
            u_mean = np.mean(U_mean)
            u_max = max(U_max)

            Re_mean = u_mean * D_mitral / nu
            Re_max = u_max * D_mitral / nu

            CFL_mean = u_mean * dt / h
            CFL_max = u_max * dt / h
            info = 'Time = {0:2.4e}, timestep = {1:6d}, max Reynolds number={2:2.3f}, mean Reynolds number={3:2.3f}' + \
                   ', max velocity={4:2.3f}, mean velocity={5:2.3f}, max CFL={6:2.3f}, mean CFL={7:2.3f}'
            info_green(info.format(t, tstep, Re_max, Re_mean, u_max, u_mean, CFL_max, CFL_mean))

    # Sample velocity and pressure in points/probes
    eval_dict["centerline_u_x_probes"](u_[0])
    eval_dict["centerline_u_y_probes"](u_[1])
    eval_dict["centerline_u_z_probes"](u_[2])
    eval_dict["centerline_p_probes"](p_)

    # Store sampled velocity and pressure
    if tstep % dump_probe_frequency == 0:
        dump_probes(eval_dict, newfolder, tstep)

    # Save velocity and pressure for post-processing
    if tstep % save_solution_frequency == 0 and tstep >= save_solution_at_tstep:
        store_velocity_and_pressure_h5(NS_parameters, U, p_, tstep, u_, u_mean0, u_mean1, u_mean2)


# Oasis hook called after the simulation has finished
def theend_hook(u_mean, u_mean0, u_mean1, u_mean2, T, dt, save_solution_at_tstep, save_solution_frequency,
                **NS_namespace):
    store_u_mean(T, dt, save_solution_at_tstep, save_solution_frequency, u_mean, u_mean0, u_mean1, u_mean2,
                 NS_parameters)
