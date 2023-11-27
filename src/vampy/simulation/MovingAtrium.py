import json
import pickle
from dolfin import set_log_level, UserExpression
from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers
from os import makedirs
from scipy.interpolate import splrep, splev
from scipy.spatial import KDTree
from vampy.simulation.Probe import Probes  # type: ignore
from vampy.simulation.Womersley import make_womersley_bcs, compute_boundary_geometry_acrn
from vampy.simulation.simulation_common import get_file_paths, print_mesh_information, \
    store_velocity_and_pressure_h5, dump_probes

# FEniCS specific command to control the desired level of logging, here set to critical errors
set_log_level(50)


def problem_parameters(commandline_kwargs, scalar_components, NS_parameters, **NS_namespace):
    """
    Problem file for running CFD simulation in left atrial models consisting of arbitrary number of pulmonary veins (PV)
    (normally 3 to 7), and one outlet being the mitral valve. A Womersley velocity profile is applied at the inlets,
    where the total flow rate is split between the area ratio of the PVs. The mitral valve is considered open with a
    constant pressure of p=0. Flow rate and flow split values for the inlet condition are computed from the
    pre-processing script automatedPreProcessing.py. The simulation is run for two cycles (adjustable), but only the
    results/solutions from the second cycle are stored to avoid non-physiological effects from the first cycle.
    One cardiac cycle is set to 0.951 s from [1], and scaled by a factor of 1000, hence all parameters are in
    [mm] or [ms].
    The script has been adjusted to handle moving domains.

    [1] Hoi, Yiemeng, et al. "Characterization of volumetric flow rate waveforms at the carotid bifurcations of older
        adults." Physiological measurement 31.3 (2010): 291.
    """

    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        mesh_path = commandline_kwargs["mesh_path"]
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
        NS_parameters['mesh_path'] = mesh_path
        track_blood = NS_parameters['track_blood']
        globals().update(NS_parameters)
    else:
        # Override some problem specific parameters
        # Parameters are in mm and ms
        cardiac_cycle = float(commandline_kwargs.get("cardiac_cycle", 1000))
        number_of_cycles = int(commandline_kwargs.get("number_of_cycles", 2))
        track_blood = bool(commandline_kwargs.get("track_blood", True))

        NS_parameters.update(
            # Moving atrium parameters
            use_supg=True,  # SUPG for transport equation
            dynamic_mesh=True,  # Run moving mesh simulation
            compute_velocity_and_pressure=True,  # Only solve mesh equations (For volume computation)
            # Blood residence time
            track_blood=track_blood,
            # Backflow parameters
            backflow_beta=0.4,
            backflow_facets=[],
            # Fluid parameters
            nu=3.3018868e-3,  # Viscosity [nu: 0.0035 Pa-s / 1060 kg/m^3 = 3.3018868E-6 m^2/s == 3.3018868E-3 mm^2/ms]
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            # Simulation parameters
            cardiac_cycle=cardiac_cycle,  # Run simulation for 1 cardiac cycles [ms]
            T=cardiac_cycle * number_of_cycles,  # Total simulation length
            dt=0.1,  # # Time step size [ms]
            # Frequencies to save data
            dump_probe_frequency=500,  # Dump frequency for sampling velocity & pressure at probes along the centerline
            save_solution_frequency=20,  # Save frequency for velocity and pressure field
            save_solution_frequency_xdmf=1e10,  # Save frequency for velocity and pressure field
            save_solution_after_cycle=0,  # Store solution after 1 cardiac cycle
            save_volume_frequency=1e10,  # Save frequency for storing volume
            save_flow_metrics_frequency=50,  # Frequency for storing flow metrics
            # Oasis specific parameters
            checkpoint=1000,  # Overwrite solution in Checkpoint folder each checkpoint
            print_intermediate_info=500,
            folder="results_moving_atrium",
            mesh_path=commandline_kwargs["mesh_path"],
            # Solver parameters
            max_iter=2,
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True,
            krylov_solvers={'monitor_convergence': False,
                            'report': False,
                            'relative_tolerance': 1e-8,
                            'absolute_tolerance': 1e-8}
        )

    if track_blood:
        if MPI.rank(MPI.comm_world) == 0:
            print("-- Computing blood residence time --")
        scalar_components += ["blood"]


def scalar_source(scalar_components, **NS_namespace):
    """Return a dictionary of scalar sources."""
    return dict((ci, Constant(1)) for ci in scalar_components)


def mesh(mesh_path, **NS_namespace):
    # Read mesh and print mesh information
    if mesh_path.endswith(".xml") or mesh_path.endswith(".xml.gz"):
        mesh = Mesh(mesh_path)
    elif mesh_path.endswith(".h5"):
        mesh = Mesh()
        with HDF5File(MPI.comm_world, mesh_path, "r") as f:
            f.read(mesh, "mesh", False)

    print_mesh_information(mesh)

    return mesh


class MeshMotionMapping(UserExpression):
    def __init__(self, points, cycle, **kwargs):
        self.N = 100
        self.initial_points = points
        self.points = self.interpolate_points()
        self.motion_mapping = {}
        self.counter = -1
        self.num_points = self.points.shape[0]
        self.num_samples = self.points.shape[-1]
        self.time = np.linspace(0, cycle, self.num_samples)
        self.tree = self.create_kd_tree()
        super().__init__(**kwargs)

    def get_motion_mapping(self):
        return self.motion_mapping

    def create_kd_tree(self):
        tree = KDTree(self.points[:, :, 0])
        return tree

    def interpolate_points(self):
        points = self.initial_points
        time = np.linspace(0, 1, points.shape[2])
        time_r = np.linspace(0, 1, self.N)
        move = np.zeros((points.shape[0], points.shape[1], number_of_points + 1))

        # Use interp1d if smooth displacement
        x = interp1d(time, points[:, 0, :], axis=1)
        y = interp1d(time, points[:, 1, :], axis=1)
        z = interp1d(time, points[:, 2, :], axis=1)

        move[:, 0, :] = x(time_r)
        move[:, 1, :] = y(time_r)
        move[:, 2, :] = z(time_r)

        return move

    def eval(self, _, x):
        self.counter += 1
        _, index = self.tree.query(x)
        # FIXME: Set spline parameter objectively
        s = 0.5
        x_ = splrep(self.time, self.points[index, 0, :], s=s, per=True)
        y_ = splrep(self.time, self.points[index, 1, :], s=s, per=True)
        z_ = splrep(self.time, self.points[index, 2, :], s=s, per=True)

        self.motion_mapping[self.counter] = [x_, y_, z_]


class WallMotion(UserExpression):
    def __init__(self, t, motion, max_counter, direction, cycle, **kwargs):
        self.t = t
        self.motion = motion
        self.max_counter = max_counter
        self.counter = -1
        self.direction = direction
        self.cycle = cycle
        super().__init__(**kwargs)

    def eval(self, values, _):
        self.counter += 1
        # No motion for inlet/outlet
        values[:] = splev(self.t % self.cycle, self.motion[self.counter][self.direction], der=1)
        if self.counter == self.max_counter:
            self.counter = -1


def pre_boundary_condition(boundary, **NS_namespace):
    return dict(boundary=boundary)


def create_bcs(NS_expressions, dynamic_mesh, x_, cardiac_cycle, backflow_facets, mesh, mesh_path, nu, t,
               V, Q, id_in, id_out, track_blood, boundary, **NS_namespace):
    if mesh_path.endswith(".xml") or mesh_path.endswith(".xml.gz"):
        mesh_filename = ".xml"
    elif mesh_path.endswith(".h5"):
        mesh_filename = ".h5"

    rank = MPI.rank(MPI.comm_world)
    coords = ['x', 'y', 'z']
    # Variables needed during the simulation
    if boundary is None:
        boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
        if mesh_path.endswith(".h5"):
            with HDF5File(MPI.comm_world, mesh_path, "r") as f:
                f.read(boundary, "boundary")

    # Get IDs for inlet(s) and outlet(s)
    info_path = mesh_path.split(mesh_filename)[0] + "_info.json"
    with open(info_path) as f:
        info = json.load(f)

    id_in[:] = info['inlet_ids']
    id_out[:] = info['outlet_id']
    id_wall = min(id_in + id_out) - 1

    # Apply backflow at the outlet (Mitral valve)
    backflow_facets[:] = info['outlet_id']

    # Find corresponding areas
    ds_new = Measure("ds", domain=mesh, subdomain_data=boundary)
    area_total = 0
    for ID in id_in:
        area_total += assemble(Constant(1.0) * ds_new(ID))

    setup = "moving" if dynamic_mesh else "rigid"
    # Look for case specific flow rate
    flow_rate_path = mesh_path.split(mesh_filename)[0] + f"_flowrate_{setup}.txt"
    if path.exists(flow_rate_path):
        t_values, Q_ = np.loadtxt(flow_rate_path).T
    else:
        # Use default flow rate
        flow_rate_path = path.join(path.dirname(path.abspath(__file__)), "PV_values")
        # Load normalized time and flow rate values
        Q_mean = info['mean_flow_rate']
        t_values, Q_ = np.loadtxt(flow_rate_path).T
        Q_ *= Q_mean

    t_values *= 1000  # Scale time in normalised flow wave form to [ms]

    for ID in id_in:
        tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, ID, boundary)
        Q_scaled = Q_ * tmp_area / area_total

        # Create Womersley boundary condition at inlets
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

    # Set wall boundary conditions
    if dynamic_mesh:
        # Moving walls
        if rank == 0:
            print("Loading displacement points")
        points = np.load(mesh_path.split(mesh_filename)[0] + "_points.np", allow_pickle=True)
        if rank == 0:
            print("Creating splines for displacement")

        # Define wall movement
        motion_mapping = MeshMotionMapping(points, cardiac_cycle, element=V.ufl_element())
        bc_tmp = DirichletBC(V, motion_mapping, boundary, id_wall)
        bc_tmp.apply(x_["u0"])
        x_["u0"].zero()

        wall_counter_max = motion_mapping.counter
        wall_motion = motion_mapping.get_motion_mapping()

        # Remove explicitly from memory.
        del motion_mapping
        del points

        if rank == 0:
            print("Creating wall boundary conditions")

        for i, coord in enumerate(coords):
            wall_ = WallMotion(t, wall_motion, wall_counter_max, i, cardiac_cycle, element=V.ufl_element())
            NS_expressions["wall_%s" % coord] = wall_

        bc_wall = [DirichletBC(V, NS_expressions["wall_%s" % coord], boundary, id_wall) for coord in coords]
    else:
        # No slip on walls
        bc_wall = [DirichletBC(V, Constant(0.0), boundary, id_wall)] * len(coords)

    # Create lists with all velocity boundary conditions
    bc_u0 = [bc_inlets[ID][0] for ID in id_in] + [bc_wall[0]]
    bc_u1 = [bc_inlets[ID][1] for ID in id_in] + [bc_wall[1]]
    bc_u2 = [bc_inlets[ID][2] for ID in id_in] + [bc_wall[2]]

    if track_blood:
        bc_blood = [DirichletBC(V, Constant(0.0), boundary, ID) for ID in id_in]

        return dict(u0=bc_u0, u1=bc_u1, u2=bc_u2, p=bc_p, blood=bc_blood)
    else:
        return dict(u0=bc_u0, u1=bc_u1, u2=bc_u2, p=bc_p)


def pre_solve_hook(u_components, id_in, id_out, dynamic_mesh, V, Q, cardiac_cycle, dt,
                   save_solution_after_cycle, mesh_path, mesh, newfolder, velocity_degree,
                   restart_folder, **NS_namespace):
    id_wall = min(id_in + id_out) - 1
    # Extract diameter at mitral valve
    if mesh_path.endswith(".xml") or mesh_path.endswith(".xml.gz"):
        mesh_filename = ".xml"
    elif mesh_path.endswith(".h5"):
        mesh_filename = ".h5"

    info_path = mesh_path.split(mesh_filename)[0] + "_info.json"
    with open(info_path) as f:
        info = json.load(f)

    outlet_area = info['outlet_area']
    D_mitral = np.sqrt(4 * outlet_area / np.pi)

    # Create point for evaluation
    n = FacetNormal(mesh)
    eval_dict = {}
    rel_path = mesh_path.split(mesh_filename)[0] + "_probe_point.json"
    with open(rel_path, 'r') as infile:
        probe_points = np.array(json.load(infile))

    # Store points file in checkpoint
    if MPI.rank(MPI.comm_world) == 0:
        probe_points.dump(path.join(newfolder, "Checkpoint", "points"))

    eval_dict["centerline_u_x_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_y_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_z_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_p_probes"] = Probes(probe_points.flatten(), Q)

    if restart_folder is None:
        # Get files to store results
        files = get_file_paths(newfolder, ['brt', 'd'])
        NS_parameters.update(dict(files=files))
    else:
        files = NS_namespace["files"]

    # Save mesh as HDF5 file for post-processing
    if not path.exists(files['half']['mesh']):
        with HDF5File(MPI.comm_world, files['half']["mesh"], "w") as mesh_file:
            mesh_file.write(mesh, "mesh")

    # Create Probes path
    probes_folder = path.join(newfolder, "Probes")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(probes_folder):
            makedirs(probes_folder)

    # Define xdmf writers
    viz_U, viz_blood = get_visualization_writers(newfolder, ['velocity', 'blood'])

    # Extract dof map and coordinates
    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    # Create vector function for storing velocity and deformation
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    U = Function(Vv)
    D = Function(Vv)
    u_mean = Function(Vv)
    u_mean0 = Function(V)
    u_mean1 = Function(V)
    u_mean2 = Function(V)

    # Time step when solutions for post-processing should start being saved
    save_solution_at_tstep = int(cardiac_cycle * save_solution_after_cycle / dt)

    # Set mesh equation boundary condition
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
    if dynamic_mesh:
        # Set wall motion BCS
        bc_mesh = dict((ui, []) for ui in u_components)

        # Zero wall velocity at inlet and outlet
        noslip = Constant(0.0)
        bc_out = [DirichletBC(V, noslip, boundary, ID) for ID in id_out]
        bcu_in = [DirichletBC(V, noslip, boundary, ID) for ID in id_in]

        # Add wall movement to wall
        bcu_wall_x = DirichletBC(V, NS_expressions["wall_x"], boundary, id_wall)
        bcu_wall_y = DirichletBC(V, NS_expressions["wall_y"], boundary, id_wall)
        bcu_wall_z = DirichletBC(V, NS_expressions["wall_z"], boundary, id_wall)

        bc_mesh["u0"] = [bcu_wall_x] + bc_out + bcu_in
        bc_mesh["u1"] = [bcu_wall_y] + bc_out + bcu_in
        bc_mesh["u2"] = [bcu_wall_z] + bc_out + bcu_in
    else:
        bc_mesh = {}

    return dict(outlet_area=outlet_area, id_wall=id_wall, D_mitral=D_mitral, n=n, eval_dict=eval_dict, U=U, D=D,
                u_mean=u_mean, u_mean0=u_mean0, u_mean1=u_mean1, u_mean2=u_mean2, boundary=boundary, bc_mesh=bc_mesh,
                coordinates=coordinates, viz_U=viz_U, save_solution_at_tstep=save_solution_at_tstep, dof_map=dof_map,
                probes_folder=probes_folder, viz_blood=viz_blood)


def update_boundary_conditions(t, dynamic_mesh, NS_expressions, id_in, **NS_namespace):
    # Update inlet condition
    for ID in id_in:
        for i in [0, 1, 2]:
            NS_expressions[f"inlet_{ID}"][i].set_t(t)

    if dynamic_mesh:
        # Update wall motion BCs
        for coord in ["x", "y", "z"]:
            NS_expressions[f"wall_{coord}"].t = t


def temporal_hook(mesh, id_wall, id_out, cardiac_cycle, dt, t, save_solution_frequency, u_, id_in, tstep, newfolder,
                  eval_dict, dump_probe_frequency, p_, save_solution_at_tstep, nu, D_mitral, U, u_mean0, u_mean1,
                  u_mean2, save_flow_metrics_frequency, save_volume_frequency, save_solution_frequency_xdmf, viz_U,
                  boundary, outlet_area, q_, w_, viz_blood, track_blood, D, du_, T, **NS_namespace):
    if tstep % save_volume_frequency == 0:
        compute_volume(mesh, t, newfolder)

    if tstep % save_solution_frequency_xdmf == 0:
        for i in range(3):
            assign(U.sub(i), u_[i])

        viz_U.write(U, t)
        if track_blood:
            viz_blood.write(q_['blood'], t)

    if tstep % save_flow_metrics_frequency == 0:
        h = mesh.hmin()
        compute_flow_quantities(u_, D_mitral, nu, mesh, t, tstep, dt, h, outlet_area, boundary, id_out, id_in, id_wall,
                                period=cardiac_cycle, newfolder=newfolder, dynamic_mesh=False, write_to_file=True)

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
        files = NS_parameters['files']
        if t <= (T / 2):
            files = files['half']
        else:
            files = files['full']

        store_velocity_and_pressure_h5(files, U, p_, tstep, u_, u_mean0, u_mean1, u_mean2, D, du_, q_['blood'])
