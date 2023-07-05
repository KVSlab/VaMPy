import json
import pickle
from os import makedirs
from pprint import pprint

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers
from scipy.interpolate import splrep, splev
from scipy.spatial import distance

from vampy.simulation.Probe import Probes  # type: ignore
from vampy.simulation.Womersley import make_womersley_bcs, compute_boundary_geometry_acrn
from vampy.simulation.simulation_common import store_u_mean, get_file_paths, print_mesh_information

# FEniCS specific command to control the desired level of logging, here set to critical errors
set_log_level(50)


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
    else:
        # Override some problem specific parameters
        # Parameters are in mm and ms
        cardiac_cycle = float(commandline_kwargs.get("cardiac_cycle", 1000))
        number_of_cycles = float(commandline_kwargs.get("number_of_cycles", 1))

        NS_parameters.update(
            # Moving atrium parameters
            dynamic_mesh=True,  # Run moving mesh simulation
            compute_velocity_and_pressure=True,  # Only solve mesh equations
            # Backflow parameters
            backflow_beta=0.2,
            backflow_facets=[],
            # Fluid parameters
            nu=3.3018868e-3,  # Viscosity [nu: 0.0035 Pa-s / 1060 kg/m^3 = 3.3018868E-6 m^2/s == 3.3018868E-3 mm^2/ms]
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            # Simulation parameters
            cardiac_cycle=cardiac_cycle,  # Run simulation for 1 cardiac cycles [ms]
            # FIXME: For scaling only
            T=10,  # Total simulation length
            dt=0.2,  # # Time step size [ms]
            # Frequencies to save data
            dump_probe_frequency=500,  # Dump frequency for sampling velocity & pressure at probes along the centerline
            save_solution_frequency=5,  # Save frequency for velocity and pressure field
            save_solution_frequency_xdmf=5,  # Save frequency for velocity and pressure field
            save_solution_after_cycle=0,  # Store solution after 1 cardiac cycle
            save_volume_frequency=1e10,  # Save frequency for storing volume
            save_flow_metrics_frequency=4e10,  # Frequency for storing flow metrics
            # Oasis specific parameters
            checkpoint=50000,  # Overwrite solution in Checkpoint folder each checkpoint
            print_intermediate_info=200,
            folder="results_moving_atrium",
            mesh_path=commandline_kwargs["mesh_path"],
            # Solver parameters
            max_iter=2,
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True
        )

    if MPI.rank(MPI.comm_world) == 0:
        print("=== Starting simulation for MovingAtrium.py ===")
        print("Running with the following parameters:")
        pprint(NS_parameters)


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


class Surface_counter(UserExpression):
    def __init__(self, points, cycle, **kwargs):
        self.motion = {}
        self.motion_mapping = {}
        self.counter = -1
        self.points = points
        self.num_points = self.points.shape[0]
        self.num_samples = self.points.shape[-1]
        self.time = np.linspace(0, cycle, self.num_samples)
        self.motion = self.precompute_splines()
        super().__init__(**kwargs)

    def get_motion_mapping(self):
        return self.motion_mapping

    def precompute_splines(self):
        motion = {}
        s = 2.5  # Inria Atrium
        for index in range(self.num_points):
            x_ = splrep(self.time, self.points[index, 0, :], s=s, per=True)
            y_ = splrep(self.time, self.points[index, 1, :], s=s, per=True)
            z_ = splrep(self.time, self.points[index, 2, :], s=s, per=True)
            motion[index] = [x_, y_, z_]
        return motion

    def eval(self, _, x):
        self.counter += 1
        index = distance.cdist([x], self.points[:, :, 0]).argmin()
        self.motion_mapping[self.counter] = self.motion[index]


class Wall_motion(UserExpression):
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


def create_bcs(NS_expressions, dynamic_mesh, x_, cardiac_cycle, backflow_facets, mesh, mesh_path, nu, t,
               V, Q, id_in, id_out, **NS_namespace):
    if mesh_path.endswith(".xml") or mesh_path.endswith(".xml.gz"):
        mesh_filename = ".xml"
    elif mesh_path.endswith(".h5"):
        mesh_filename = ".h5"

    rank = MPI.rank(MPI.comm_world)
    coords = ['x', 'y', 'z']
    # Variables needed during the simulation
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
    backflow_facets[:] = info['outlet_id']

    # Find corresponding areas
    ds_new = Measure("ds", domain=mesh, subdomain_data=boundary)
    area_total = 0
    for ID in id_in:
        area_total += assemble(Constant(1.0) * ds_new(ID))

    # Load normalized time and flow rate values
    if dynamic_mesh:
        flow_rate_path = mesh_path.split(mesh_filename)[0] + "_flowrate_moving.txt"
    else:
        flow_rate_path = mesh_path.split(mesh_filename)[0] + "_flowrate_rigid.txt"

    t_values, Q_ = np.loadtxt(flow_rate_path).T
    t_values *= 1000  # Scale time in normalised flow wave form to [ms]

    for ID in id_in:
        tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, ID, boundary)
        Q_scaled = Q_ * tmp_area / area_total

        # Create Womersley boundary condition at inlet
        inlet = make_womersley_bcs(t_values, Q_scaled, mesh, nu, tmp_area, tmp_center, tmp_radius, tmp_normal,
                                   V.ufl_element())
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
        wall_counter = Surface_counter(points, cardiac_cycle, element=V.ufl_element())
        bc_tmp = DirichletBC(V, wall_counter, boundary, id_wall)
        bc_tmp.apply(x_["u0"])
        x_["u0"].zero()

        wall_counter_max = wall_counter.counter
        wall_motion = wall_counter.get_motion_mapping()

        # Remove explicitly from memory.
        del wall_counter

        if rank == 0:
            print("Creating wall boundary conditions")

        for i, coord in enumerate(coords):
            wall_ = Wall_motion(t, wall_motion, wall_counter_max, i, cardiac_cycle, element=V.ufl_element())
            NS_expressions["wall_%s" % coord] = wall_

        bc_wall = [DirichletBC(V, NS_expressions["wall_%s" % coord], boundary, id_wall) for coord in coords]
    else:
        # No slip on walls
        bc_wall = [DirichletBC(V, Constant(0.0), boundary, id_wall)] * len(coords)

    # Create lists with all velocity boundary conditions
    bc_u0 = [bc_inlets[ID][0] for ID in id_in] + [bc_wall[0]]
    bc_u1 = [bc_inlets[ID][1] for ID in id_in] + [bc_wall[1]]
    bc_u2 = [bc_inlets[ID][2] for ID in id_in] + [bc_wall[2]]

    return dict(u0=bc_u0, u1=bc_u1, u2=bc_u2, p=bc_p)


def pre_solve_hook(u_components, id_in, id_out, dynamic_mesh, V, Q, cardiac_cycle, dt,
                   save_solution_after_cycle, mesh_path, mesh, newfolder, velocity_degree,
                   restart_folder, **NS_namespace):
    id_wall = min(id_in + id_out) - 1
    # Extract diameter at mitral valve
    mesh_filename=".h5"
    info_path = mesh_path.split(mesh_filename)[0] + "_info.json"
    with open(info_path) as f:
        info = json.load(f)

    outlet_area = info['outlet_area']
    D_mitral = np.sqrt(4 * outlet_area / np.pi)

    # Create point for evaluation
    n = FacetNormal(mesh)
    eval_dict = {}
    rel_path = mesh_path.split(mesh_filename)[0] + "_probe_point"
    probe_points = np.load(rel_path, encoding='latin1', fix_imports=True, allow_pickle=True)

    # Define xdmf writers
    viz_U = get_visualization_writers(newfolder, ['velocity'])[0]

    # Extract dof map and coordinates
    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

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

    # Create Probes path
    probes_folder = path.join(newfolder, "Probes")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(probes_folder):
            makedirs(probes_folder)

    # Create vector function for storing velocity
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    U = Function(Vv)
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

    return dict(outlet_area=outlet_area, id_wall=id_wall, D_mitral=D_mitral, n=n, eval_dict=eval_dict, U=U,
                u_mean=u_mean, u_mean0=u_mean0, u_mean1=u_mean1, u_mean2=u_mean2, boundary=boundary, bc_mesh=bc_mesh,
                coordinates=coordinates, viz_U=viz_U, save_solution_at_tstep=save_solution_at_tstep,
                dof_map=dof_map, probes_folder=probes_folder)


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
                  boundary, outlet_area, probes_folder, **NS_namespace):
    if tstep % save_volume_frequency == 0:
        compute_volume(mesh, t, newfolder)

    if tstep % save_solution_frequency_xdmf == 0:
        for i in range(3):
            assign(U.sub(i), u_[i])

        viz_U.write(U, t)

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
        # Save variables along the centerline for CFD simulation
        # diagnostics and light-weight post-processing
        for component in ["u_x", "u_y", "u_z", "p"]:
            arr = eval_dict[f"centerline_{component}_probes"].array()
            if MPI.rank(MPI.comm_world) == 0:
                # Dump stats
                arr.dump(path.join(probes_folder, f"{component}_{str(tstep)}.probes"))

        # Clear stats
        MPI.barrier(MPI.comm_world)
        for component in ["u_x", "u_y", "u_z", "p"]:
            eval_dict[f"centerline_{component}_probes"].clear()

    # Save velocity and pressure for post-processing
    if tstep % save_solution_frequency == 0 and tstep >= save_solution_at_tstep:
        # Assign velocity components to vector solution
        for i in range(3):
            assign(U.sub(i), u_[i])

        # Get save paths
        files = NS_parameters['files']
        p_path = files['p']
        u_path = files['u']
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
def theend_hook(u_mean, u_mean0, u_mean1, u_mean2, T, dt, save_solution_at_tstep, save_solution_frequency,
                **NS_namespace):
    store_u_mean(T, dt, save_solution_at_tstep, save_solution_frequency, u_mean, u_mean0, u_mean1, u_mean2,
                 NS_parameters)
