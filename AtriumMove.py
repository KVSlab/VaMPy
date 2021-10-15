import json
import pickle

import numpy as np
from oasis.problems.NSfracStep.MovingCommon import mesh_velocity_setup, get_visualization_files, mesh_velocity_solve
from oasis.problems.NSfracStep.Womersley import compute_boundary_geometry_acrn
from scipy.interpolate import splrep, splev

from ..NSfracStep import *

set_log_level(99)


def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
    else:
        # Override some problem specific parameters
        # parameters are in mm and ms
        NS_parameters.update(
            compute_volume=False,
            mean_velocity=0.1,  # Simulate with constant mean velocity
            rigid_walls=True,  # Run rigid wall simulation
            # Fluid parameters
            # Viscosity [nu_inf: 0.0035 Pa-s / 1060 kg/m^3 = 3.3018868E-6 m^2/s == 3.3018868E-3 mm^2/ms]
            nu=3.3018868e-3,
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            # Simulation parameters
            T=1000,  # 1045,  # Run simulation for 1 cardiac cycles [ms]
            dt=50,  # 0.05,  # 10 000 steps per cycle [ms]
            no_of_cycles=5.0,  # Number of cycles
            dump_stats=100,
            store_data=2000000,
            store_data_tstep=5,  # Start storing data at 1st cycle
            save_step=100,
            save_step_problem=1,
            checkpoint=20000,
            print_intermediate_info=100,
            tstep_print=1000,
            folder="results_atrium_move",
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


def mesh(mesh_path, **NS_namespace):
    # Read mesh and print mesh information
    mesh = Mesh(mesh_path)
    print_mesh_information(mesh)

    return mesh


def pre_boundary_condition(mesh, mesh_path, id_out, id_in, T, **NS_namespace):
    # Variables needed during the simulation
    D = mesh.geometry().dim()
    boundary = MeshFunction("size_t", mesh, D - 1, mesh.domains())
    boundary.set_values(boundary.array() + 1)
    ds_new = Measure("ds", domain=mesh, subdomain_data=boundary)

    # Get IDs for inlet(s) and outlet(s)
    info_path = mesh_path.split(".")[0] + "_info.json"
    with open(info_path) as f:
        info = json.load(f)

    id_wall = 1
    id_out[:] = info['inlet_id']
    id_in[:] = info['outlet_ids']

    # Find corresponding areas
    area_total = 0
    for ID in id_in:
        area_total += assemble(Constant(1.0) * ds_new(ID))

    # Compute flow rate
    flow_rate_data = np.loadtxt(mesh_path.split(".")[0] + "_flux.txt")
    t = np.linspace(0, T, len(flow_rate_data))
    flow_rate = splrep(t, flow_rate_data, s=2, per=True)

    return dict(boundary=boundary, ds_new=ds_new, area_total=area_total, flow_rate=flow_rate, id_wall=id_wall)


class Surface_counter(UserExpression):
    def __init__(self, points, cycle, **kwargs):
        self.motion = {}
        self.counter = -1
        self.points = points
        self.time = np.linspace(0, cycle, self.points.shape[-1])
        super().__init__(**kwargs)

    def get_motion(self):
        return self.motion

    def eval(self, _, x):
        self.counter += 1
        index = np.argmin(np.sqrt(np.sum((self.points[:, :, 0] - np.array(x)) ** 2, axis=1)))
        s = 0.01
        x_ = splrep(self.time, self.points[index, 0, :], s=s, per=True)
        y_ = splrep(self.time, self.points[index, 1, :], s=s, per=True)
        z_ = splrep(self.time, self.points[index, 2, :], s=s, per=True)

        self.motion[self.counter] = [x_, y_, z_]


class InletParabolic(UserExpression):
    def __init__(self, tstep, dt, N, n, Q_profile, center, R2, area, area_total, mean_velocity, **kwargs):
        self.tstep = tstep
        self.dt = dt
        self.N = N
        self.Q_profile = Q_profile
        self.Q = 0
        self.center = center
        self.R2 = R2
        self.area = area
        self.area_total = area_total
        self.normal_component = n
        self.mean_velocity = mean_velocity
        super().__init__(**kwargs)

    def update(self, tstep):
        self.tstep = tstep
        tstep = self.tstep % self.N
        tmp_scale = 10  # TODO: Only for testing
        self.Q = splev(tstep * self.dt, self.Q_profile) / tmp_scale

    def eval(self, value, x):
        if self.mean_velocity is not None:
            U0 = self.mean_velocity
        else:
            U0 = self.Q * self.area / self.area_total

        x0 = self.center[0]
        x1 = self.center[1]
        x2 = self.center[2]
        R2 = self.R2
        parabolic = U0 * (1 - ((x0 - x[0]) ** 2 + (x1 - x[1]) ** 2 + (x2 - x[2]) ** 2) / R2)

        value[:] = - self.normal_component * parabolic


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


def create_bcs(NS_expressions, x_, boundary, V, Q, T, mesh, mesh_path, id_in, id_out, id_wall,
               tstep, dt, area_total, flow_rate, no_of_cycles, t, mean_velocity, **NS_namespace):
    # Create inlet boundary conditions
    N = int(T / dt)
    bc_inlets = {}
    for i, ID in enumerate(id_in):
        tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, id_in[i], boundary)
        inlet, coeffs = [], []
        for normal_component in tmp_normal:
            R2 = tmp_area / np.pi  # R**2
            in_ = InletParabolic(tstep, dt, N, normal_component, flow_rate, tmp_center, R2, tmp_area, area_total,
                                 mean_velocity, element=V.ufl_element())
            inlet.append(in_)

        NS_expressions[ID] = inlet
        bc_inlet = [DirichletBC(V, inlet[i], boundary, ID) for i in range(3)]
        bc_inlets[ID] = bc_inlet

    # Set outlet boundary conditions
    bc_p = []
    for i, ID in enumerate(id_out):
        bc = DirichletBC(Q, Constant(0.0), boundary, ID)
        bc_p.append(bc)
        NS_expressions['P'] = bc

    print("Loading displacement points")
    # Get mesh motion
    points = np.load(mesh_path.split(".")[0] + "_points.np", allow_pickle=True)
    cycle = T * no_of_cycles

    # Define wall movement
    surface_counter = Surface_counter(points, cycle, element=V.ufl_element())
    bc_tmp = DirichletBC(V, surface_counter, boundary, id_wall)
    bc_tmp.apply(x_["u0"])
    x_["u0"].zero()

    counter_max = surface_counter.counter
    motion = surface_counter.get_motion()

    # Remove explicitly from memory.
    del surface_counter
    del points

    # Outlet
    noslip = Constant(0.0)
    outlet_bc = DirichletBC(V, noslip, boundary, id_out[0])

    # Wall (mesh) motion
    for i, coor in enumerate(['x', 'y', 'z']):
        wall_motion = Wall_motion(t, motion, counter_max, i, cycle, element=V.ufl_element())
        NS_expressions["wall_%s" % coor] = wall_motion
        NS_expressions["outlet_%s" % coor] = outlet_bc

    #  Fluid velocity at walls
    bcu_wall_x = bcu_wall_y = bcu_wall_z = DirichletBC(V, noslip, boundary, id_wall)

    # Create lists with all boundary conditions
    bc_u0 = []
    bc_u1 = []
    bc_u2 = []
    for ID in id_in:
        bc_u0.append(bc_inlets[ID][0])
        bc_u1.append(bc_inlets[ID][1])
        bc_u2.append(bc_inlets[ID][2])
    bc_u0.append(bcu_wall_x)
    bc_u1.append(bcu_wall_y)
    bc_u2.append(bcu_wall_z)

    return dict(u0=bc_u0, u1=bc_u1, u2=bc_u2, p=bc_p)


def pre_solve_hook(mesh, V, id_in, id_out, id_wall, newfolder, u_components, velocity_degree, x_, boundary,
                   **NS_namespace):
    # Create point for evaluation
    n = FacetNormal(mesh)
    viz_p, viz_u = get_visualization_files(newfolder)

    print("Setup mesh velocity")
    L_mesh, a_mesh, coordinates, dof_map, mesh_sol, u_vec \
        = mesh_velocity_setup(V, mesh, u_components, velocity_degree, x_)

    # Set wall motion BCS
    bc_mesh = dict((ui, []) for ui in u_components)

    bcu_in = []
    noslip = Constant(0.0)
    for i, ID in enumerate(id_in):
        bcu = DirichletBC(V, noslip, boundary, ID)
        bcu_in.append(bcu)

    bcu_wall_x = DirichletBC(V, NS_expressions["wall_x"], boundary, id_wall)
    bcu_wall_y = DirichletBC(V, NS_expressions["wall_y"], boundary, id_wall)
    bcu_wall_z = DirichletBC(V, NS_expressions["wall_z"], boundary, id_wall)

    bc_out_x = DirichletBC(V, noslip, boundary, id_out[0])
    bc_out_y = DirichletBC(V, noslip, boundary, id_out[0])
    bc_out_z = DirichletBC(V, noslip, boundary, id_out[0])

    bc_mesh["u0"] = [bcu_wall_x] + [bc_out_x] + bcu_in
    bc_mesh["u1"] = [bcu_wall_y] + [bc_out_y] + bcu_in
    bc_mesh["u2"] = [bcu_wall_z] + [bc_out_z] + bcu_in

    h = EdgeLength(mesh)

    return dict(boundary=boundary, viz_u=viz_u, viz_p=viz_p, u_vec=u_vec, n=n, h=h, bc_mesh=bc_mesh, mesh_sol=mesh_sol,
                a_mesh=a_mesh, L_mesh=L_mesh, coordinates=coordinates, dof_map=dof_map)


def update_prescribed_motion(t, dt, wx_, rigid_walls, w_, u_components, bc_mesh, coordinates, dof_map, mesh_sol,
                             NS_expressions, A_cache, a_mesh, L_mesh, id_in, tstep, **NS_namespace):
    if rigid_walls:
        update_prescribed_motion.not_implemented = True
        return False

    # Update wall motion BCs
    for dir in ["x", "y", "z"]:
        NS_expressions["wall_{}".format(dir)].t = t
        NS_expressions["outlet_{}".format(dir)].t = t

    for ID in id_in:
        for i in [0, 1, 2]:
            NS_expressions[ID][i].update(tstep)

    move = mesh_velocity_solve(A_cache, L_mesh, a_mesh, bc_mesh, coordinates, dof_map, dt, mesh_sol, u_components, w_,
                               wx_)
    return move


def temporal_hook(u_, mesh, newfolder, tstep, t, save_step_problem, u_vec, viz_u, nu, dt, compute_volume,
                  **NS_namespace):
    if compute_volume:
        ComputeVolume(mesh, t, newfolder)
    if tstep % save_step_problem == 0:
        assign(u_vec.sub(0), u_[0])
        assign(u_vec.sub(1), u_[1])
        assign(u_vec.sub(2), u_[2])

        viz_u.write(u_vec, t)

    if tstep % save_step_problem * 5 == 0:
        L_mitral = 20
        h = mesh.hmin()
        ComputeReynoldsNumberAndCFL(u_, L_mitral, nu, mesh, t, tstep, dt, h, dynamic_mesh=False)


def ComputeVolume(mesh, t, newfolder):
    volume = assemble(Constant(1.0) * dx(mesh))
    data = [t, volume]
    write_data_to_file(newfolder, data, "atrium_volume.txt")


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
