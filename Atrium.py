import json
import os
import pickle
from os import path, makedirs
import csv

import numpy as np
from numpy.core.fromnumeric import mean
from simulation.Womersley import make_womersley_bcs, compute_boundary_geometry_acrn

from oasis.problems.NSfracStep import *
from scipy.interpolate import UnivariateSpline
from scipy.integrate import simps, romberg
from IPython import embed
from ufl.measure import measure_names

set_log_level(50)

def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    backflow = bool(commandline_kwargs.get("backflow", True))
    backflow_beta = float(commandline_kwargs.get("backflow_beta", 0.5))
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
    else:
        # Override some problem specific parameters
        # Parameters are in mm and ms
        cardiac_cycle = 1000
        number_of_cycles = 1
        NS_parameters.update(
            # Fluid parameters
            nu=3.3018868e-3,  # Viscosity [nu_inf: 0.0035 Pa-s / 1060 kg/m^3 = 3.3018868E-6 m^2/s == 3.3018868E-3 mm^2/ms]
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            area_ratio=[],
            # Simulation parameters
            cardiac_cycle=cardiac_cycle,
            T=cardiac_cycle * number_of_cycles,  # Simulation end time [ms]# Run simulation for 1 cardiac cycles [ms]
            dt= 0.1,  # 10 000 steps per cycle [ms]
            dump_stats=100,
            store_data=5e6,
            store_data_tstep=50,  # Start storing data at 1st cycle
            save_step=50e8,
            save_step_problem=10,
            save_solution_frequency = 5,
            save_solution_after_cycle=1,  # Store solution after 1 cardiac cycle
            checkpoint=20000,
            print_intermediate_info=100,
            tstep_print = 10e5,
            folder="results_atrium",
            mesh_path=commandline_kwargs["mesh_path"],
            # Solver parameters
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True,
            krylov_solvers=dict(monitor_convergence=False),
            dim_MV = [],
            dim_PV = []
        )
        if backflow:
            NS_parameters.update(
                # Backflow
                backflow_facets=[],
                backflow_beta=backflow_beta
            )

    mesh_file = NS_parameters["mesh_path"].split("/")[-1]
    case_name = mesh_file.split(".")[0]
    NS_parameters["folder"] = path.join(NS_parameters["folder"], case_name)


def mesh(mesh_path, **NS_namespace):
    
    mesh = Mesh(mesh_path)
    #mesh_ = Mesh(mesh_path)
    #mesh = refine(mesh_)

    return mesh

class boundary_expression(UserExpression):

    def __init__(self, coeffs, omega, period, normal, normal_component, area, centers, radius, profile_type, **kwargs):
        
        self.Cn = coeffs
        self.omega = omega
        self.period = period
        self.normal = normal
        self.normal_component = normal_component
        self.area = area
        self.centers = centers
        self.radius = radius
        self.profile_type = profile_type
        self.num_fourier_coefficients = len(coeffs)
        self.t = None
        self.ns = ns = np.arange(1,self.num_fourier_coefficients)

        super().__init__(**kwargs)

    
    def update(self, t):
        self.t = float(t) % self.period
        self._expnt = np.exp((self.omega * self.t * 1j) * self.ns)

    def parabolic(self, x):
        r2 = self.radius**2
        x0=self.centers[0]
        x1=self.centers[1]
        x2=self.centers[2]
        return 1 - ( pow(x0 - x[0], 2) + pow(x1 - x[1], 2) + pow(x2 - x[2], 2) ) / r2

    def eval(self,values,x):

        # Profile
        if self.profile_type == 'uniform':
            par = (self.Cn[0] + np.dot(self.Cn[1:], self._expnt)).real
        elif self.profile_type == 'parabolic':
            par = (self.Cn[0]*self.parabolic(x) + np.dot(self.Cn[1:], self._expnt)).real
        
        # Scale by negative normal direction 
        values[0] = -self.normal_component * par
       

class InletParabolic(UserExpression):
    def __init__(self, n, center, area, mean_velocity, Q_profile, area_total, **kwargs):

        self.center = center
        self.R2 = area / np.pi
        self.area = area
        self.normal_component = n
        self.mean_velocity = mean_velocity
        self.area_total = area_total
        self.Q = Q_profile

        super().__init__(**kwargs)

    def eval(self, value, x):
        
        if self.mean_velocity is not None:
            U0 = self.mean_velocity
        else:
            Q0 = self.Q * self.area / self.area_total
            U0 = 2*(Q0 / self.area)

        x0 = self.center[0]
        x1 = self.center[1]
        x2 = self.center[2]
        R2 = self.R2
        parabolic = U0 * (1 - ((x0 - x[0]) ** 2 + (x1 - x[1]) ** 2 + (x2 - x[2]) ** 2) / R2)

        value[:] = - self.normal_component * parabolic
        
def create_bcs(t, NS_expressions, V, Q, area_ratio, mesh, mesh_path, nu, backflow, backflow_facets, id_in, id_out,velocity_degree, pressure_degree, **NS_namespace):
    # Mesh function
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
    boundary.set_values(boundary.array() + 1)
    ds_new = Measure("ds", domain=mesh, subdomain_data=boundary)

    # print(boundary.array())

    # Get IDs for inlet(s) and outlet(s)
    info_path = mesh_path.split(".")[0] + "_info.json"

    with open(info_path) as f:
        info = json.load(f)

    id_wall = 1
    id_out[:] = info['outlet_id']
    id_in[:] = info['inlet_ids']
    if backflow:
        backflow_facets[:] = id_out
        
    area_ratio = info['area_ratio']

    # Find corresponding areas
    area_total = 0
    for ID in id_in:
        area_total += assemble(Constant(1.0) * ds_new(ID))
        
    # No slip condition at wall
    wall = Constant(0.0)
    # Create Boundary conditions for the wall
    bc_wall = DirichletBC(V, wall, boundary, id_wall)

    # Get inlet and outlet areas from json file
    area_in = [info[key] for key in info.keys() if (key.startswith('inlet') and key.endswith('area'))]
    area_out = [info[key] for key in info.keys() if (key.startswith('outlet') and key.endswith('area'))]
    # Transform into np.array
    area_in = np.array(area_in)
    area_out = np.array(area_out)

    dim_MV = np.sqrt(4*area_in.max()/np.pi)  #[mm]
    dim_PV = np.sqrt(4*area_out/np.pi)
    NS_parameters['dim_MV'] = dim_MV
    NS_parameters['dim_PV'] = dim_PV
 
    #Will implement Analytical profile
    t_values , V_values = [], [] 
    """try:
        t_values, V_values = np.loadtxt(path.join(path.dirname(path.abspath(__file__)), "pv_velocity.txt")).T
        t_values *= 1000
        V_values *=1
    except ValueError:
        raise"""
    
    flow_rate = 5.5 * 50/3  #mm3/ms  # 1 [l/min] = 50 / 3 [mm3/ms]
    mean_velocity = 0.3 #m/s
    bc_inlets = {}
    for i, ID in enumerate(id_in):
        tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, id_in[i], boundary)  
        inlet, coeffs = [], []
        #coeffs, omega, period = get_coeffients(t_values, V_values)
        for normal_component in tmp_normal:
            # _in = boundary_expression(coeffs, omega, period, tmp_normal, normal_component, tmp_area, tmp_center, tmp_radius,'parabolic',element = V.ufl_element())
            _in = InletParabolic( normal_component, tmp_center, tmp_area, None, flow_rate, area_total, element=V.ufl_element()) #parabolic profile
            inlet.append(_in)
        
        NS_expressions[ID] = inlet
        bc_inlet = [DirichletBC(V, inlet[i], boundary, ID) for i in range(3)]
        bc_inlets[ID] = bc_inlet
        # Set start time equal to t_0
        # for uc in inlet:
        #     uc.update(t)

    bc_p = []   
    for i, ID in enumerate(id_out):
        bc = DirichletBC(Q, Constant(0.0), boundary, ID)
        bc_p.append(bc)
        NS_expressions['P'] = bc
   

    # Create lists with all boundary conditions
    bc_u0 = []
    bc_u1 = []
    bc_u2 = []
    for ID in id_in:
        bc_u0.append(bc_inlets[ID][0])
        bc_u1.append(bc_inlets[ID][1])
        bc_u2.append(bc_inlets[ID][2])
    bc_u0.append(bc_wall)
    bc_u1.append(bc_wall)
    bc_u2.append(bc_wall)
    
    return dict(u0=bc_u0, u1=bc_u1, u2=bc_u2, p=bc_p)

def get_file_paths(folder):
    # Create folder where data and solutions (velocity, mesh, pressure) is stored
    common_path = path.join(folder, "PostProc")
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
    
    # Mesh function
    boundary = MeshFunction("size_t", mesh, 2, mesh.domains())
    boundary.set_values(boundary.array() + 1)
    n = FacetNormal(mesh)

    # Create point for evaluation
    """eval_dict = {}
    rel_path = mesh_path.split(".")[0] + "_probe_point"
    probe_points = np.load(rel_path, encoding='latin1', fix_imports=True, allow_pickle=True)

    # Store points file in checkpoint
    if MPI.rank(MPI.comm_world) == 0:
       probe_points.dump(path.join(newfolder, "Checkpoint", "points"))

    eval_dict["centerline_u_x_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_y_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_u_z_probes"] = Probes(probe_points.flatten(), V)
    eval_dict["centerline_p_probes"] = Probes(probe_points.flatten(), Q)"""

    if restart_folder is None:
        # Get files to store results
        files = get_file_paths(newfolder)
        NS_parameters.update(dict(files=files))
    else:
        files = NS_namespace["files"]

    # Save mesh as HDF5 file for post processing
    with HDF5File(MPI.comm_world, files["mesh"], "w") as mesh_file:
        mesh_file.write(mesh, "mesh")
    
    h = EdgeLength(mesh)
    
    # Create vector function for storing velocity
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    U = Function(Vv)
    u_mean = Function(Vv)
    u_mean0 = Function(V)
    u_mean1 = Function(V)
    u_mean2 = Function(V)
    u_vec = Function(Vv, name="u")
    viz_p, viz_u = get_visualization_files(newfolder)

    # Tstep when solutions for post processing should start being saved
    save_solution_at_tstep = int(cardiac_cycle * save_solution_after_cycle / dt)

    return dict(boundary=boundary, n=n, h=h, viz_uu=viz_u, viz_pp=viz_p, u_vec=u_vec, U=U, 
                u_mean=u_mean, u_mean0=u_mean0, u_mean1=u_mean1, u_mean2=u_mean2, save_solution_at_tstep=save_solution_at_tstep)


def u_dot_n(u, n):
    return (dot(u, n) - abs(dot(u, n))) / 2


def velocity_tentative_hook(mesh, boundary, u_ab, x_1, b, A, ui, u, v, backflow_facets, backflow_beta,
                            **NS_namespace):
    boundary = MeshFunction("size_t", mesh, 2, mesh.domains())
    boundary.set_values(boundary.array() + 1)

    if backflow_facets != []:
        ds = Measure("ds", domain=mesh, subdomain_data=boundary)
        n = FacetNormal(mesh)
        K = assemble(u_dot_n(u_ab, n) * dot(u, v) * ds(backflow_facets[0]))
        A.axpy(-backflow_beta * 0.5, K, True)
        b[ui].axpy(backflow_beta * 0.5, K * x_1[ui])


def temporal_hook(h, u_, q_, p_, mesh, tstep, save_step_problem,dump_stats, newfolder, id_in, id_out, boundary, n, store_data,
                  NS_parameters, NS_expressions, area_ratio, dt, t, tstep_print, save_solution_frequency, save_solution_at_tstep, u_vec,
                   U, viz_uu, viz_pp, u_mean0, u_mean1, u_mean2, **NS_namespace):
   
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
    boundary.set_values(boundary.array() + 1)

   
    # Update boundary condition
    # for i, in_id in enumerate(id_in):
    #     for uc in NS_expressions[in_id]:
    #         uc.update(t)
                
    # Compute flux and update pressure condition
    if tstep >=1 and tstep %tstep_print == 0:
        #Q_in, Q_ins, Q_out, V_out, V_ins = compute_flow_rates(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n, tstep, u_, newfolder, t)
        Q_ins, Q_outs, V_ins, V_outs = compute_flow_rates(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n, tstep, u_, newfolder, t)

    if MPI.rank(MPI.comm_world) == 0 and tstep >=1 and tstep %tstep_print == 0:
        velocity_path = path.join(newfolder, "Solutions", "velocity.txt")
        flowrate_path = path.join(newfolder, "Solutions", "flowrate.txt")
        # if not path.isdir(path.join(newfolder, "Solutions")):
        #     os.mkdir(path.join(newfolder, "Solutions"))
        s = "{:2.4e} "
        for i in range(len(V_outs) + len(V_ins)):
            s = s + "{:.4f} "
        s = s + "\n"

        with open(velocity_path, 'a') as filename:
            filename.write(s.format(*[t, *V_outs, *V_ins]))
            #filename.write("{:2.4e}, {:.4f}, {:.4f}, {:.4f}, {:.4f}, {:.4f} \n".format(t, V_out, V_ins[0], V_ins[1], V_ins[2], V_ins[3]))
        with open(flowrate_path, 'a') as filename:
            #fname.write("{:2.4e}, {:.4f}, {:.4f}, {:.4f}, {:.4f}, {:.4f}, {:.4f} \n".format(t, Q_out, Q_ins[0], Q_ins[1], Q_ins[2], Q_ins[3], sum(Q_ins)))  
            filename.write(s.format(*[t, *Q_outs, *Q_ins]))

    if tstep % tstep_print == 0:
        DG = FunctionSpace(mesh, "DG", 0)
        U_ = project(sqrt(inner(u_, u_)), DG, solver_type='cg')

        cfl = U_.vector().get_local() * dt / h   
        
        max_cfl = cfl.max()
        min_cfl = cfl.min() 

        dim_MV = NS_parameters['dim_MV']  
        dim_PV = NS_parameters['dim_PV'] 
        
        Re =  U_.vector().get_local() * dim_MV / NS_parameters['nu']
        Re_MV = Re.max()
        Re_ =  U_.vector().get_local() * dim_PV / NS_parameters['nu']
        Re_PV = Re_.max()

    if MPI.rank(MPI.comm_world) == 0 and tstep % tstep_print == 0:  #print_intermediate_info       
        info_green('Time = {:0.4f}, timestep = {:0.4f}, max_CFL={:0.4f}, min_CFL={:0.4f}, Re_PV={:0.4f}, Re_MV={:0.4f}'
                    .format(t, tstep, max_cfl, min_cfl, Re_PV, Re_MV))
        print("Sum of Q_in = {:0.4f} Q_out = {:0.4f}".format(sum(Q_ins), Q_out))

    if tstep % save_step_problem == 0:
        assign(u_vec.sub(0), u_[0])
        assign(u_vec.sub(1), u_[1])
        assign(u_vec.sub(2), u_[2])

        viz_uu.write(u_vec, t)
        viz_pp.write(p_, t)
        
    # # Sample velocity in points
    # eval_dict["centerline_u_x_probes"](u_[0])
    # eval_dict["centerline_u_y_probes"](u_[1])
    # eval_dict["centerline_u_z_probes"](u_[2])
    # eval_dict["centerline_p_probes"](p_)

    # # Store sampled velocity
    # if tstep % dump_stats == 0:
    #     filepath = path.join(newfolder, "Stats")
    #     if MPI.rank(MPI.comm_world) == 0:
    #         if not path.exists(filepath):
    #             makedirs(filepath)

    #     arr_u_x = eval_dict["centerline_u_x_probes"].array()
    #     arr_u_y = eval_dict["centerline_u_y_probes"].array()
    #     arr_u_z = eval_dict["centerline_u_z_probes"].array()
    #     arr_p = eval_dict["centerline_p_probes"].array()

    #     #Dump stats
    #     if MPI.rank(MPI.comm_world) == 0:
    #        arr_u_x.dump(path.join(filepath, "u_x_%s.probes" % str(tstep)))
    #        arr_u_y.dump(path.join(filepath, "u_y_%s.probes" % str(tstep)))
    #        arr_u_z.dump(path.join(filepath, "u_z_%s.probes" % str(tstep)))
    #        arr_p.dump(path.join(filepath, "p_%s.probes" % str(tstep)))

    #     #Clear stats
    #     MPI.barrier(MPI.comm_world)
    #     eval_dict["centerline_u_x_probes"].clear()
    #     eval_dict["centerline_u_y_probes"].clear()
    #     eval_dict["centerline_u_z_probes"].clear()

    # Save velocity and pressure
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

def get_visualization_files(newfolder):
    viz_u = XDMFFile(MPI.comm_world, path.join(newfolder, "Solutions", "u.xdmf"))
    viz_p = XDMFFile(MPI.comm_world, path.join(newfolder, "Solutions", "p.xdmf"))
    for viz in [viz_u, viz_p]:
        viz.parameters["rewrite_function_mesh"] = True
        viz.parameters["flush_output"] = True
        viz.parameters["functions_share_mesh"] = True
    return viz_p, viz_u

def compute_flow_rates(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n, tstep, u_, newfolder,t):
    
    V = FunctionSpace(mesh, "DG", 1)
    f = Function(V)
    f.vector()[:] = 1.

    Q_outs = []
    V_outs = []
    timer = Timer("compute_flow_rates - Compute V_out")
    for i, out_id in enumerate(id_out):
        Q_out = abs(assemble(dot(u_, n) * ds(out_id, domain=mesh, subdomain_data=boundary)))
        dso = assemble(f*ds(out_id, domain=mesh, subdomain_data=boundary))
        # print('id_out:', id_out[0], dso)
        V_out = Q_out/ dso
        Q_outs.append(Q_out)
        V_outs.append(V_out)
    timer.stop()
    print("Time compute_flow_rates - Compute V_out = ", timer.elapsed()[0])
    
    f.vector()[:] = 1.
    Q_ins = []
    Q_ideals = []
    V_ins = []
    timer = Timer("compute_flow_rates - Compute V_in")
    for i, in_id in enumerate(id_in):
        Q_in = abs(assemble(dot(u_, n) * ds(in_id, domain=mesh, subdomain_data=boundary)))
        dsi = assemble(f*ds(in_id, domain=mesh, subdomain_data=boundary))
        # print(in_id, dsi)
        V_in = Q_in / dsi
        Q_ins.append(Q_in)
        V_ins.append(V_in)
        #Q_ideal = area_ratio[i] * Q_in  #FIXIT: area_ratio
        #Q_ideals.append(Q_ideal)
    timer.stop()
    print("Time compute_flow_rates - Compute V_in = ", timer.elapsed()[0])

    return Q_ins, Q_outs, V_ins, V_outs
    #return Q_in, Q_ins, Q_out, V_out, V_ins 
    #return Q_ideals, Q_in, Q_ins, Q_out, V_out, V_ins

def EdgeLength(mesh):
    # Compute edge length
    DG = FunctionSpace(mesh, "DG", 0)
    circumradius = Circumradius(mesh)
    circumradius_DG = project(circumradius, DG).vector()[:]

    if mesh.geometry().dim() == 2:
        edge_length = circumradius_DG * 2
    elif mesh.geometry().dim() == 3:
        edge_length = circumradius_DG * sqrt(8 / 3)

    return edge_length


def get_coeffients(x,y):

    num_fourier_coefficients = 20
    
    N = len(x) + 1
    # Compute transient profile as interpolation of given coefficients
    period = max(x)
    transient_profile = UnivariateSpline(x, y, s=0, k=1)

    # Compute fourier coefficients of transient profile
    timedisc = np.linspace(0, period, N)

    Cn, omega = fourier_coefficients(timedisc, transient_profile, period, num_fourier_coefficients)

    return Cn, omega, period

def fourier_coefficients(x, y, T, N):

    '''From x-array and y-spline and period T, calculate N complex Fourier coefficients.'''
    omega = 2*np.pi/T
    ck = []
    ck.append(1/T*simps(y(x), x))
    for n in range(1,N):
        c = 1/T*simps(y(x)*np.exp(-1j*n*omega*x), x)

        # Clamp almost zero real and imag components to zero
        if 1:
            cr = c.real
            ci = c.imag
            if abs(cr) < 1e-14: cr = 0.0
            if abs(ci) < 1e-14: ci = 0.0
            c = cr + ci*1j

        ck.append(2*c)

    return ck, omega
