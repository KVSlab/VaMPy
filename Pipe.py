import pickle
from os import path, makedirs, getcwd
import pprint
import matplotlib.pyplot as plt
import numpy as np
from mshr import *
from dolfin import *
import numpy as np

from oasis.problems.NSfracStep import *
#from ..NSfracStep import *

set_log_level(50)

def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder

    else:
        T = 5
        nu = 0.005  # 
        dt = 0.05
        NS_parameters.update(
            
            # Fluid parameters
            nu=nu,  # Kinematic viscosity
            # Simulation parameters
            T= T,
            dt=dt,  # Time step size [ms]
            U0 = 0.50,   # m/s
            save_step=1,
            checkpoint=500,  
            print_intermediate_info=100,
            folder="results_pipe_parabolic",
            mesh_path=commandline_kwargs["mesh_path"],
            # Solver parameters
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True,
            krylov_solvers=dict(monitor_convergence=False)
        )


def mesh(mesh_path, **NS_namespace):

    mesh = Mesh()
    mesh_file =  XDMFFile(MPI.comm_world, mesh_path)
    mesh_file.read(mesh)

    return mesh

    
def create_bcs(V, Q, sys_comp, nu, U0, mesh, mesh_path, **NS_namespace):
    info_red("Creating boundary conditions")
    
    # Geometry
    L = mesh.coordinates()[:, 0].max() - mesh.coordinates()[:, 0].min()
    D = mesh.coordinates()[:, 2].max() - mesh.coordinates()[:, 2].min()
    R = D/2
    
    # Mark geometry
    inlet = AutoSubDomain(lambda x, b: b and x[0] < 0.0 + DOLFIN_EPS * 1)
    wall = AutoSubDomain(lambda x, b: b and ((x[1]**2 + x[2]**2) >= R**2 - 0.01 ))
    outlet = AutoSubDomain(lambda x, b: b and (x[0] >= L - DOLFIN_EPS * 1000))
    
    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, mesh.domains())
    boundary.set_all(0)
    inlet.mark(boundary, 1)
    outlet.mark(boundary, 2)
    wall.mark(boundary, 3)

    # File("boundary.pvd") << boundary

    id_inlet = 1
    id_outlet = 2
    id_wall = 3

    # Reynolds
    print("U0 is {}".format(U0))
    Re = U0 * 2 * R / nu
    print("Reynolds number is {}".format(Re))

    noslip = Constant(0.0)
    inlet_in = Constant(U0) 
    
    # parabolic
    R2 = R**2
    inlet_in = Expression('U*( 1 - ( pow(x[0] - x0, 2) + pow(x[1] - x1, 2) + pow(x[2] -x2, 2) ) / R2 )',
                degree=2, U=U0, R2=R2, x0=0.0, x1=0.0, x2=0.0)

    bcu_wall = DirichletBC(V, Constant(0.0), boundary, id_wall)
    bcu_inlet_x = DirichletBC(V, inlet_in, boundary, id_inlet) 
    bcu_inlet_y = DirichletBC(V, Constant(0.0), boundary, id_inlet)
    bcu_inlet_z = DirichletBC(V, Constant(0.0), boundary, id_inlet)
    bcp_outlet = DirichletBC(Q, Constant(0.0), boundary, id_outlet)

    bcs = dict((ui, []) for ui in sys_comp)
    bcs['u0'] = [bcu_inlet_x, bcu_wall]
    bcs['u1'] = [bcu_inlet_y, bcu_wall]
    bcs['u2'] = [bcu_inlet_z, bcu_wall]
    bcs["p"] = [bcp_outlet]

    return bcs


def temporal_hook(mesh, u_, p_, q_, t, tstep, save_step, T, dt, **NS_namespace):
        
    if tstep % save_step == 0:
        info_green('Time = {0:2.4e}, timestep = {1:6d}, End time = {2:2.4e}'.format(t, tstep, T))
