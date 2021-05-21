import json
import os
import pickle
from os import path, makedirs
import csv

import numpy as np
from numpy.core.fromnumeric import mean
from Womersley import make_womersley_bcs, compute_boundary_geometry_acrn
#from Parabolic import make_parabolic_bcs, compute_boundary_geometry_acrn

from fenicstools import Probes
from oasis.problems.NSfracStep import *

set_log_level(50)


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
            # Fluid parameters
            nu=3.5e-3,  # Viscosity
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            area_ratio=[],
            # Simulation parameters
            T=1045,  # Run simulation for 1 cardiac cycles [ms]
            dt=0.1, # 10 000 steps per cycle [ms]
            no_of_cycles=1,
            dump_stats=100,
            store_data=5,
            store_data_tstep=1000,  # Start storing data at 2nd cycle
            save_step=200,
            checkpoint=5000,
            print_intermediate_info=100,
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

    mesh_file = NS_parameters["mesh_path"].split("/")[-1]
    case_name = mesh_file.split(".")[0]
    NS_parameters["folder"] = path.join(NS_parameters["folder"], case_name)


def mesh(mesh_path, **NS_namespace):
    # Read mesh
    return Mesh(mesh_path)
       

def create_bcs(t, NS_expressions, V, Q, area_ratio, mesh, mesh_path, nu, id_in, id_out, pressure_degree, **NS_namespace):
    # Mesh function
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
    boundary.set_values(boundary.array() + 1)

    id_in = [2,3,5,6] # Hardcoded. FIXIT: automated prepocessing
    id_out = [4]

    # Read case parameters
    info_path = mesh_path.split(".")[0] + ".json"
    with open(info_path) as f:
        info = json.load(f)

    # Extract flow split ratios and inlet/outlet IDs
    id_info = info['idFileLine'].split()
    #id_out.append(int(id_info[1]))
    #id_in[:] = [int(p) for p in id_info[2].split(",")]
    area_ratio[:] = [float(p) for p in info['areaRatioLine'].split()[-1].split(",")]

    area_in = np.array([info['outlet0_area'],info['outlet1_area'],info['outlet2_area'],info['outlet3_area']])
    area_out = np.array([info['inlet_area']])
    dim_MV = np.sqrt(4*area_in.max()/np.pi)  #[mm]
    dim_PV = np.sqrt(4*area_out/np.pi)
    NS_parameters['dim_MV'] = dim_MV
    NS_parameters['dim_PV'] = dim_PV
 
    Q_mean = 20     # [ml]
    t_values, Q_ = np.load(path.join(path.dirname(path.abspath(__file__)), "ICA_values"))
    Q_values = Q_mean * Q_
    t_values *= 1000

    # t_values , V_values = [], [] 
    # try:
    #     with open(path.join(path.dirname(path.abspath(__file__)), "Dataset.csv")) as csvfile:    
    #         csvReader = csv.reader(csvfile)    
    #         for row in csvReader:   
    #             t_values.append(float(row[0])*1000) #convert to ms 
    #             V_values.append(float(row[1])) 
    # except ValueError:
    #     raise

    bc_inlets = {}
    for i, ID in enumerate(id_in):
        tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, id_in[i], boundary)
        inlet = make_womersley_bcs(t_values, Q_values, mesh, nu, tmp_area, tmp_center, tmp_radius, tmp_normal, V.ufl_element())
        #inlet = make_parabolic_bcs(t_values, V_values, mesh, nu, tmp_area, tmp_center, tmp_radius, tmp_normal, V.ufl_element(),coeffstype="V")
        NS_expressions[ID] = inlet
        bc_inlet = [DirichletBC(V, inlet[i], boundary, ID) for i in range(3)]
        bc_inlets[ID] = bc_inlet
        # Set start time equal to t_0
        for uc in inlet:
            uc.set_t(t)

    bc_p = []
    for i, ID in enumerate(id_out):
        bc = DirichletBC(Q, Constant(0.0), boundary, ID)
        bc_p.append(bc)
        NS_expressions['P'] = bc
   
    # No slip condition at wall
    wall = Constant(0.0)
    # Create Boundary conditions for the wall
    bc_wall = DirichletBC(V, wall, boundary, 1)

    # Need to be updated if there are more than four inlets
    return dict(u0=[bc_inlets[id_in[0]][0], bc_inlets[id_in[1]][0], bc_inlets[id_in[2]][0], bc_inlets[id_in[3]][0], bc_wall],
                u1=[bc_inlets[id_in[0]][1], bc_inlets[id_in[1]][1], bc_inlets[id_in[2]][1], bc_inlets[id_in[3]][1], bc_wall],
                u2=[bc_inlets[id_in[0]][2], bc_inlets[id_in[1]][2], bc_inlets[id_in[2]][2], bc_inlets[id_in[3]][2], bc_wall],
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
    
    h = EdgeLength(mesh)

    return dict(eval_dict=eval_dict, boundary=boundary, n=n, h=h)


def temporal_hook(h, u_, q_, p_, mesh, tstep, dump_stats, eval_dict, newfolder, id_in, id_out, boundary, n, store_data,
                  NS_parameters, NS_expressions, area_ratio, dt, t, store_data_tstep, **NS_namespace):
    
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
    boundary.set_values(boundary.array() + 1)

    id_in = [2,3,5,6] # Hardcoded. FIXIT: automated prepocessing
    id_out = [4]

    # Update boundary condition
    for i, in_id in enumerate(id_in):
        for uc in NS_expressions[in_id]:
            uc.set_t(t)
                
    # Compute flux and update pressure condition
    if tstep > 2:
        Q_ideals, Q_in, Q_ins, Q_out, V_out, V_ins = compute_flow_rates(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n,
                                                           tstep, u_, newfolder, t)

    DG = FunctionSpace(mesh, "DG", 0)
    U = project(sqrt(inner(u_, u_)), DG)

    cfl = U.vector().get_local() * dt / h   #Check the used unit in calculating cfl before running simulation
    
    max_cfl = cfl.max()
    min_cfl = cfl.min() 

    dim_MV = NS_parameters['dim_MV']  
    dim_PV = NS_parameters['dim_PV'] 
    
    # Re =  U.vector().get_local() * dim_MV / NS_parameters['nu']
    # Re_MV = Re.max()*1000
    # Re_ =  U.vector().get_local() * dim_PV / NS_parameters['nu']
    # Re_PV = Re_.max()*1000

    if MPI.rank(MPI.comm_world) == 0 and tstep % 10 == 0:
        print("=" * 10, tstep, "=" * 10)
        print()
        # info_green('Time = {0:2.4e}, timestep = {1:6d}, max_CFL={2:2.3f}, min_CFL={3:2.3f}, Re_PV={0:2.4e}, Re_MV={0:2.4e}'
        #        .format(t, tstep, max_cfl, min_cfl, Re_PV, Re_MV))
        info_green('Time = {0:2.4e}, timestep = {1:6d}, max_CFL={2:2.3f}, min_CFL={3:2.3f}'
               .format(t, tstep, max_cfl, min_cfl))
        print("Sum of Q_in = {:0.4f} Q_out = {:0.4f}".format(sum(Q_ins), Q_out))
        #for i, in_id in enumerate(id_in):
        #    print(("({:d}) New pressure {:0.4f}").format(out_id, NS_expressions[out_id].p))
        # for i, in_id in enumerate(id_out):
        #     print(("({:d}) area ratio {:0.4f}, ideal: {:0.4f} actual:" +
        #            " {:0.4f}").format(in_id, area_ratio[i], Q_ideals[i], Q_ins[i]))      
        print()
    

    if MPI.rank(MPI.comm_world) == 0 and tstep > 2:
        velocity_path = path.join(newfolder, "Timeseries", "velocity.txt")
        if not path.isdir(path.join(newfolder, "Timeseries")):
            os.mkdir(path.join(newfolder, "Timeseries"))
        with open(velocity_path, 'a') as filename:
            filename.write("{:2.4e}, {:.4f}, {:.4f}, {:.4f}, {:.4f}, {:.4f} \n".format(t, V_out, V_ins[0], V_ins[1], V_ins[2], V_ins[3]))  

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

        #Dump stats
        if MPI.rank(MPI.comm_world) == 0:
           arr_u_x.dump(path.join(filepath, "u_x_%s.probes" % str(tstep)))
           arr_u_y.dump(path.join(filepath, "u_y_%s.probes" % str(tstep)))
           arr_u_z.dump(path.join(filepath, "u_z_%s.probes" % str(tstep)))
           arr_p.dump(path.join(filepath, "p_%s.probes" % str(tstep)))

        #Clear stats
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


def compute_flow_rates(NS_expressions, area_ratio, boundary, id_in, id_out, mesh, n, tstep, u_, newfolder,t):
    
    V = FunctionSpace(mesh, "DG", 1)
    f = Function(V)
    f.vector()[:] = 1.

    Q_out = abs(assemble(dot(u_, n) * ds(id_out[0], domain=mesh, subdomain_data=boundary)))
    dso = assemble(f*ds(id_out[0], domain=mesh, subdomain_data=boundary))
    V_out = Q_out/ dso
    
    f.vector()[:] = 1.
    Q_ins = []
    Q_ideals = []
    V_ins = []
    for i, in_id in enumerate(id_in):
        Q_in = abs(assemble(dot(u_, n) * ds(in_id, domain=mesh, subdomain_data=boundary)))
        dsi = assemble(f*ds(in_id, domain=mesh, subdomain_data=boundary))
        V_in = Q_in / dsi
        Q_ins.append(Q_in)
        V_ins.append(V_in)
        Q_ideal = area_ratio[i] * Q_in  #FIXIT: area_ratio
        Q_ideals.append(Q_ideal)

    return Q_ideals, Q_in, Q_ins, Q_out, V_out, V_ins

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

