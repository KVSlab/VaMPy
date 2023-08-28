from __future__ import print_function

from pathlib import Path

from dolfin import *
from dolfin import Function, VectorFunctionSpace, FunctionSpace, parameters, MPI, HDF5File, Mesh, XDMFFile, BoundaryMesh, project, inner

from postprocessing_common_phaseaveraged import STRESS, read_command_line, get_dataset_names
from os import getcwd, makedirs, path

import json
import os
import pickle
import csv

import numpy as np
from numpy.core.fromnumeric import mean

from scipy.interpolate import UnivariateSpline
from scipy.integrate import simps, romberg
from IPython import embed
from ufl.measure import measure_names

import pickle
import random
from pprint import pprint

import fenics as fe
from numpy import tanh
from ufl.algorithms import extract_coefficients, extract_arguments
import mpi4py as pyMPI
import dolfin as df

try:
    parameters["reorder_dofs_serial"] = False
except NameError:
    pass

def compute_hemodynamic_indices(folder, nu, rheology_model, rho, dt, velocity_degree, T, save_frequency, start_cycle, step, average_over_cycles):
    """
    Loads velocity fields from completed CFD simulation,
    and computes and saves the following hemodynamic quantities:
    (1) WSS - Wall shear stress
    (2) TAWSS - Time averaged wall shear stress
    (3) TWSSG - Temporal wall shear stress gradient
    (4) OSI - Oscillatory shear index
    (5) RRT - Relative residence time
    (6) ECAP - Endothelial cell activation potential

    The resulting wall shear stress will be in units Pascal [Pa], given that the provided
    density (rho) is in [kg/m^3], the time step (dt) is in [ms], and viscosity (nu) is in [mm^2/ms].

    Args:
        velocity_degree (int): Finite element degree of velocity
        case_path (Path): Path to results from simulation
        rheology_model (str): "Newtonian" is the default case, while "non-Newtonian" is set for non-Newtonian simulations
        nu (float): Kinematic viscosity
        rho (float): Fluid density
        dt (float): Time step of simulation
    """
    # File paths
    file_path_u = path.join(folder, "u.h5")
    mesh_path = path.join(folder, "mesh.h5")
    file_u = HDF5File(MPI.comm_world, file_path_u, "r")

    # Determine what time step to start post-processing from
    start = int(T / dt / save_frequency * (start_cycle - 1))
    print("start=", start)

    # Get names of data to extract
    if MPI.rank(MPI.comm_world) == 0:
        print("Reading dataset names")

    dataset = get_dataset_names(file_u, start=start, step=step)

    
    if (rheology_model != "Newtonian"):
        file_path_nunn = path.join(folder, "nu.h5")
    #-----------------------------------------------------------------------------------------------
    # Read mesh saved as HDF5 format
    mesh = Mesh()
    with HDF5File(MPI.comm_world, mesh_path.__str__(), "r") as mesh_file:
        mesh_file.read(mesh, "mesh", False)

    # Load mesh
    bm = BoundaryMesh(mesh, 'exterior')

    if MPI.rank(MPI.comm_world) == 0:
        print("Define function spaces")
    V_b1 = VectorFunctionSpace(bm, "CG", 1)
    U_b1 = FunctionSpace(bm, "CG", 1)
    if (rheology_model != "Newtonian"):
        FunctionSpace_nunn = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)

    if MPI.rank(MPI.comm_world) == 0:
        print("Define functions")
    u = Function(V)
    if (rheology_model != "Newtonian"):
        nunn = Function(FunctionSpace_nunn)

    # RRT
    RRT = Function(U_b1)
    RRT_avg = Function(U_b1)
    
    # OSI
    OSI = Function(U_b1)
    OSI_avg = Function(U_b1)

    # ECAP
    ECAP = Function(U_b1)
    ECAP_avg = Function(U_b1)

    # WSS_mean
    WSS_mean = Function(V_b1)
    WSS_mean_avg = Function(V_b1)
    ###wss_mean = Function(U_b1)
    
    # TAWSS
    TAWSS = Function(U_b1)
    TAWSS_avg = Function(U_b1)
    ##tawss = Function(U_b1)

    # TWSSG
    TWSSG = Function(U_b1)
    TWSSG_avg = Function(U_b1)
    ###twssg_ = Function(U_b1)
    twssg = Function(V_b1)
    tau_prev = Function(V_b1)

    if (rheology_model == "Newtonian"):
        stress = STRESS(u, 0.0,   nu, rheology_model, mesh)
        print("rheology_model=", rheology_model)
    else:
        stress = STRESS(u, 0.0, nunn, rheology_model, mesh)
        print("rheology_model=", rheology_model)

    dabla = get_dabla_function()

    # Get number of saved steps and cycles
    saved_time_steps_per_cycle = int(T / dt / save_frequency / step)
    print("saved_time_steps_per_cycle =", saved_time_steps_per_cycle)
    print("len(dataset)=",len(dataset))
    n_cycles = int(len(dataset) / saved_time_steps_per_cycle)
    print("n_cycles=", n_cycles)

    # Set number of cycles to average over
    cycles_to_average = list(range(1, n_cycles + 1)) if average_over_cycles else []
    counters_to_save = [saved_time_steps_per_cycle * cycle for cycle in cycles_to_average]
    cycle_names = [""] + ["_cycle_{:02d}".format(cycle) for cycle in cycles_to_average]
    print("cycles_to_average=", cycles_to_average)
    print("counters_to_save=", counters_to_save)
    print("cycle_names=", cycle_names)
    

    # Create XDMF files for saving indices
    fullname = file_path_u.replace("u.h5", "%s%s.xdmf")
    fullname = fullname.replace("Solutions", "Hemodynamics")
    index_names = ["RRT", "OSI", "ECAP", "TAWSS", "TWSSG"]
    index_variables = [RRT, OSI, ECAP, TAWSS, TWSSG]
    index_variables_avg = [RRT_avg, OSI_avg, ECAP_avg, TAWSS_avg, TWSSG_avg]

    index_dict = dict(zip(index_names, index_variables))
    index_dict_cycle = dict(zip(index_names, index_variables_avg))

    indices = {}
    for cycle_name in cycle_names:
        for index in index_names + ["WSS"]:
            indices[index + cycle_name] = XDMFFile(MPI.comm_world, fullname % (index, cycle_name))
            indices[index + cycle_name].parameters["rewrite_function_mesh"] = False
            indices[index + cycle_name].parameters["flush_output"] = True
            indices[index + cycle_name].parameters["functions_share_mesh"] = True

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Start post processing", "=" * 10)

    counter = start
    for data in dataset:
        # Update file_counter
        #counter += 1
        counter += step

        file_u.read(u, data)

        tstep = file_u.attributes(data)["timestamp"]
        if MPI.rank(MPI.comm_world) == 0:
            print("=" * 10, "Timestep: {}".format(tstep), "=" * 10)

        # Compute WSS
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS (mean)")
        tau = stress()
        tau.vector()[:] = tau.vector()[:] * rho
        WSS_mean_avg.vector().axpy(1, tau.vector())

        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS (absolute value)")
        tawss = project(inner(tau, tau) ** (1 / 2), U_b1)
        TAWSS_avg.vector().axpy(1, tawss.vector())
        ##dabla(tau.vector(), TAWSS.vector())
        ##TAWSS_avg.vector().axpy(1, TAWSS.vector())

        # Compute TWSSG
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute TWSSG")
        twssg.vector().set_local((tau.vector().get_local() - tau_prev.vector().get_local()) / dt)
        twssg.vector().apply("insert")
        twssg_ = project(inner(twssg, twssg) ** (1 / 2), U_b1)
        TWSSG_avg.vector().axpy(1, twssg_.vector())

        # Update tau
        if MPI.rank(MPI.comm_world) == 0:
            print("Update WSS \n")
        tau_prev.vector().zero()
        tau_prev.vector().axpy(1, tau.vector())

        # Save instantaneous WSS
        tau.rename("WSS", "WSS")
        indices["WSS"].write(tau, dt * tstep)

        if len(cycles_to_average) != 0 and counter == counters_to_save[0]:
            # Get cycle number
            cycle = int(counters_to_save[0] / saved_time_steps_per_cycle)
            print("cycle=", cycle)
            if MPI.rank(MPI.comm_world) == 0:
                print("=" * 10, "Storing cardiac cycle {}".format(cycle), "=" * 10)

            # Get average over sampled time steps
            for index in [TWSSG_avg, TAWSS_avg, WSS_mean_avg]:
                index.vector()[:] = index.vector()[:] / saved_time_steps_per_cycle

            # Compute TAWSS
            wss_mean = project(inner(WSS_mean_avg, WSS_mean_avg) ** (1 / 2), U_b1)
            ###wss_mean = project(inner(WSS_mean_avg, WSS_mean_avg) ** (1 / 2), V_b1)
            #dabla(WSS_mean_avg.vector(), wss_mean.vector())
            wss_mean_vec = wss_mean.vector().get_local()
            tawss_vec = TAWSS_avg.vector().get_local()

            # Compute RRT, OSI, and ECAP based on mean and absolute WSS
            RRT_avg.vector().set_local(1 / wss_mean_vec)
            OSI_avg.vector().set_local(0.5 * (1 - wss_mean_vec / tawss_vec))
            ECAP_avg.vector().set_local(OSI_avg.vector().get_local() / tawss_vec)

            for index in [RRT_avg, OSI_avg, ECAP_avg]:
                index.vector().apply("insert")

            # Rename displayed variable names
            for var, name in zip(index_variables_avg, index_names):
                var.rename(name, name)

            # Store solution
            for name, index in index_dict_cycle.items():
                indices[name + "_cycle_{:02d}".format(cycle)].write(index)

            # Append solution to total solution
            for index, index_avg in zip(index_dict.values(), index_dict_cycle.values()):
                index_avg.vector().apply("insert")
                index.vector().axpy(1, index_avg.vector())

            WSS_mean_avg.vector().apply("insert")
            WSS_mean.vector().axpy(1, WSS_mean_avg.vector())

            # Reset tmp solutions
            for index_avg in index_dict_cycle.values():
                index_avg.vector().zero()

            WSS_mean_avg.vector().zero()

            counters_to_save.pop(0)

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Saving hemodynamic indices", "=" * 10)

    # Time average computed indices
    if average_over_cycles!=0:
        n=n_cycles
    else:
        n=(counter - start) // step
        
    if len(cycles_to_average) != 0:
        index_dict = index_dict
    else:
        index_dict=index_dict_cycle
        
    if len(cycles_to_average) != 0:
        WSS_mean = WSS_mean
    else:
        WSS_mean = WSS_mean_avg

    print("n=", n)
    print("index_dict=", index_dict)
    print("WSS_mean=", WSS_mean)
    

    index_dict['TWSSG'].vector()[:] = index_dict['TWSSG'].vector()[:] / n
    index_dict['TAWSS'].vector()[:] = index_dict['TAWSS'].vector()[:] / n
    WSS_mean.vector()[:] = WSS_mean.vector()[:] / n
    wss_mean = project(inner(WSS_mean, WSS_mean) ** (1 / 2), U_b1)
    wss_mean_vec = wss_mean.vector().get_local()
    tawss_vec = index_dict['TAWSS'].vector().get_local()

    # Compute RRT, OSI, and ECAP based on mean and absolute WSS
    index_dict['RRT'].vector().set_local(1 / wss_mean_vec)
    index_dict['OSI'].vector().set_local(0.5 * (1 - wss_mean_vec / tawss_vec))
    index_dict['ECAP'].vector().set_local(index_dict['OSI'].vector().get_local() / tawss_vec)

    for index in ['RRT', 'OSI', 'ECAP']:
        index_dict[index].vector().apply("insert")

    # Rename displayed variable names
    for name, var in index_dict.items():
        var.rename(name, name)

    # Write indices to file
    for name, index in index_dict.items():
        indices[name].write(index)

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Post processing finished", "=" * 10)
        print("Results saved to: {}".format(folder))



def get_dabla_function():
    """                                                                                                                   
    Compiles a string in C++ and expose as a Python object (dabla),                                                       
    used to compute several hemodynamic quantities.                                                                       
                                                                                                                          
    Returns:                                                                                                              
        dabla: C++ compiled function                                                                                      
    """

    cpp_code = """                                                                                                        
    #include <pybind11/pybind11.h>                                                                                        
    #include <dolfin.h>                                                                                                   
    namespace dolfin                                                                                                      
    {                                                                                                                     
        void dabla(dolfin::GenericVector& a, dolfin::GenericVector& b) {                                                  
            for (unsigned int i=0; i < b.size(); i++) {                                                                   
                b.setitem(i, pow((pow(a[i], 2) + pow(a[b.size() + i], 2) + pow(a[2 * b.size() + i], 2) ), 0.5));          
            }                                                                                                             
        }                                                                                                                 
    }                                                                                                                     
    PYBIND11_MODULE(SIGNATURE, m)                                                                                         
    {                                                                                                                     
        m.def("dabla", &dolfin::dabla);                                                                                   
    }                                                                                                                     
    """

    dabla = compile_cpp_code(cpp_code).dabla



if __name__ == '__main__':
    folder, nu, rheology_model, rho, dt, velocity_degree, _, T, save_frequency, start_cycle, step, average_over_cycles = read_command_line()
    compute_hemodynamic_indices(folder, nu, rheology_model, rho, dt, velocity_degree, T, save_frequency, start_cycle, step, average_over_cycles)
