from __future__ import print_function

from pathlib import Path

import numpy as np
from dolfin import *

from postprocessing_common import STRESS, read_command_line

parameters["reorder_dofs_serial"] = False


def compute_wss(case_path, nu, dt):
    """
    Loads velocity fields from completed CFD simulation,
    and computes and saves the following hemodynamic quantities:
    (1) WSS - Wall shear stress
    (2) TWSSG - Temporal wall shear stress gradient
    (3) OSI - Oscillatory shear index
    (4) RRT - Relative residence time

    Args:
        case_path (Path): Path to results from simulation
        nu (float): Viscosity
        dt (float): Time step of simulation
    """
    # File paths
    file_path_x = case_path / "u0.h5"
    file_path_y = case_path / "u1.h5"
    file_path_z = case_path / "u2.h5"
    mesh_path = case_path / "mesh.xdmf" #generated from serial

    # Read headers in HDF5 files
    #f_0 = XDMFFile(MPI.comm_world, file_path_x.__str__())
    #f_1 = XDMFFile(MPI.comm_world, file_path_y.__str__())
    #f_2 = XDMFFile(MPI.comm_world, file_path_z.__str__())

    # Start post-processing from 2nd cycle using every 10th time step, or 2000 time steps per cycle
    start = 0  # save_data = 5 -> 10000 / 5 = 2000
    step = 2  # save_data = 5 ->    10 / 5 = 2

    # Read mesh
    mesh = Mesh()
    with XDMFFile(MPI.comm_world, mesh_path.__str__()) as mesh_file:
        mesh_file.read(mesh)

    # Load mesh
    bm = BoundaryMesh(mesh, 'exterior')

    if MPI.rank(MPI.comm_world) == 0:
        print("Define function spaces")
    velocity_degree = 1
    V_b1 = VectorFunctionSpace(bm, "CG", 1)
    U_b1 = FunctionSpace(bm, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)
    U = FunctionSpace(mesh, "CG", velocity_degree)

    if MPI.rank(MPI.comm_world) == 0:
        print("Define functions")
    u = Function(V)
    u0 = Function(U)
    u1 = Function(U)
    u2 = Function(U)

    # RRT
    rrt_ = Function(U_b1)
    RRT = Function(U_b1)

    # WSS
    WSS = Function(V_b1)

    # OSI
    OSI = Function(U_b1)
    WSS_new = Function(U_b1)
    osi_ = Function(U_b1)

    # TWSSG
    TWSSG = Function(U_b1)
    twssg_ = Function(U_b1)
    twssg = Function(V_b1)
    tau_prev = Function(V_b1)

    stress = STRESS(u, 0.0, nu, mesh)
    dabla = get_dabla_function()

    if MPI.rank(MPI.comm_world) == 0:
        print("Start 'simulation'")

    file_counter = start
    #for file_counter in range(0,11):
    while True:
        if MPI.rank(MPI.comm_world) == 0:
            print("Timestep", file_counter * 5)

        # Read velocity components to respective functions and assign them to vector function u
        try:
            #f1 = HDF5File(mesh.mpi_comm(), file_path_x.__str__(), "r")
            #f2 = HDF5File(mesh.mpi_comm(), file_path_y.__str__(), "r")
            #f3 = HDF5File(mesh.mpi_comm(), file_path_z.__str__(), "r")
            f1 = HDF5File(MPI.comm_world, file_path_x.__str__(), "r")
            f2 = HDF5File(MPI.comm_world, file_path_y.__str__(), "r")
            f3 = HDF5File(MPI.comm_world, file_path_z.__str__(), "r")
            vec_name = "/velocity/vector_%d"%file_counter
            timestamp = f1.attributes(vec_name)["timestamp"]
            print("Timestep",timestamp)
            f1.read(u0,vec_name)
            f2.read(u1,vec_name)
            f3.read(u2,vec_name)
        except:
            print("Failed to read in velocity at file_counter={}".format(file_counter))
            break
        assign(u.sub(0), u0)
        assign(u.sub(1), u1)
        assign(u.sub(2), u2)

        # Compute WSS
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS")
        tau = stress()
        tau.vector()[:] = tau.vector()[:] * 1000
        WSS.vector().axpy(1, tau.vector())

        if MPI.rank(MPI.comm_world) == 0:
            print("Compute OSI")
        dabla(tau.vector(), osi_.vector())
        OSI.vector().axpy(1, osi_.vector())

        # Compute TWSSG
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute TWSSG")
        twssg.vector().set_local((tau.vector().get_local() - tau_prev.vector().get_local()) / dt)
        twssg.vector().apply("insert")
        dabla(twssg.vector(), twssg_.vector())
        TWSSG.vector().axpy(1, twssg_.vector())

        # Update tau
        if MPI.rank(MPI.comm_world) == 0:
            print("Update WSS")
        tau_prev.vector().zero()
        tau_prev.vector().axpy(1, tau.vector())

        # Update file_counter
        file_counter += step

    n = (file_counter - start) // step
    TWSSG.vector()[:] = TWSSG.vector()[:] / n
    WSS.vector()[:] = WSS.vector()[:] / n
    OSI.vector()[:] = OSI.vector()[:] / n
    WSS_new.vector()[:] = OSI.vector()[:]

    #try:
        #func(WSS.vector(), rrt_.vector())
        #rrt_arr = rrt_.vector().get_local()
        #rrt_arr[rrt_arr.nonzero()[0]] = 1e-6
        #rrt_arr[np.isnan(rrt_arr)] = 1e-6
        #RRT.vector().apply("insert")

        #OSI_arr = OSI.vector().get_local()
        #OSI_arr[OSI_arr.nonzero()[0]] = 1e-6
        #OSI_arr[np.isnan(OSI_arr)] = 1e-6
        #OSI.vector().set_local(0.5*(1 - rrt_arr / OSI_arr))
        #OSI.vector().apply("insert")
        #save = True
    #except:
        #print("Fail for OSI and RRT")
        #save = False

    #if save:
        #osi_file = File((case_path / "OSI.xml.gz").__str__())
        #osi_file << OSI
        #del osi_file

        #osi_file = File((case_path / "OSI.pvd").__str__())
        #osi_file << OSI
        #del osi_file

        #rrt_file = File((case_path / "RRT.xml.gz").__str__())
        #rrt_file << RRT
        #del rrt_file

        #rrt_file = File((case_path / "RRT.pvd").__str__())
        #rrt_file << RRT
        #del rrt_file

    twssg_file = File((case_path / "TWSSG.xml.gz").__str__())
    twssg_file << TWSSG
    del twssg_file

    twssg_file = File((case_path / "TWSSG.pvd").__str__())
    twssg_file << TWSSG
    del twssg_file

    wss_file = File((case_path / "WSS.xml.gz").__str__())
    wss_file << WSS_new
    del wss_file

    wss_file = File((case_path / "WSS.pvd").__str__())
    wss_file << WSS_new
    del wss_file

    value2 = (case_path / "TWSSG1.xdmf").__str__()
    twssg1 = XDMFFile(MPI.comm_world, value2)
    twssg1.parameters["flush_output"] = True
    twssg1.parameters["functions_share_mesh"] = True
    twssg1.parameters["rewrite_function_mesh"] = False
    twssg1.write(TWSSG)


    value1 = (case_path / "WSS1.xdmf").__str__()
    wss1 = XDMFFile(MPI.comm_world, value1)
    wss1.parameters["flush_output"] = True
    wss1.parameters["functions_share_mesh"] = True
    wss1.parameters["rewrite_function_mesh"] = False
    wss1.write(WSS_new)



def get_dabla_function():
    """
    Compiles a string in C++ and expose as a Python object (dabla),
    used the compute several hemodynamic quantities.

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

    return dabla


if __name__ == '__main__':
    folder, nu, dt = read_command_line()
    compute_wss(folder, nu, dt)
