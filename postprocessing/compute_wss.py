from __future__ import print_function
import os
from time import time
from hashlib import sha1
from pathlib import Path
import sys

from dolfin import *
import numpy as np

from stress import STRESS

parameters["reorder_dofs_serial"] = False


def main(case_path):
    # File paths
    file_path_x = case_path / "u0.xdmf"
    file_path_y = case_path / "u1.xdmf"
    file_path_z = case_path / "u2.xdmf"
    mesh_path = case_path / "mesh.xdmf"

    # Read headers in HDF5 files
    f_0 = XDMFFile(MPI.comm_world, file_path_x.__str__())
    f_1 = XDMFFile(MPI.comm_world, file_path_y.__str__())
    f_2 = XDMFFile(MPI.comm_world, file_path_z.__str__())

    # Start post-processing from 2nd cycle using every 10th time step, or 2000 time steps per cycle
    start = 0   # save_data = 5 -> 10000 / 5 = 2000
    step = 1    # save_data = 5 ->    10 / 5 = 2
    dt = 0.951

    # Read mesh
    mesh = Mesh()
    with XDMFFile(MPI.comm_world, mesh_path.__str__()) as mesh_file:
        mesh_file.read(mesh)

    # Load mesh
    bm = BoundaryMesh(mesh, 'exterior')
    f = File((case_path / "Boundary_mesh.pvd").__str__())
    f << bm
    del f

    if MPI.rank(MPI.comm_world) == 0:
        print("Define spaces")
    uorder = 1
    V_b1 = VectorFunctionSpace(bm, "CG", 1)
    U_b1 = FunctionSpace(bm, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", uorder)
    U = FunctionSpace(mesh, "CG", uorder)

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

    mu = 0.0035
    stress = STRESS(u, 0.0, mu, mesh)

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
    func = compile_cpp_code(cpp_code).dabla

    if MPI.rank(MPI.comm_world) == 0:
        print("Start 'simulation'")

    file_counter = start
    while True:
        if MPI.rank(MPI.comm_world) == 0:
            print("Timestep", file_counter*5)

        # Get u
        try:
            f_0.read_checkpoint(u0, "u0", file_counter)
            f_1.read_checkpoint(u1, "u1", file_counter)
            f_2.read_checkpoint(u2, "u2", file_counter)
        except:
            break
        assign(u.sub(0), u0)
        assign(u.sub(1), u1)
        assign(u.sub(2), u2)

        # Compute WSS
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute WSS")
        tau = stress(u)
        tau.vector()[:] = tau.vector()[:] * 1000
        WSS.vector().axpy(1, tau.vector())

        if MPI.rank(MPI.comm_world) == 0:
            print("Compute OSI")
        func(tau.vector(), osi_.vector())
        OSI.vector().axpy(1, osi_.vector())

        # Compute TWSSG
        if MPI.rank(MPI.comm_world) == 0:
            print("Compute TWSSG")
        twssg.vector().set_local((tau.vector().get_local() - tau_prev.vector().get_local()) / dt)
        twssg.vector().apply("insert")
        func(twssg.vector(), twssg_.vector())
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

    try:
        func(WSS.vector(), rrt_.vector())
        rrt_arr = rrt_.vector().get_local()
        rrt_arr[rrt_arr.nonzero()[0]] = 1e-6
        rrt_arr[np.isnan(rrt_arr)] = 1e-6
        RRT.vector().apply("insert")

        OSI_arr = OSI.vector().get_local()
        OSI_arr[OSI_arr.nonzero()[0]] = 1e-6
        OSI_arr[np.isnan(OSI_arr)] = 1e-6
        OSI.vector().set_local(0.5*(1 - rrt_arr / OSI_arr))
        OSI.vector().apply("insert")
        save = True
    except:
        print("Fail for OSI and RRT")
        save = False

    if save:
        osi_file = File((case_path / "OSI.xml.gz").__str__())
        osi_file << OSI
        del osi_file

        osi_file = File((case_path / "OSI.pvd").__str__())
        osi_file << OSI
        del osi_file

        rrt_file = File((case_path / "RRT.xml.gz").__str__())
        rrt_file << RRT
        del rrt_file

        rrt_file = File((case_path / "RRT.pvd").__str__())
        rrt_file << RRT
        del rrt_file

    twssg_file = File((case_path / "TWSSG.xml.gz").__str__())
    twssg_file << TWSSG
    del twssg_file

    twssg_file = File((case_path / "TWSSG.pvd").__str__())
    twssg_file << TWSSG
    del twssg_file

    wss_file = File((case_path /  "WSS.xml.gz").__str__())
    wss_file << WSS_new
    del wss_file

    wss_file = File((case_path / "WSS.pvd").__str__())
    wss_file << WSS_new
    del wss_file

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("Run program as 'python path/to/compute_wss.py path/to/results/run_number/VTK'")
        sys.exit(0)

    main(Path.cwd() / sys.argv[1])
