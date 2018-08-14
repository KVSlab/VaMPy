from __future__ import print_function
import os
from time import time
from hashlib import sha1
from dolfin import *
from os import listdir, path, sep, system
from stress import STRESS
import sys
import numpy as np

parameters["reorder_dofs_serial"] = False

#def load_velocity(u, files_u0, files_u1, files_u2, i):
#    if len(files_u1) == 0:
#	u.vector()[:] = Function(files_u0[i])


def get_dataset_names(f, num_files=3000000, step=1, start=1, print_=True,
                      name="velocity%s"):
    check = True

    # Find start file
    t0 = time()
    while check:
        if f.has_dataset(name % start):
            check = False
            start -= step
        start += step

    # Get names
    names = []
    for i in range(num_files):
        step = 1
        index = start + i*step
        if f.has_dataset(name % index):
            names.append(name % index)
    t1 = time()

    if MPI.rank(mpi_comm_world()) == 0 and print_:
        print("")
        print("="*6 + " Timesteps to average over " + "="*6)
        print("Length on data set names:", len(names))
        print("Start index:", start)
        print("Wanted num files:", start)
        print("Step between files:", step)
        print("Time used:", t1 - t0, "s")

    return names

def main(case_path):
    file_path_x = path.join(case_path, "u0.h5")
    file_path_y = path.join(case_path, "u1.h5")
    file_path_z = path.join(case_path, "u2.h5")

    f_0 = HDF5File(mpi_comm_world(), file_path_x, "r")
    f_1 = HDF5File(mpi_comm_world(), file_path_y, "r")
    f_2 = HDF5File(mpi_comm_world(), file_path_z, "r")

    start = 10000
    if MPI.rank(mpi_comm_world()) == 0:
        print("The postprocessing start from", start)
    dataset_names = get_dataset_names(f_0 , start=start)[::10]

    mesh = Mesh()
    f_0.read(mesh, "Mesh", False)
    dt = 0.951 # (time2 - time1) / (timestep2 - timestep1)

    # Load mesh
    bm = BoundaryMesh(mesh, 'exterior')
    f = File(path.join(case_path, "Boundary_mesh.pvd"))
    f << bm
    del f

    if MPI.rank(mpi_comm_world()) == 0:
        print("Define spaces")
    uorder = 1
    V_b1 = VectorFunctionSpace(bm, "CG", 1)
    U_b1 = FunctionSpace(bm, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", uorder)
    U = FunctionSpace(mesh, "CG", uorder)

    if MPI.rank(mpi_comm_world()) == 0:
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
namespace dolfin {
    void dabla(dolfin::GenericVector& a, dolfin::GenericVector& b) {
        for (unsigned int i=0; i < b.size(); i++) {
            b.setitem(i, pow( (pow(a[i],2) + pow(a[b.size()+i],2) +pow(a[2*b.size()+i],2) ) ,0.5 ));
       }
   }
}
"""
    func = getattr(compile_extension_module(cpp_code), 'dabla')

    if MPI.rank(mpi_comm_world()) == 0:
        print("Start 'simulation'")
    for data in dataset_names:
        if MPI.rank(mpi_comm_world()) == 0:
            print("Timestep", data[8:])

        # Get u
        f_0.read(u0, data)
        f_1.read(u1, data)
        f_2.read(u2, data)
        assign(u.sub(0), u0)
        assign(u.sub(1), u1)
        assign(u.sub(2), u2)

        # Compute WSS
        if MPI.rank(mpi_comm_world()) == 0:
            print("Compute WSS")
        tau = stress(u)
        tau.vector()[:] = tau.vector()[:] * 1000
        WSS.vector().axpy(1, tau.vector()[:])

        if MPI.rank(mpi_comm_world()) == 0:
            print("Compute OSI")
        func(tau.vector(), osi_.vector())
        OSI.vector().axpy(1, osi_.vector())

        # Compute TWSSG
        if MPI.rank(mpi_comm_world()) == 0:
            print("Compute TWSSG")
        twssg.vector().set_local((tau.vector().array() - tau_prev.vector().array()) / dt)
        twssg.vector().apply("insert")
        func(twssg.vector(), twssg_.vector())
        TWSSG.vector().axpy(1, twssg_.vector())

        # Update tau
        if MPI.rank(mpi_comm_world()) == 0:
            print("Update WSS")
        tau_prev.vector().zero()
        tau_prev.vector().axpy(1, tau.vector())


    TWSSG.vector()[:] = TWSSG.vector()[:] / len(dataset_names)
    WSS.vector()[:] = WSS.vector()[:] / len(dataset_names)
    OSI.vector()[:] = OSI.vector()[:] / len(dataset_names)
    WSS_new.vector()[:] = OSI.vector()[:]

    try:
        func(WSS.vector(), rrt_.vector())
        rrt_arr = rrt_.vector().array()
        rrt_arr[rrt_arr.nonzero()[0]] = 1e-6
        rrt_arr[np.nan(rrt_arr)] = 1e-6
        RRT.vector().apply("insert")

        OSI_arr = OSI.vector().array()
        OSI_arr[OSI_arr.nonzero()[0]] = 1e-6
        OSI_arr[np.isnan(OSI_arr)] = 1e-6
        OSI.vector().set_local(0.5*(1 - rrt_arr / OSI_arr))
        OSI.vector().apply("insert")
        save = True
    except:
        print("Fail for OSI and RRT")
        save = False

    if save:
        osi_file = File(path.join(case_path, "OSI.xml.gz"))
        osi_file << OSI
        del osi_file

        osi_file = File(path.join(case_path, "OSI.pvd"))
        osi_file << OSI
        del osi_file

        rrt_file = File(path.join(case_path, "RRT.xml.gz"))
        rrt_file << RRT
        del rrt_file

        rrt_file = File(path.join(case_path, "RRT.pvd"))
        rrt_file << RRT
        del rrt_file

    twssg_file = File(path.join(case_path, "TWSSG.xml.gz"))
    twssg_file << TWSSG
    del twssg_file

    twssg_file = File(path.join(case_path, "TWSSG.pvd"))
    twssg_file << TWSSG
    del twssg_file

    wss_file = File(path.join(case_path, "WSS.xml.gz"))
    wss_file << WSS_new
    del wss_file

    wss_file = File(path.join(case_path, "WSS.pvd"))
    wss_file << WSS_new
    del wss_file

if __name__ == '__main__':
    main(sys.argv[1])
