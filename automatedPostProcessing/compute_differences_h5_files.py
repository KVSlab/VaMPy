from __future__ import print_function

from pathlib import Path

from dolfin import *
from os import getcwd, makedirs, path
from argparse import ArgumentParser
import numpy as np

try:
    parameters["reorder_dofs_serial"] = False
except NameError:
    pass

try:
    from dolfin import *

    try:
        parameters["allow_extrapolation"] = True
    except NameError:
        pass
except ImportError:
    pass


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--case', type=str, default="simulation/results_folder/model/data/1/Solutions",
                        help="Path to simulation results (PostProc folder for executing compute_hemodynamic_indices.py and h5_files folder for executing compute_differences_h5_files.py): mesh.h5, u.h5, p.h5, nu.h5", metavar="PATH")    
    parser.add_argument('--dt', type=float, default=0.0951, help="Time step of simulation, measured in [ms]")
    parser.add_argument('--velocity-degree', type=int, default=1, help="Degree of velocity element")
    parser.add_argument('--casenonN', type=str, default="simulation/results_folder/model/data/1/Solutions",
                        help="Path to non-N simulation results (h5_files forder for executing compute_differences_h5_files.py: WSS.h5, TAWSS.h5, OSI.h5, RRT.h5, ECAP.h5, TWSSG.h5)", metavar="PATH")

    args = parser.parse_args()

    return args.case, args.dt, args.velocity_degree, args.casenonN

def compute_differences_hemodynamic_indices(case_path_N, case_path_nonN, dt, velocity_degree, hi, function_space_for_hi):
    """
    Loads velocity fields from completed CFD simulation,
    and computes and saves the following hemodynamic quantities:
    (1) WSS - Wall shear stress
    (2) TAWSS - Time averaged wall shear stress
    (3) TWSSG - Temporal wall shear stress gradient
    (4) OSI - Oscillatory shear index
    (5) RRT - Relative residence time
    (6) ECAP - Endothelial cell activation potential

    Args:
        velocity_degree (int): Finite element degree of velocity
        case_path_N (Path): Path to results from Newtonian results
        case_path_nonN (Path): Path to results from non-Newtonian results
        dt (float): Time step of simulation
    """
    # File paths Newtonian
    original_case_path_N = case_path_N
    case_path_N = Path(case_path_N)
    #file_path_u = case_path_N / "u.h5"
    #mesh_path = case_path_N / "mesh.h5"
    common_path_xdmf_N = path.join(case_path_N, "h5_files")
    case_path_N = Path(common_path_xdmf_N)
    #file_path_WSS = case_path_N / "WSS.h5"
    file_path_hi = case_path_N / "{}.h5".format(str(hi))

    # File paths non-Newtonian
    case_path_nonN = Path(case_path_nonN)
    #file_path_u_nonN = case_path_nonN / "u.h5"
    common_path_xdmf_nonN = path.join(case_path_nonN, "h5_files")
    case_path_nonN = Path(common_path_xdmf_nonN)
    # file_path_WSS_nonN = case_path_nonN / "WSS.h5"
    file_path_hi_nonN = case_path_nonN / "{}.h5".format(str(hi))
    #------------------------------------------------------------------------------------------------------------
    # Create folder to store h5 files
    common_path_diff = path.join(original_case_path_N, "Diff_h5_files")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path_diff):
            makedirs(common_path_diff)
    common_path_diff = Path(common_path_diff)
    #------------------------------------------------------------------------------------------------------------
    """# WSS_mean
    WSS_mean = Function(V_b1)
    wss_mean = Function(function_space_for_hi)"""

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Start post processing of {}".format(str(hi)), "=" * 10)

    # HI Newtonian case
    x_N = Function(function_space_for_hi) # It is where I save the variables
    f = HDF5File(MPI.comm_world, file_path_hi.__str__(), "r")
    f.read(x_N, "/{}".format(str(hi)))

    # HI non-Newtonian case
    x_nonN = Function(function_space_for_hi) # It is where I save the variables
    f_nonN = HDF5File(MPI.comm_world, file_path_hi_nonN.__str__(), "r")
    f_nonN.read(x_nonN, "/{}".format(str(hi)))

    DIFF = Function(function_space_for_hi) # It is where I save the variables
    diff = Function(function_space_for_hi) # Output file 

    #DIFF.vector()[:] = np.absolute( (x_nonN.vector()[:] - x_N.vector()[:]) * 100 / (x_N.vector()[:] + 1E-13) )
    DIFF.vector()[:] = np.abs( (x_nonN.vector().get_local() - x_N.vector().get_local()) * 100 / (x_N.vector().get_local() + 1E-13) )
    DIFF.rename("{}".format(str(hi)), "{}".format(str(hi)))

    # Save difference between Newtonian and non-Newtonian results
    diff_path = (common_path_diff  / "{}.xdmf".format(str(hi))).__str__()
    diff      = XDMFFile(MPI.comm_world, diff_path)

    for f in [diff]:
        f.parameters["flush_output"] = True
        f.parameters["functions_share_mesh"] = True
        f.parameters["rewrite_function_mesh"] = False

    diff.write(DIFF)
    print("========== Post processing finished ==========")
    print("Results saved to: {}".format(common_path_diff))


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

    return dabla

def compute_differences_instantaneous_WSS(case_path_N, case_path_nonN, dt, velocity_degree, function_space_for_WSS):
    """
    (1) WSS - Wall shear stress Newtonian
    (2) WSS - Wall shear stress non-Newtonian
    Args:
        velocity_degree (int): Finite element degree of velocity
        case_path_N (Path): Path to results from Newtonian results
        case_path_nonN (Path): Path to results from non-Newtonian results
        dt (float): Time step of simulation
    """
    # File paths Newtonian
    original_case_path_N = case_path_N
    case_path_N = Path(case_path_N)
    common_path_xdmf_N = path.join(case_path_N, "h5_files")
    case_path_N = Path(common_path_xdmf_N)
    print("path_N = ", case_path_N)
    file_path_WSS = case_path_N / "WSS.h5"

    # File paths non-Newtonian
    case_path_nonN = Path(case_path_nonN)
    common_path_xdmf_nonN = path.join(case_path_nonN, "h5_files")
    case_path_nonN = Path(common_path_xdmf_nonN)
    print("path_nonN = ", case_path_nonN)
    file_path_WSS_nonN = case_path_nonN / "WSS.h5"
    #------------------------------------------------------------------------------------------------------------
    # Create folder to store h5 files
    common_path_diff = path.join(original_case_path_N, "Diff_h5_files")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path_diff):
            makedirs(common_path_diff)
    common_path_diff = Path(common_path_diff)
    #------------------------------------------------------------------------------------------------------------
    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Start post processing of instantaneous WSS", "=" * 10)
    #------------------------------------------------------------------------------------------------------------
    # WSS_mean Newtonian case
    WSS_mean = Function(function_space_for_WSS)
    # WSS_mean non-Newtonian case
    WSS_mean_nonN = Function(function_space_for_WSS)
    # WSS_mean non-Newtonian case
    DIFF = Function(function_space_for_WSS) # It is where I save the variables
    diff = Function(function_space_for_WSS) # Output file 
    #------------------------------------------------------------------------------------------------------------
    # Create writer for saving the difference between Newtonian and non-Newtonian results
    diff_path = (common_path_diff  / "diff_WSS.xdmf").__str__()
    diff      = XDMFFile(MPI.comm_world, diff_path)

    diff.parameters["flush_output"] = True
    diff.parameters["functions_share_mesh"] = True
    diff.parameters["rewrite_function_mesh"] = True
    #------------------------------------------------------------------------------------------------------------
    # Start post-processing from 2nd cycle using every 10th time step, or 2000 time steps per cycle
    start = 0  # save_data = 5 -> 10000 / 5 = 2000
    step = 1  # save_data = 5 ->    10 / 5 = 2

    file_counter = start
    while True:
        try:
            # Read WSS Newtonian case
            f = HDF5File(MPI.comm_world, file_path_WSS.__str__(), "r")
            vec_name = "/instantaneous_WSS/vector_%d" % file_counter
            tstep = f.attributes(vec_name)["timestamp"]
            print("=" * 10, "Timestep Newtonian: {}".format(tstep), "=" * 10)
            f.read(WSS_mean, vec_name)
            # Read WSS non-Newtonian case
            f_nonN = HDF5File(MPI.comm_world, file_path_WSS_nonN.__str__(), "r")
            vec_name_nonN = "/instantaneous_WSS/vector_%d" % file_counter
            tstep_nonN = f_nonN.attributes(vec_name_nonN)["timestamp"]
            print("=" * 10, "Timestep non-Newtonian: {}".format(tstep_nonN), "=" * 10)
            f.read(WSS_mean_nonN, vec_name_nonN)
        except:
            print("=" * 10, "Finished reading solutions", "=" * 10)
            break

        # Save instantaneous differences of WSS as XDMF file
        #DIFF.vector()[:] = np.absolute( (WSS_mean_nonN.vector()[:] - WSS_mean.vector()[:]) * 100 / (WSS_mean.vector()[:] + 1E-13) )
        DIFF.vector()[:] = np.abs( (WSS_mean_nonN.vector().get_local() - WSS_mean.vector().get_local() ) * 100 / (WSS_mean.vector().get_local() + 1E-13) )
        #DIFF.vector()[:] = (WSS_mean.vector()[:] - WSS_mean.vector()[:]) * 100 / (WSS_mean.vector()[:] + 1E-5)
        #DIFF.vector()[:] = WSS_mean_nonN.vector()[:]
        #DIFF.vector()[:] = WSS_mean.vector()[:]
        print("WSS_min=", WSS_mean.vector().get_local().min())
        print("WSS_min_nonN=", WSS_mean_nonN.vector().get_local().min())
        print("WSS_max=", WSS_mean.vector().get_local().max())
        print("WSS_max_nonN=", WSS_mean_nonN.vector().get_local().max())
        DIFF.rename("WSS", "WSS")
        diff.write(DIFF, tstep*dt)
        #------------------------------------------------------------------------------------------------------------
        # Update file_counter
        file_counter += step
        print("file_counter=", file_counter)

    print("=" * 10, "Saving differences of WSS", "=" * 10)
    n = (file_counter - start) // step
    print("n=", n)

    print("========== Post processing finished ==========")
    print("Results saved to: {}".format(common_path_diff))


if __name__ == '__main__':
    folder, dt, velocity_degree, folder_nonN = read_command_line()

    # File paths Newtonian
    original_case_path_N = folder
    folder = Path(folder)
    #file_path_u = case_path_N / "u.h5"
    mesh_path = folder / "mesh.h5"

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
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)

    if MPI.rank(MPI.comm_world) == 0:
        print("Define functions")
    u = Function(V)


    """HI = ["TAWSS", "OSI", "RRT", "ECAP", "TWSSG"]
    for hi in HI:
        compute_differences_hemodynamic_indices(folder, folder_nonN, dt, velocity_degree, hi, U_b1)"""
        
    # There is some problem with the computation of the instantaneous differences of WSS!!! 
    # It shows that the difference is 0, but TAWSS shows that there are big differences!
    compute_differences_instantaneous_WSS(original_case_path_N, folder_nonN, dt, velocity_degree, V_b1)
