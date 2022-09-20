from __future__ import print_function

from pathlib import Path

from dolfin import *

from postprocessing_common import read_command_line

try:
    parameters["reorder_dofs_serial"] = False
except NameError:
    pass


def compute_velocity_and_pressure(case_path, dt, velocity_degree, pressure_degree):
    """
    Loads velocity and pressure from compressed .h5 CFD solution and
    converts and saves to .xdmf format for visualization (in e.g. ParaView).

    Args:
        case_path (str): Path to results from simulation
        dt (float): Time step of simulation
        velocity_degree (int): Finite element degree of velocity
        pressure_degree (int): Finite element degree of pressure
    """
    # File paths
    case_path = Path(case_path)
    file_path_u = case_path / "u.h5"
    file_path_p = case_path / "p.h5"
    mesh_path = case_path / "mesh.h5"

    # Read mesh saved as HDF5 format
    mesh = Mesh()
    with HDF5File(MPI.comm_world, mesh_path.__str__(), "r") as mesh_file:
        mesh_file.read(mesh, "mesh", False)

    # Define functionspaces and functions
    if MPI.rank(MPI.comm_world) == 0:
        print("Define function spaces")
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)
    Q = FunctionSpace(mesh, "CG", pressure_degree)

    if MPI.rank(MPI.comm_world) == 0:
        print("Define functions")

    u = Function(V)
    p = Function(Q)

    # Create writer for velocity and pressure
    u_path = (case_path / "velocity.xdmf").__str__()
    p_path = (case_path / "pressure.xdmf").__str__()

    u_writer = XDMFFile(MPI.comm_world, u_path)
    p_writer = XDMFFile(MPI.comm_world, p_path)

    for writer in [u_writer, p_writer]:
        writer.parameters["flush_output"] = True
        writer.parameters["functions_share_mesh"] = True
        writer.parameters["rewrite_function_mesh"] = False

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Start post processing", "=" * 10)

    file_counter = 0
    while True:
        # Read in velocity solution to vector function u and pressure to function p
        try:
            u_h5 = HDF5File(MPI.comm_world, file_path_u.__str__(), "r")
            p_h5 = HDF5File(MPI.comm_world, file_path_p.__str__(), "r")
            u_name = "/velocity/vector_%d" % file_counter
            p_name = "/pressure/vector_%d" % file_counter
            timestamp = u_h5.attributes(u_name)["timestamp"]
            print("=" * 10, "Timestep: {}".format(timestamp), "=" * 10)
            u_h5.read(u, u_name)
            p_h5.read(p, p_name)
        except:
            print("=" * 10, "Finished reading solutions", "=" * 10)
            break

        # Store velocity
        u.rename("velocity", "velocity")
        u_writer.write(u, dt * file_counter)

        # Store pressure
        p.rename("pressure", "pressure")
        p_writer.write(p, dt * file_counter)

        # Update file_counter
        file_counter += 1

    print("========== Post processing finished ==========")
    print("Results saved to: {}".format(case_path))


if __name__ == '__main__':
    folder, _, _, dt, velocity_degree, pressure_degree, _, _, _, _, _ = read_command_line()
    compute_velocity_and_pressure(folder, dt, velocity_degree, pressure_degree)
