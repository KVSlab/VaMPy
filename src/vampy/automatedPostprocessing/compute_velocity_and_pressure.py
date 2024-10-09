from os import path

from dolfin import (
    MPI,
    Function,
    FunctionSpace,
    HDF5File,
    Mesh,
    VectorFunctionSpace,
    XDMFFile,
    parameters,
)

from vampy.automatedPostprocessing.postprocessing_common import read_command_line, get_dataset_names

try:
    parameters["reorder_dofs_serial"] = False
except NameError:
    pass


def compute_velocity_and_pressure(folder, dt, save_frequency, velocity_degree, pressure_degree, step):
    """
    Loads velocity and pressure from compressed .h5 CFD solution and
    converts and saves to .xdmf format for visualization (in e.g. ParaView).

    Args:
        folder (str): Path to results from simulation
        dt (float): Time step of simulation
        save_frequency (int): Frequency that velocity and pressure has been stored
        velocity_degree (int): Finite element degree of velocity
        pressure_degree (int): Finite element degree of pressure
        step (int): Step size determining number of times data is sampled
    """
    # File paths
    file_path_u = path.join(folder, "u.h5")
    file_path_p = path.join(folder, "p.h5")
    file_path_d = path.join(folder, "d.h5")
    file_path_mesh = path.join(folder, "mesh.h5")

    # Define HDF5Files
    file_u = HDF5File(MPI.comm_world, file_path_u, "r")
    file_p = HDF5File(MPI.comm_world, file_path_p, "r")
    file_d = None
    if path.exists(file_path_d):
        file_d = HDF5File(MPI.comm_world, file_path_d, "r")

    # Read in datasets
    dataset_u = get_dataset_names(file_u, step=step, vector_filename="/velocity/vector_%d")
    dataset_p = get_dataset_names(file_p, step=step, vector_filename="/pressure/vector_%d")
    if file_d is not None:
        dataset_d = get_dataset_names(file_d, step=step, vector_filename="/deformation/vector_%d")

    # Read mesh saved as HDF5 format
    mesh = Mesh()
    with HDF5File(MPI.comm_world, file_path_mesh, "r") as mesh_file:
        mesh_file.read(mesh, "mesh", False)

    # Define functionspaces and functions
    if MPI.rank(MPI.comm_world) == 0:
        print("Define function spaces")
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)
    Q = FunctionSpace(mesh, "CG", pressure_degree)

    if MPI.rank(MPI.comm_world) == 0:
        print("Define functions")

    u = Function(V)
    d = Function(V)
    p = Function(Q)

    # Create writer for velocity and pressure
    u_path = path.join(folder, "velocity.xdmf")
    p_path = path.join(folder, "pressure.xdmf")

    u_writer = XDMFFile(MPI.comm_world, u_path)
    p_writer = XDMFFile(MPI.comm_world, p_path)

    for writer in [u_writer, p_writer]:
        writer.parameters["flush_output"] = True
        writer.parameters["functions_share_mesh"] = True
        writer.parameters["rewrite_function_mesh"] = False

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Start post processing", "=" * 10)

    counter = 1
    for i in range(len(dataset_u)):
        # Set physical time (in [ms])
        t = dt * counter * save_frequency

        file_u.read(u, dataset_u[i])
        file_p.read(p, dataset_p[i])

        if MPI.rank(MPI.comm_world) == 0:
            timestamp = file_u.attributes(dataset_u[i])["timestamp"]
            print("=" * 10, "Timestep: {}".format(timestamp), "=" * 10)

        # Store velocity
        u.rename("velocity", "velocity")
        u_writer.write(u, t)

        # Store pressure
        p.rename("pressure", "pressure")
        p_writer.write(p, t)

        # Store deformation
        # NB: Storing together with velocity.
        if file_d is not None:
            file_d.read(d, dataset_d[i])
            d.rename("deformation", "deformation")
            u_writer.write(d, t)

        # Update file_counter
        counter += step

    print("========== Post processing finished ==========")
    print("Results saved to: {}".format(folder))


def main_convert():
    folder, _, _, dt, velocity_degree, pressure_degree, _, _, save_frequency, _, _, step, _ = read_command_line()
    compute_velocity_and_pressure(folder, dt, save_frequency, velocity_degree, pressure_degree, step)


if __name__ == '__main__':
    folder, _, _, dt, velocity_degree, pressure_degree, _, _, save_frequency, _, _, step, _ = read_command_line()
    compute_velocity_and_pressure(folder, dt, save_frequency, velocity_degree, pressure_degree, step)
