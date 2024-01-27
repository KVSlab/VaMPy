from os import path, listdir

from dolfin import parameters, FunctionSpace, XDMFFile, MPI, VectorFunctionSpace, HDF5File, Mesh, Function
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
    folders = [path.join(folder, f) for f in listdir(folder) if "SolutionsFull_" in f]
    mesh_path = path.join(folders[0], "mesh.h5")
    file_us = [HDF5File(MPI.comm_world, path.join(f, "u.h5"), "r") for f in folders]
    file_ps = [HDF5File(MPI.comm_world, path.join(f, "p.h5"), "r") for f in folders]
    file_ds = [HDF5File(MPI.comm_world, path.join(f, "d.h5"), "r") for f in folders]

    # Define HDF5Files
    dataset_us = []
    dataset_ps = []
    dataset_ds = []
    counters = []
    for i in range(len(file_us)):
        file_u = file_us[i]
        file_p = file_ps[i]
        file_d = file_ds[i]

        # Read in datasets
        dataset_u = get_dataset_names(file_u, step=step, vector_filename="/velocity/vector_%d")
        dataset_p = get_dataset_names(file_p, step=step, vector_filename="/pressure/vector_%d")
        dataset_d = get_dataset_names(file_d, step=step, vector_filename="/deformation/vector_%d")
        counters.append(len(dataset_u))

        dataset_us += dataset_u
        dataset_ps += dataset_p
        dataset_ds += dataset_d

    # file_path_mesh = "/Users/henriakj/PhD/Code/VaMPy/src/vampy/simulation/MyMesh/mesh.h5"

    # Read mesh saved as HDF5 format
    mesh = Mesh()
    with HDF5File(MPI.comm_world, mesh_path, "r") as mesh_file:
        mesh_file.read(mesh, "mesh", True)

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
    k = 0
    for i in range(len(dataset_us)):
        # Set physical time (in [ms])
        t = dt * counter * save_frequency
        file_u, file_p, file_d = file_us[k], file_ps[k], file_ds[k]
        if counter in counters:
            k += 1

        file_u.read(u, dataset_us[i])
        file_p.read(p, dataset_ps[i])

        if MPI.rank(MPI.comm_world) == 0:
            timestamp = file_u.attributes(dataset_us[i])["timestamp"]
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
            file_d.read(d, dataset_ds[i])
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
