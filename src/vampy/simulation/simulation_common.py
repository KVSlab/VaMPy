from os import path, makedirs

import numpy as np
from dolfin import MPI, assemble, Constant, assign, HDF5File, Measure


def get_file_paths(folder, additional_variables=[]):
    """
    Create folder where data and solutions (velocity, mesh, pressure) is stored

    Args:
        folder (str): Path to data storage location
        additional_variables (list): List of additional variables to store

    Returns:
        files (dict): Contains filepaths for respective solution files
    """
    common_path_half = path.join(folder, "Solutions")
    common_path_full = path.join(folder, "SolutionsFull")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path_half):
            makedirs(common_path_half)
        if not path.exists(common_path_full):
            makedirs(common_path_full)
    variables = ["p", "u", "u_mean", "mesh"] + additional_variables
    files = {
        'half': {},
        'full': {}
    }

    for variable in variables:
        files['half'][variable] = path.join(common_path_half, f"{variable}.h5")
        files['full'][variable] = path.join(common_path_full, f"{variable}.h5")

    return files


def print_mesh_information(mesh):
    """
    Print geometric information about the volumetric mesh

    Args:
        mesh (dolfin.cpp.mesh.Mesh): Volumetric mesh
    """
    dx = Measure("dx", domain=mesh)
    comm = MPI.comm_world
    local_xmin = mesh.coordinates()[:, 0].min()
    local_xmax = mesh.coordinates()[:, 0].max()
    local_ymin = mesh.coordinates()[:, 1].min()
    local_ymax = mesh.coordinates()[:, 1].max()
    local_zmin = mesh.coordinates()[:, 2].min()
    local_zmax = mesh.coordinates()[:, 2].max()
    xmin = comm.gather(local_xmin, 0)
    xmax = comm.gather(local_xmax, 0)
    ymin = comm.gather(local_ymin, 0)
    ymax = comm.gather(local_ymax, 0)
    zmin = comm.gather(local_zmin, 0)
    zmax = comm.gather(local_zmax, 0)

    local_num_cells = mesh.num_cells()
    local_num_edges = mesh.num_edges()
    local_num_faces = mesh.num_faces()
    local_num_facets = mesh.num_facets()
    local_num_vertices = mesh.num_vertices()
    num_cells = comm.gather(local_num_cells, 0)
    num_edges = comm.gather(local_num_edges, 0)
    num_faces = comm.gather(local_num_faces, 0)
    num_facets = comm.gather(local_num_facets, 0)
    num_vertices = comm.gather(local_num_vertices, 0)
    volume = assemble(Constant(1) * dx)

    if MPI.rank(MPI.comm_world) == 0:
        print("=== Mesh information ===")
        print(f"X range: {min(xmin)} to {max(xmax)} (delta: {max(xmax) - min(xmin):.4f})")
        print(f"Y range: {min(ymin)} to {max(ymax)} (delta: {max(ymax) - min(ymin):.4f})")
        print(f"Z range: {min(zmin)} to {max(zmax)} (delta: {max(zmax) - min(zmin):.4f})")
        print(f"Number of cells: {sum(num_cells)}")
        print(f"Number of cells per processor: {int(np.mean(num_cells))}")
        print(f"Number of edges: {sum(num_edges)}")
        print(f"Number of faces: {sum(num_faces)}")
        print(f"Number of facets: {sum(num_facets)}")
        print(f"Number of vertices: {sum(num_vertices)}")
        print(f"Volume: {volume:.4f}")
        print(f"Number of cells per volume: {sum(num_cells) / volume:.4f}")


def store_u_mean(T, dt, save_solution_at_tstep, save_solution_frequency, u_mean, u_mean0, u_mean1, u_mean2,
                 NS_parameters):
    """
    Store the time averaged velocity into file u_mean

    Args:
        T (float): End time
        dt (float): Time step size
        save_solution_at_tstep (int): Time step when solution is to be saved
        save_solution_frequency (int): Frequency of solution saving
        u_mean (Function): Vector function for storing vector solution of average velocity
        u_mean0 (Function): Function storing x-component
        u_mean1 (Function): Function storing y-component
        u_mean2 (Function): Function storing z-component
    """
    # Get the file path
    files = NS_parameters['files']
    u_mean_path = files["u_mean"]

    # Divide the accumulated velocity by the number of steps
    NumTStepForAverage = (T / dt - save_solution_at_tstep) / save_solution_frequency + 1
    u_mean0.vector()[:] = u_mean0.vector()[:] / NumTStepForAverage
    u_mean1.vector()[:] = u_mean1.vector()[:] / NumTStepForAverage
    u_mean2.vector()[:] = u_mean2.vector()[:] / NumTStepForAverage
    assign(u_mean.sub(0), u_mean0)
    assign(u_mean.sub(1), u_mean1)
    assign(u_mean.sub(2), u_mean2)

    # Save u_mean
    with HDF5File(MPI.comm_world, u_mean_path, "w") as u_mean_file:
        u_mean_file.write(u_mean, "u_mean")


def store_velocity_and_pressure_h5(files, U, p_, tstep, u_, u_mean0, u_mean1, u_mean2, D=None, du_=None, blood=None):
    """
    Store the velocity and pressure values to an HDF5 file.

    Args:
        NS_parameters (dict): A dictionary containing the parameters for Navier-Stokes equations.
        U (Function): A vector function space to assign the velocity components.
        p_ (Function): The pressure function.
        tstep (int): The current time step.
        u_ (List[Function]): A list containing the velocity components.
        u_mean0 (Function): The accumulated x-component of the velocity.
        u_mean1 (Function): The accumulated y-component of the velocity.
        u_mean2 (Function): The accumulated z-component of the velocity.
        du_ (List[Function]): List of deformation components
        D (Function): A vector function space to assign the deformation components.

    Returns:
        None
    """
    # Assign velocity components to vector solution
    for i in range(3):
        assign(U.sub(i), u_[i])

    # Get save paths
    p_path = files['p']
    u_path = files['u']
    brt_path = files['brt']
    file_mode = "w" if not path.exists(p_path) else "a"

    # Save pressure
    with HDF5File(MPI.comm_world, p_path, file_mode=file_mode) as viz_p:
        viz_p.write(p_, "/pressure", tstep)

    # Save velocity
    with HDF5File(MPI.comm_world, u_path, file_mode=file_mode) as viz_u:
        viz_u.write(U, "/velocity", tstep)

    # Save residence time
    with HDF5File(MPI.comm_world, brt_path, file_mode=file_mode) as viz_brt:
        viz_brt.write(blood, "/blood", tstep)

    # Accumulate velocity
    u_mean0.vector().axpy(1, u_[0].vector())
    u_mean1.vector().axpy(1, u_[1].vector())
    u_mean2.vector().axpy(1, u_[2].vector())

    # Save deformation if present
    if D is not None and du_ is not None:
        for i in range(3):
            assign(D.sub(i), du_[i])

        # Save path to deformation
        d_path = files['d']

        # Save deformation
        with HDF5File(MPI.comm_world, d_path, file_mode=file_mode) as viz_d:
            viz_d.write(D, "/deformation", tstep)


def dump_probes(eval_dict, newfolder, tstep):
    """
    Save variables along the centerline for CFD simulation diagnostics and light-weight post-processing.

    Args:
        eval_dict (dict): Dictionary with probe arrays for velocity components and pressure along the centerline.
        newfolder (str): The path to the folder where the probe data will be saved.
        tstep (float): The current time step.

    Returns:
        None
    """
    # Construct the file path for probes
    filepath = path.join(newfolder, "Probes")

    # Ensure the directory exists on the master process
    if MPI.rank(MPI.comm_world) == 0 and not path.exists(filepath):
        makedirs(filepath)

    # Extract probe arrays for each variable
    variables = ["centerline_u_x_probes", "centerline_u_y_probes", "centerline_u_z_probes", "centerline_p_probes"]
    arrs = {var: eval_dict[var].array() for var in variables}

    # Dump stats on the master process
    if MPI.rank(MPI.comm_world) == 0:
        filenames = {
            "centerline_u_x_probes": f"u_x_{tstep}.probes",
            "centerline_u_y_probes": f"u_y_{tstep}.probes",
            "centerline_u_z_probes": f"u_z_{tstep}.probes",
            "centerline_p_probes": f"p_{tstep}.probes"
        }
        for var, arr in arrs.items():
            arr.dump(path.join(filepath, filenames[var]))

    # Ensure all processes have dumped data before clearing
    MPI.barrier(MPI.comm_world)

    # Clear stats
    for var in variables:
        eval_dict[var].clear()
