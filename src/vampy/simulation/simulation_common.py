from os import makedirs, path

import numpy as np
from dolfin import MPI, Constant, HDF5File, Measure, assemble, assign


def get_file_paths(folder):
    """
    Create folder where data and solutions (velocity, mesh, pressure) is stored

    Args:
        folder (str): Path to data storage location

    Returns:
        files (dict): Contains filepaths for respective solution files
    """
    common_path = path.join(folder, "Solutions")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path):
            makedirs(common_path)

    file_p = path.join(common_path, "p.h5")
    file_u = path.join(common_path, "u.h5")
    file_u_mean = path.join(common_path, "u_mean.h5")
    file_mesh = path.join(common_path, "mesh.h5")
    files = {"u": file_u, "u_mean": file_u_mean, "p": file_p, "mesh": file_mesh}

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
        print(
            "X range: {} to {} (delta: {:.4f})".format(
                min(xmin), max(xmax), max(xmax) - min(xmin)
            )
        )
        print(
            "Y range: {} to {} (delta: {:.4f})".format(
                min(ymin), max(ymax), max(ymax) - min(ymin)
            )
        )
        print(
            "Z range: {} to {} (delta: {:.4f})".format(
                min(zmin), max(zmax), max(zmax) - min(zmin)
            )
        )
        print("Number of cells: {}".format(sum(num_cells)))
        print("Number of cells per processor: {}".format(int(np.mean(num_cells))))
        print("Number of edges: {}".format(sum(num_edges)))
        print("Number of faces: {}".format(sum(num_faces)))
        print("Number of facets: {}".format(sum(num_facets)))
        print("Number of vertices: {}".format(sum(num_vertices)))
        print("Volume: {:.4f}".format(volume))
        print("Number of cells per volume: {:.4f}".format(sum(num_cells) / volume))


def store_u_mean(
    T,
    dt,
    save_solution_at_tstep,
    save_solution_frequency,
    u_mean,
    u_mean0,
    u_mean1,
    u_mean2,
    NS_parameters,
):
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
    files = NS_parameters["files"]
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
