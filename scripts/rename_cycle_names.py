from os import path, listdir

from dolfin import MPI, HDF5File
from vampy.automatedPostprocessing.postprocessing_common import get_dataset_names


def main(folder, step=2):
    folders_sorted = [path.join(folder, f) for f in sorted(listdir(folder)) if "SolutionsFull_" in f]
    folders_unsorted = [path.join(folder, f) for f in sorted(folder) if "SolutionsFull_" in f]

    file_us = [HDF5File(MPI.comm_world, path.join(f, "u.h5"), "r") for f in folders_unsorted]
    mesh_path = path.join(folder, "mesh.h5")
    file_path_u_avg = path.join(folder, "u_avg.h5")

    saved_time_steps_per_cycle = 250
    for i in range(len(file_us)):
        file_u = file_us[i]
        dataset_u = get_dataset_names(file_u, step=step, vector_filename="/velocity/vector_%d")
        slice_id = len(dataset_u) % saved_time_steps_per_cycle
        if slice_id != 0:
            dataset_u_sliced = dataset_u[:-slice_id]
        else:
            dataset_u_sliced = dataset_u


if __name__ == '__main__':
    cases = ["0003"]
    conditions = ["SR", "AF"][-1:]

    for case in cases:
        for condition in conditions:
            folder = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1"
            main(folder)
