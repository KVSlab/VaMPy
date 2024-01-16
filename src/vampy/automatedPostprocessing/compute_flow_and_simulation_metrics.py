from os import path, listdir

import numpy as np
from dolfin import FunctionSpace, Function, MPI, VectorFunctionSpace, Timer, project, HDF5File, XDMFFile, \
    assign, CellDiameter, Mesh, TestFunction, list_timings, TimingClear, TimingType
from vampy.automatedPostprocessing.postprocessing_common import read_command_line, get_dataset_names


def compute_flow_and_simulation_metrics(folder, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step, average_over_cycles):
    """
    Computes several flow field characteristics and metrics for the velocity field.
    Reads in the .h5 soution stored in the 'Solution' folder within the results' folder,
    and stores the post-processing files in the `FlowMetrics` folder.

    Args:
        folder (str): The folder path where the data files are located.
        nu (float): The kinematic viscosity value.
        dt (float): The time step size.
        velocity_degree (int): The degree of the velocity function space.
        T (float): The total simulation time.
        times_to_average (list): A list of time points to perform phase averaging. Needs to be in the inverval [0,T).
        save_frequency (int): The frequency of saving data during the simulation.
        start_cycle (int): The starting cycle number for averaging.
        step (int): The step size for extracting data.
        average_over_cycles (bool): A flag indicating whether to perform cycle averaging.
    """
    # File paths
    folders = [path.join(folder, f) for f in sorted(listdir(folder)) if "SolutionsFull_" in f]
    file_us = [HDF5File(MPI.comm_world, path.join(f, "u.h5"), "r") for f in folders]
    mesh_path = path.join(folder, "mesh.h5")
    file_path_u_avg = path.join(folder, "u_avg.h5")

    dataset_us = []
    file_counters = []
    saved_time_steps_per_cycle = int(T / dt / save_frequency / step)
    print(f"-- Reading from folders: {folders}")
    for i in range(len(file_us)):
        file_u = file_us[i]
        dataset_u = get_dataset_names(file_u, step=step, vector_filename="/velocity/vector_%d")

        if "0021" in folder and i == 0:
            dataset_u_sliced = dataset_u[-250:]
        else:
            slice_id = len(dataset_u) % saved_time_steps_per_cycle
            if slice_id != 0:
                dataset_u_sliced = dataset_u[:-slice_id]
            else:
                dataset_u_sliced = dataset_u

        # Add to collective dataset
        dataset_us += dataset_u_sliced
        file_counters += [i] * len(dataset_u_sliced)

    # Get names of data to extract
    if MPI.rank(MPI.comm_world) == 0:
        print("Reading dataset names")

    # Extract specific time steps if phase averaging
    number_of_cycles = int(len(dataset_us) / saved_time_steps_per_cycle)
    cycles_to_average = None

    # Get mesh information
    mesh = Mesh()
    with HDF5File(MPI.comm_world, mesh_path, "r") as mesh_file:
        mesh_file.read(mesh, "mesh", False)

    # Function space
    DG = FunctionSpace(mesh, 'DG', 0)
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)
    Vv = FunctionSpace(mesh, "CG", velocity_degree)

    if MPI.rank(MPI.comm_world) == 0:
        print("Computing average velocity: u_avg(x,t)")

    # Define u_avg(x,t)
    u = Function(V)
    u_avg = Function(V)

    # Read velocity and compute cycle averaged velocity
    if not path.exists(file_path_u_avg):
        compute_u_avg(dataset_us, file_counters, file_us, file_path_u_avg, number_of_cycles, saved_time_steps_per_cycle,
                      start_cycle, u, u_avg)

    # Perform phase averaging (Average over cycles at specified time point(s))
    if len(times_to_average) != 0:
        dataset_dict, dataset_dict_avg, number_of_files = \
            get_files_for_phase_averaging(times_to_average, dt, save_frequency, step, number_of_cycles, start_cycle,
                                          saved_time_steps_per_cycle, dataset_us)
    # Perform cycle averaging (Average per cycle) and time averaging (Average of all cycles)
    else:
        dataset_dict, dataset_dict_avg, number_of_files = \
            get_files_for_cycle_averaging(dataset_us, file_counters, file_us, file_path_u_avg, number_of_cycles,
                                          saved_time_steps_per_cycle,
                                          start_cycle)
        if average_over_cycles:
            cycles_to_average = [cycle + 1 for cycle in range(number_of_cycles)]

    # Compute flow and simulation metrics
    for time_to_average, dataset in dataset_dict.items():
        if len(times_to_average) != 0 and MPI.rank(MPI.comm_world) == 0:
            print(f"Phase averaging results over {number_of_files} cycles at t={time_to_average} ms")

        define_functions_and_iterate_dataset(folder, file_counters, file_us, time_to_average, dataset,
                                             dataset_dict_avg[time_to_average], dt,
                                             file_path_u_avg, mesh, nu, number_of_files, DG, V, Vv,
                                             cycles_to_average, saved_time_steps_per_cycle)


# Function to reshape the array
def reshape_array(arr, rows, cols):
    if len(arr) != rows * cols:
        raise ValueError("The size of the original array does not match the specified dimensions.")
    return [arr[i::cols] for i in range(cols)]


def compute_u_avg(dataset_us, file_counters, file_us, file_path_u_avg, n_cycles, saved_time_steps_per_cycle,
                  start_cycle, u, u_avg):
    """
    Iterate over saved time steps and compute average velocity based on save sampling parameters

    Args:
        dataset_names (dict): Contains timestep and respective names of solution timesteps
        file_path_u_avg (str): File path to average velocity solution
        file_u (str): File path to velocity solution
        n_cycles (int): Number of cardiac cycles
        saved_time_steps_per_cycle (int): Determines number of saved time steps per cycle
        start_cycle (int): Determines which cardiac cycle to start sampling from
        u (Function): Function storing velocity
        u_avg (Function): Function for storing average velocity
        start_cycle (int): Determines which cardiac cycle to start from for post-processing
    """
    files_u = reshape_array(dataset_us, n_cycles, saved_time_steps_per_cycle)
    files_counter = reshape_array(file_counters, n_cycles, saved_time_steps_per_cycle)

    for files, counters in zip(files_u, files_counter):
        time = -1
        for k, data in zip(counters, files):
            file_u = file_us[k]

            if time == -1:
                time = int(file_u.attributes(data)["timestamp"])

            # Accumulate velocity
            file_u.read(u, data)
            u_avg.vector().axpy(1, u.vector())

        # Average over pre-defined amount of cycles
        u_avg.vector()[:] /= (n_cycles - start_cycle + 1)

        if MPI.rank(MPI.comm_world) == 0:
            print("=" * 10, f"Computing average velocity at time: {time} ms", "=" * 10)

        # Save average velocity to HDF5 format
        file_mode = "w" if not path.exists(file_path_u_avg) else "a"
        viz_u_avg = HDF5File(MPI.comm_world, file_path_u_avg, file_mode=file_mode)
        viz_u_avg.write(u_avg, "/velocity", time)
        viz_u_avg.close()

        # Reset u_avg vector
        u_avg.vector().zero()


def get_files_for_cycle_averaging(dataset_us, file_counters, file_us, file_path_u_avg, number_of_cycles,
                                  saved_time_steps_per_cycle, start_cycle):
    """
    Retrieves the dataset dictionaries for cycle averaging.

    Args:
        dataset_names (list): A list of dataset names.
        file_path_u_avg (str): The file path for averaging the dataset.
        number_of_cycles (int): The total number of cycles.
        saved_time_steps_per_cycle (int): The number of time steps saved per cycle.
        start_cycle (int): The starting cycle number for averaging.

    Returns:
        dataset_dict (dict): A dictionary containing dataset names for cycle averaging.
        dataset_dict_avg (dict): A dictionary containing averaged dataset names.
        number_of_files (int): The total number of files.
    """
    # Create a dataset with dataset names from the starting index
    id_start = (start_cycle - 1) * saved_time_steps_per_cycle
    dataset_dict = {"": dataset_us[id_start:]}

    # Create similar dataset for the mean velocity
    file_u_avg = HDF5File(MPI.comm_world, file_path_u_avg, "r")
    dataset_dict_avg = {"": get_dataset_names(file_u_avg, start=0, step=1) * (number_of_cycles - start_cycle + 1)}

    # Get number of files based on the sliced dataset names
    number_of_files = len(dataset_dict[""])

    return dataset_dict, dataset_dict_avg, number_of_files


def get_files_for_phase_averaging(times_to_average, dt, save_frequency, step, number_of_cycles, start_cycle,
                                  saved_time_steps_per_cycle, dataset_names):
    """
    Generates dataset dictionaries for phase averaging based on the selected times.

    Args:
        times_to_average (list): A list of time points to perform phase averaging.
        dt (float): The time step size.
        save_frequency (int): The frequency of saving data during the simulation.
        step (int): The step size for extracting data.
        number_of_cycles (int): The total number of cycles.
        start_cycle (int): The starting cycle number for averaging.
        saved_time_steps_per_cycle (int): The number of time steps saved per cycle.
        dataset_names (list): A list of dataset names.

    Returns:
        dataset_dict (dict): A dictionary mapping time points to corresponding dataset names.
        dataset_dict_avg (dict): A dictionary mapping time points to averaged dataset names.
        number_of_files (int): The total number of files.
    """
    dataset_dict = {}
    dataset_dict_avg = {}

    # Iterate over selected times to average over
    for t in times_to_average:
        # Find corresponding IDs for dataset_name list
        time_id = int(t / dt / save_frequency / step)
        time_ids = [time_id + saved_time_steps_per_cycle * k for k in range(number_of_cycles)][start_cycle - 1:]

        # Get dataset_names at corresponding times
        dataset_dict[f"_{t}"] = [dataset_names[max(0, j - 1)] for j in time_ids]
        dataset_dict_avg[f"_{t}"] = [dataset_names[max(0, time_id - 1)]] * len(time_ids)

    # Get the number of files for the last time point
    number_of_files = len(dataset_dict[f"_{times_to_average[-1]}"])

    return dataset_dict, dataset_dict_avg, number_of_files


def define_functions_and_iterate_dataset(folder, file_counters, file_us, time_to_average, dataset, dataset_avg, dt,
                                         file_path_u_avg, mesh, nu, number_of_files, DG, V, Vv, cycles_to_average,
                                         saved_time_steps_per_cycle):
    """
    Defines functions and vector functions for all metrics to be computed, iterates through dataset and computes
    several flow and simulation metrics.
    After iteration, saves the  metrics to .xdmf file format.

    Args:
        time_to_average (int): Time in ms to average over
        number_of_files (int): Number of dataset files to iterate over
        DG (FunctionSpace): Discontinous Galerkin function space
        V (VectorFunctionSpace): Continous Galerkin function space for vector solutions
        Vv (FunctionSpace): Continous Galerkin function space for scalar solutions
        mesh (Mesh): Volumetric mesh of computational domain
        dataset (dict): Contains velocity solution dict keys
        dataset_avg (dict): Contains average velocity solution dict keys
        file_path_u (str): File path to velocity solution
        file_path_u_avg (str): File path to average velocity solution
        file_u (str): File path to velocity solution
        saved_time_steps_per_cycle (int): Determines number of saved time steps per cycle
        folder (str): Path to simulation results
        nu (float): Viscosity
        dt (float): Time step in [ms]
        cycles_to_average (list): List of which cycles to average
    """
    # Functions for storing values
    v = TestFunction(DG)
    u = Function(V)
    u_avg = Function(V)
    u_time_avg = Function(V)
    u_time_cycle_avg = Function(V)
    u_prime = Function(V)

    # Plus-values
    l_plus_avg = Function(DG)
    l_plus_cycle_avg = Function(DG)
    l_plus = Function(DG)
    t_plus_avg = Function(DG)
    t_plus_cycle_avg = Function(DG)
    t_plus = Function(DG)

    # Kolmogorov scales
    length_scale = Function(DG)
    length_scale_avg = Function(DG)
    length_scale_cycle_avg = Function(DG)
    time_scale = Function(DG)
    time_scale_avg = Function(DG)
    time_scale_cycle_avg = Function(DG)
    velocity_scale = Function(DG)
    velocity_scale_avg = Function(DG)
    velocity_scale_cycle_avg = Function(DG)

    # Inner grad(u), grad(u)
    turbulent_dissipation = Function(DG)
    turbulent_dissipation_avg = Function(DG)
    turbulent_dissipation_cycle_avg = Function(DG)
    strain = Function(DG)
    strain_avg = Function(DG)
    strain_cycle_avg = Function(DG)
    dissipation = Function(DG)
    dissipation_avg = Function(DG)
    dissipation_cycle_avg = Function(DG)

    # Energy
    kinetic_energy = Function(Vv)
    kinetic_energy_avg = Function(Vv)
    kinetic_energy_cycle_avg = Function(Vv)
    turbulent_kinetic_energy = Function(Vv)
    turbulent_kinetic_energy_avg = Function(Vv)
    turbulent_kinetic_energy_cycle_avg = Function(Vv)

    # Velocity
    u0 = Function(Vv)
    u1 = Function(Vv)
    u2 = Function(Vv)
    u0_prime = Function(Vv)
    u1_prime = Function(Vv)
    u2_prime = Function(Vv)

    # CFL
    CFL = Function(DG)
    CFL_avg = Function(DG)
    CFL_cycle_avg = Function(DG)

    # Characteristic edge length
    h = CellDiameter(mesh)
    characteristic_edge_length = project(h, DG)

    # Create XDMF files for saving metrics
    if cycles_to_average is None:
        save_path = path.join(folder, "FlowMetrics", "%s{}%s.xdmf".format(time_to_average))
        cycles_to_average = []
    else:
        save_path = path.join(folder, "FlowMetrics", "%s%s.xdmf")
        number_of_files = len(cycles_to_average)

    metric_names = ["kinetic_energy", "turbulent_kinetic_energy"]

    metric_variables_cycle_avg = [kinetic_energy_cycle_avg, turbulent_kinetic_energy_cycle_avg]

    metric_variables_avg = [kinetic_energy_avg, turbulent_kinetic_energy_avg]

    metric_dict_cycle = dict(zip(metric_names, metric_variables_cycle_avg))
    metric_dict = dict(zip(metric_names, metric_variables_avg))

    counters_to_save = [saved_time_steps_per_cycle * cycle for cycle in cycles_to_average]
    cycle_names = [""] + ["_cycle_{:02d}".format(cycle) for cycle in cycles_to_average]
    metrics = {}
    for cycle_name in cycle_names:
        for vn in metric_dict.keys():
            metrics[vn + cycle_name] = XDMFFile(MPI.comm_world, save_path % (vn, cycle_name))
            metrics[vn + cycle_name].parameters["rewrite_function_mesh"] = False
            metrics[vn + cycle_name].parameters["flush_output"] = True

    # Get u average
    file_u_avg = HDF5File(MPI.comm_world, file_path_u_avg, "r")

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Starting post processing", "=" * 10)

    counter = 0
    for k, data, data_avg in zip(file_counters, dataset, dataset_avg):
        counter += 1

        file_u = file_us[k]

        # Read velocity and cycle averaged velocity
        file_u_avg.read(u_avg, data_avg)
        file_u.read(u, data)

        if MPI.rank(MPI.comm_world) == 0:
            time = file_u.attributes(data)["timestamp"]
            print("=" * 10, f"Time: {time} ms", "=" * 10)

        # Compute u_prime
        t0 = Timer("u prime")
        u_prime.vector().set_local(u.vector().get_local() - u_avg.vector().get_local())
        u_prime.vector().apply("insert")
        t0.stop()

        t0 = Timer("Kinetic energy")
        assign(u0, u.sub(0))
        assign(u1, u.sub(1))

        if mesh.geometry().dim() == 3:
            assign(u2, u.sub(2))

        kinetic_energy.vector().set_local(
            0.5 * (u0.vector().get_local() ** 2 + u1.vector().get_local() ** 2 + u2.vector().get_local() ** 2))
        kinetic_energy.vector().apply("insert")
        kinetic_energy_cycle_avg.vector().axpy(1, kinetic_energy.vector())
        t0.stop()

        t0 = Timer("Turbulent kinetic energy")
        assign(u0_prime, u_prime.sub(0))
        assign(u1_prime, u_prime.sub(1))

        if mesh.geometry().dim() == 3:
            assign(u2_prime, u_prime.sub(2))

        turbulent_kinetic_energy.vector().set_local(
            0.5 * (u0_prime.vector().get_local() ** 2
                   + u1_prime.vector().get_local() ** 2
                   + u2_prime.vector().get_local() ** 2))
        turbulent_kinetic_energy.vector().apply("insert")
        turbulent_kinetic_energy_cycle_avg.vector().axpy(1, turbulent_kinetic_energy.vector())
        t0.stop()

        if counter % 10 == 0:
            list_timings(TimingClear.clear, [TimingType.wall])

        if len(cycles_to_average) != 0 and counter == counters_to_save[0]:
            # Get cycle number
            cycle = int(counters_to_save[0] / saved_time_steps_per_cycle)
            if MPI.rank(MPI.comm_world) == 0:
                print("========== Storing cardiac cycle {} ==========".format(cycle))

            # Get average over sampled time steps
            for metric in list(metric_dict_cycle.values()):
                metric.vector()[:] = metric.vector()[:] / saved_time_steps_per_cycle

            # Store solution
            for name, metric in metric_dict_cycle.items():
                metrics[name + "_cycle_{:02d}".format(cycle)].write_checkpoint(metric, name)

            # Append solution to total solution
            for metric_avg, metric_cycle_avg in zip(list(metric_dict.values()), list(metric_dict_cycle.values())):
                metric_cycle_avg.vector().apply("insert")
                metric_avg.vector().axpy(1, metric_cycle_avg.vector())

            # Reset tmp solutions
            for metric_cycle_avg in list(metric_dict_cycle.values()):
                metric_cycle_avg.vector().zero()

            counters_to_save.pop(0)

    # Get average over sampled time steps
    metrics_dict_to_save = metric_dict if len(cycles_to_average) != 0 else metric_dict_cycle
    if len(cycles_to_average) != 0:
        number_of_files = len(cycles_to_average)

    for metric in metrics_dict_to_save.values():
        metric.vector()[:] = metric.vector()[:] / number_of_files

    # Store average data
    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Saving flow and simulation metrics", "=" * 10)

    for name, metric in metrics_dict_to_save.items():
        metrics[name].write_checkpoint(metric, name)

    # Print summary info
    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Flow and simulation metrics summary", "=" * 10)

    for metric_name, metric_value in metrics_dict_to_save.items():
        sum_ = MPI.sum(MPI.comm_world, np.sum(metric_value.vector().get_local()))
        num = MPI.sum(MPI.comm_world, metric_value.vector().get_local().shape[0])
        mean = sum_ / num
        max_ = MPI.max(MPI.comm_world, metric_value.vector().get_local().max())
        min_ = MPI.min(MPI.comm_world, metric_value.vector().get_local().min())
        if MPI.rank(MPI.comm_world) == 0:
            print(metric_name, "mean:", mean)
            print(metric_name, "max:", max_)
            print(metric_name, "min:", min_)

    if MPI.rank(MPI.comm_world) == 0:
        print("=" * 10, "Post processing finished", "=" * 10)
        save_folder = save_path.rsplit("/", 1)[0]
        print(f"Results saved to: {save_folder}")


def main_metrics():
    folder, nu, _, dt, velocity_degree, _, _, T, save_frequency, times_to_average, start_cycle, step, \
    average_over_cycles = read_command_line()

    compute_flow_and_simulation_metrics(folder, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step, average_over_cycles)


if __name__ == '__main__':
    folder, nu, _, dt, velocity_degree, _, _, T, save_frequency, times_to_average, start_cycle, step, \
    average_over_cycles = read_command_line()

    compute_flow_and_simulation_metrics(folder, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step, average_over_cycles)
