import shutil
import sys
from os import path, remove

sys.path.append("..")
sys.path.append("../automatedPostProcessing")

from automatedPostProcessing.compute_flow_and_simulation_metrics import compute_flow_and_simulation_metrics


def test_compute_flow_and_simulation_metrics_for_full_cycle():
    # Path to test results and params
    results_path, nu, dt, velocity_degree, T, save_frequency, start_cycle, step, flow_metrics_path = get_default_parameters()
    average_over_cycles = False
    times_to_average = []

    # Run post-processing
    compute_flow_and_simulation_metrics(results_path, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step, average_over_cycles)

    # Check that output folder exists
    assert path.exists(flow_metrics_path) and path.isdir(flow_metrics_path)

    # Check that output files exist
    metric_names = ["u_time_avg", "l_plus", "t_plus", "CFL", "strain", "length_scale", "time_scale",
                    "velocity_scale", "characteristic_edge_length", "dissipation", "kinetic_energy",
                    "turbulent_kinetic_energy", "turbulent_dissipation"]

    for name in metric_names:
        xdmf_path = path.join(flow_metrics_path, "{}.xdmf".format(name))
        h5_path = path.join(flow_metrics_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0

    # Check that averaged velocity file exists
    u_avg_path = path.join(results_path, "u_avg.h5")
    assert path.exists(u_avg_path)
    assert path.getsize(u_avg_path) > 0

    # Remove generated files
    shutil.rmtree(flow_metrics_path)
    remove(u_avg_path)


def test_compute_flow_and_simulation_metrics_at_one_instance():
    # Path to test results and params
    results_path, nu, dt, velocity_degree, T, save_frequency, start_cycle, step, flow_metrics_path = get_default_parameters()

    average_over_cycles = False

    # Average over specific time
    time_in_ms = 100 * dt
    times_to_average = [time_in_ms, 2 * time_in_ms, 4 * time_in_ms]

    # Run post-processing
    compute_flow_and_simulation_metrics(results_path, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step, average_over_cycles)

    # Check that output folder exists
    assert path.exists(flow_metrics_path) and path.isdir(flow_metrics_path)

    # Check that output files exist
    metric_names = ["u_time_avg", "l_plus", "t_plus", "CFL", "strain", "length_scale", "time_scale",
                    "velocity_scale", "characteristic_edge_length", "dissipation", "kinetic_energy",
                    "turbulent_kinetic_energy", "turbulent_dissipation"]
    for time in times_to_average:
        for name in metric_names:
            xdmf_path = path.join(flow_metrics_path, "{}_{}.xdmf".format(name, time))
            h5_path = path.join(flow_metrics_path, "{}_{}.h5".format(name, time))
            assert path.exists(xdmf_path) and path.exists(h5_path)
            assert path.getsize(xdmf_path) > 0
            assert path.getsize(h5_path) > 0

    # Check that averaged velocity file exists
    u_avg_path = path.join(results_path, "u_avg.h5")
    assert path.exists(u_avg_path)
    assert path.getsize(u_avg_path) > 0

    # Remove generated files
    shutil.rmtree(flow_metrics_path)
    remove(u_avg_path)


def test_compute_flow_and_simulation_metrics_averaged_over_one_cycle():
    # Path to test results and params
    results_path, nu, dt, velocity_degree, T, save_frequency, start_cycle, step, flow_metrics_path = get_default_parameters()
    times_to_average = []

    # Average over cycle 1
    average_over_cycles = True
    cycle = 1

    # Run post-processing
    compute_flow_and_simulation_metrics(results_path, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step, average_over_cycles)

    # Check that output folder exists
    assert path.exists(flow_metrics_path) and path.isdir(flow_metrics_path)

    # Check that cycle averaged output exist
    metric_names = ["u_time_avg", "l_plus", "t_plus", "CFL", "strain", "length_scale", "time_scale",
                    "velocity_scale", "characteristic_edge_length", "dissipation", "kinetic_energy",
                    "turbulent_kinetic_energy", "turbulent_dissipation"]

    for name in metric_names:
        xdmf_path = path.join(flow_metrics_path, "{}_cycle_{:02d}.xdmf".format(name, cycle))
        h5_path = path.join(flow_metrics_path, "{}_cycle_{:02d}.h5".format(name, cycle))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0

    # Check that time averaged output exists
    for name in metric_names:
        xdmf_path = path.join(flow_metrics_path, "{}.xdmf".format(name))
        h5_path = path.join(flow_metrics_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0

    # Check that averaged velocity file exists
    u_avg_path = path.join(results_path, "u_avg.h5")
    assert path.exists(u_avg_path)
    assert path.getsize(u_avg_path) > 0

    # Remove generated files
    shutil.rmtree(flow_metrics_path)
    remove(u_avg_path)


def get_default_parameters():
    # Define default parameters for tests
    results_path = "test_results/1/Solutions"
    flow_metrics_path = "test_results/1/FlowMetrics"
    nu = 3.3018e-3
    dt = 0.951
    velocity_degree = 1
    save_frequency = 50
    start_cycle = 1
    T = 951
    step = 1

    return results_path, nu, dt, velocity_degree, T, save_frequency, start_cycle, step, flow_metrics_path


if __name__ == "__main__":
    test_compute_flow_and_simulation_metrics_for_full_cycle()
    test_compute_flow_and_simulation_metrics_at_one_instance()
    test_compute_flow_and_simulation_metrics_averaged_over_one_cycle()
