import sys
from os import path

sys.path.append("..")
sys.path.append("../automatedPostProcessing")

from automatedPostProcessing.compute_flow_and_simulation_metrics import compute_flow_and_simulation_metrics


def test_compute_flow_and_simulation_metrics_for_full_cycle():
    # Path to test results and params
    results_path = "test_results/1/Solutions"
    flow_metrics_path = "test_results/1/FlowMetrics"
    nu = 3.3018e-3
    dt = 0.0951
    T = 951
    velocity_degree = 1
    save_frequency = 5
    start_cycle = 1
    times_to_average = []
    step = 1

    # Run post-processing
    compute_flow_and_simulation_metrics(results_path, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step)

    # Check that output folder exists
    assert path.exists(flow_metrics_path) and path.isdir(flow_metrics_path)

    # Check that output files exist
    metric_names = ["time_averaged_u", "l_plus", "t_plus", "CFL", "strain", "length_scale", "time_scale",
                    "velocity_scale", "characteristic_edge_length", "dissipation", "kinetic_energy",
                    "turbulent_kinetic_energy", "turbulent_dissipation"]

    for name in metric_names:
        xdmf_path = path.join(flow_metrics_path, "{}.xdmf".format(name))
        h5_path = path.join(flow_metrics_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0


def test_compute_flow_and_simulation_metrics_at_one_instance():
    # Path to test results and params
    results_path = "test_results/1/Solutions"
    flow_metrics_path = "test_results/1/FlowMetrics"
    nu = 3.3018e-3
    dt = 0.0951
    velocity_degree = 1
    save_frequency = 5
    start_cycle = 1
    T = 10 * dt
    time = 0.4755
    times_to_average = [time]
    step = 1

    # Run post-processing
    compute_flow_and_simulation_metrics(results_path, nu, dt, velocity_degree, T, times_to_average, save_frequency,
                                        start_cycle, step)

    # Check that output folder exists
    assert path.exists(flow_metrics_path) and path.isdir(flow_metrics_path)

    # Check that output files exist
    metric_names = ["time_averaged_u", "l_plus", "t_plus", "CFL", "strain", "length_scale", "time_scale",
                    "velocity_scale", "characteristic_edge_length", "dissipation", "kinetic_energy",
                    "turbulent_kinetic_energy", "turbulent_dissipation"]

    for name in metric_names:
        xdmf_path = path.join(flow_metrics_path, "{}_{}.xdmf".format(name, time))
        h5_path = path.join(flow_metrics_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0


if __name__ == "__main__":
    test_compute_flow_and_simulation_metrics_for_full_cycle()
    test_compute_flow_and_simulation_metrics_at_one_instance()
