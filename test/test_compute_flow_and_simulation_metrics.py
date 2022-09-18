import sys
from os import path

sys.path.append("..")
sys.path.append("../automatedPostProcessing")

from automatedPostProcessing.compute_flow_and_simulation_metrics import compute_flow_and_simulation_metrics


def test_compute_flow_and_simulation_metrics():
    # Path to test results and params
    results_path = "test_results/1/Solutions"
    flow_metrics_path = "test_results/1/flow_metrics"
    nu = 3.3018e-3
    dt = 0.0951

    # Run post-processing
    compute_flow_and_simulation_metrics(results_path, nu, dt, velocity_degree=1)

    # Check that output folder exists
    assert path.exists(flow_metrics_path) and path.isdir(flow_metrics_path)

    # Check that output files exist
    metric_names = ["u_mean", "l_plus", "t_plus", "CFL", "strain", "length_scale", "time_scale", "velocity_scale",
                    "characteristic_edge_length", "dissipation", "kinetic_energy", "turbulent_kinetic_energy",
                    "turbulent_dissipation"]

    for name in metric_names:
        xdmf_path = path.join(flow_metrics_path, "{}.xdmf".format(name))
        h5_path = path.join(flow_metrics_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0



if __name__ == "__main__":
    test_compute_flow_and_simulation_metrics()
