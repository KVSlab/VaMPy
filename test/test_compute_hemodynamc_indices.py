import sys
from os import path

sys.path.append("..")
sys.path.append("../automatedPostProcessing")

from automatedPostProcessing.compute_hemodynamic_indices import compute_hemodynamic_indices


def test_compute_hemodynamic_indices():
    # Path to test results and params
    results_path = "test_results/1/Solutions"
    nu = 3.3018e-3
    rho = 1060
    dt = 0.951
    velocity_degree = 1

    # Run post-processing
    compute_hemodynamic_indices(results_path, nu, rho, dt, velocity_degree)

    # Check that output files exist
    metric_names = ["ECAP", "OSI", "RRT", "TAWSS", "TWSSG", "WSS"]

    for name in metric_names:
        xdmf_path = path.join(results_path, "{}.xdmf".format(name))
        h5_path = path.join(results_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0


if __name__ == "__main__":
    test_compute_hemodynamic_indices()
