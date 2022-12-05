import sys
from os import path

sys.path.append("..")
sys.path.append("../automatedPostProcessing")

from automatedPostProcessing.compute_velocity_and_pressure import compute_velocity_and_pressure


def test_compute_velocity_and_pressure():
    # Path to test results and params
    results_path = "test_results/1/Solutions"
    dt = 0.951
    velocity_degree = 1
    pressure_degree = 1

    # Run post-processing
    compute_velocity_and_pressure(results_path, dt, velocity_degree, pressure_degree)

    # Check that output files exist
    metric_names = ["velocity", "pressure"]

    for name in metric_names:
        xdmf_path = path.join(results_path, "{}.xdmf".format(name))
        h5_path = path.join(results_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0


if __name__ == "__main__":
    test_compute_velocity_and_pressure()
