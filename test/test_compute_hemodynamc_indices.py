import sys
from os import path

sys.path.append("..")
sys.path.append("../automatedPostProcessing")

from automatedPostProcessing.compute_hemodynamic_indices import compute_hemodynamic_indices


def test_compute_hemodynamic_indices():
    # Path to test results and params
    results_path = "test_results/1/Solutions"
    indices_path = "test_results/1/Hemodynamics"
    nu = 3.3018e-3
    rho = 1060
    T = 951
    dt = 0.951
    velocity_degree = 1
    save_frequency = 50
    start_cycle = 1
    step = 1
    average_over_cycles = False

    # Run post-processing
    compute_hemodynamic_indices(results_path, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle, step,
                                average_over_cycles)

    # Check that output folder exists
    assert path.exists(indices_path) and path.isdir(indices_path)

    # Check that output files exist
    metric_names = ["ECAP", "OSI", "RRT", "TAWSS", "TWSSG", "WSS"]

    for name in metric_names:
        xdmf_path = path.join(indices_path, "{}.xdmf".format(name))
        h5_path = path.join(indices_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0


def test_compute_hemodynamic_indices_averaged_over_one_cycle():
    # Path to test results and params
    results_path = "test_results/1/Solutions"
    indices_path = "test_results/1/Hemodynamics"
    nu = 3.3018e-3
    rho = 1060
    T = 951
    dt = 0.951
    velocity_degree = 1
    save_frequency = 50
    start_cycle = 1
    step = 1
    average_over_cycles = True

    # Average over cycle 1
    cycle = 1

    # Run post-processing
    compute_hemodynamic_indices(results_path, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle, step,
                                average_over_cycles)

    # Check that output folder exists
    assert path.exists(indices_path) and path.isdir(indices_path)

    # Check that output files exist
    metric_names = ["ECAP", "OSI", "RRT", "TAWSS", "TWSSG"]

    for name in metric_names:
        xdmf_path = path.join(indices_path, "{}_cycle_{}.xdmf".format(name, cycle))
        h5_path = path.join(indices_path, "{}_cycle_{}.h5".format(name, cycle))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0


if __name__ == "__main__":
    test_compute_hemodynamic_indices()
    test_compute_hemodynamic_indices_averaged_over_one_cycle()
