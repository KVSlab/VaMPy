import shutil
from os import path

from vampy.automatedPostprocessing.compute_hemodynamic_indices import compute_hemodynamic_indices


def test_compute_hemodynamic_indices():
    # Path to test results and params
    results_path, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle, step, indices_path = \
        get_default_parameters()
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

    # Remove generated files
    shutil.rmtree(indices_path)


def test_compute_hemodynamic_indices_averaged_over_one_cycle():
    # Path to test results and params
    results_path, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle, step, indices_path \
        = get_default_parameters()

    # Average over cycle 1
    average_over_cycles = True
    cycle = 1

    # Run post-processing
    compute_hemodynamic_indices(results_path, nu, rho, dt, T, velocity_degree, save_frequency,
                                start_cycle, step, average_over_cycles)

    # Check that output folder exists
    assert path.exists(indices_path) and path.isdir(indices_path)

    # Check that cycle averaged output exist
    metric_names = ["ECAP", "OSI", "RRT", "TAWSS", "TWSSG"]

    for name in metric_names:
        xdmf_path = path.join(indices_path, "{}_cycle_{:02d}.xdmf".format(name, cycle))
        h5_path = path.join(indices_path, "{}_cycle_{:02d}.h5".format(name, cycle))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0

    # Check that full cycle hemodynamic indices are still computed
    for name in metric_names:
        xdmf_path = path.join(indices_path, "{}.xdmf".format(name))
        h5_path = path.join(indices_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0

    # Remove generated files
    shutil.rmtree(indices_path)


def get_default_parameters():
    # Define default parameters for tests
    results_path = "tests/test_results/1/Solutions"
    indices_path = "tests/test_results/1/Hemodynamics"
    nu = 3.3018e-3
    rho = 1060
    dt = 0.951
    T = 951
    velocity_degree = 1
    save_frequency = 50
    start_cycle = 1
    step = 1

    return results_path, nu, rho, dt, T, velocity_degree, save_frequency, start_cycle, step, indices_path


if __name__ == "__main__":
    test_compute_hemodynamic_indices()
    test_compute_hemodynamic_indices_averaged_over_one_cycle()
