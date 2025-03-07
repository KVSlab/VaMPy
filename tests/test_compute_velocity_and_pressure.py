# Copyright (c) 2025 Simula Research Laboratory
# SPDX-License-Identifier: GPL-3.0-or-later
from os import path, remove

from vampy.automatedPostprocessing.compute_velocity_and_pressure import compute_velocity_and_pressure


def test_compute_velocity_and_pressure():
    # Path to test results and params
    results_path = "tests/test_data/results/Solutions"
    dt = 0.951
    save_frequency = 5
    velocity_degree = 1
    pressure_degree = 1
    step = 1

    # Run post-processing
    compute_velocity_and_pressure(results_path, dt, save_frequency, velocity_degree, pressure_degree, step)

    # Check that output files exist
    metric_names = ["velocity", "pressure"]

    for name in metric_names:
        xdmf_path = path.join(results_path, "{}.xdmf".format(name))
        h5_path = path.join(results_path, "{}.h5".format(name))
        assert path.exists(xdmf_path) and path.exists(h5_path)
        assert path.getsize(xdmf_path) > 0
        assert path.getsize(h5_path) > 0

        # Remove generated output
        remove(xdmf_path)
        remove(h5_path)


if __name__ == "__main__":
    test_compute_velocity_and_pressure()
