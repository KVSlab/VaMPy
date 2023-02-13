import re
import subprocess

import pytest

number_pattern = r"(\d+.\d+)"


@pytest.mark.parametrize("num_processors", [1, 2])
def test_run_artery_problem(num_processors):
    # Path to test mesh relative to 'simulation' folder
    mesh_path = "../../../tests/test_data/mesh/artery/mesh.xml"
    mesh_path = "tests/test_data/mesh/artery/mesh.xml"

    # Simulation parameters
    dt = 0.0951
    T = 10 * dt

    # Command to run oasis with Artery.py problem
    cmd = (
        "cd src/vampy/simulation;" +
        "mpirun -np {} oasis NSfracStep T={} problem=Artery mesh_path={}"
    )

    d = subprocess.check_output(cmd.format(num_processors, T, mesh_path), shell=True)

    # Expected pressure split
    expected_pressure_0 = 0.607958
    expected_pressure_1 = 0.392042

    # Search for pressure value in print
    pattern = re.compile(r"pressure: " + number_pattern)
    pressures = []
    for match in pattern.finditer(str(d)):
        pressures.append(eval(match.group(1)))

    # Check that computed pressure value is same as expected
    assert pressures[0] == expected_pressure_0
    assert pressures[1] == expected_pressure_1

    # Expected flow rates after 10 time steps
    expected_flow_rate_in = 0.9243
    expected_flow_rate_out = 0.9213

    # Search for flow rate (Q)
    match_in = re.search(r"Q_in = " + number_pattern, str(d))
    match_out = re.search(r"Q_out = " + number_pattern, str(d))

    flow_rate_in = eval(match_in.groups()[0])
    flow_rate_out = eval(match_out.groups()[0])

    # Check that compute flow rates are same as expected
    assert flow_rate_in == expected_flow_rate_in
    assert flow_rate_out == expected_flow_rate_out


@pytest.mark.parametrize("num_processors", [1, 2])
def test_run_atrium_problem(num_processors):
    # Path to test mesh relative to 'simulation' folder
    mesh_path = "../../../tests/test_data/mesh/atrium/mesh.xml"

    # Simulation parameters
    dt = 0.0951
    T = 10 * dt

    # Command to run oasis with Artery.py problem
    cmd = (
        "cd src/vampy/simulation;" +
        "mpirun -np {} oasis NSfracStep T={} problem=Atrium mesh_path={}"
    )

    d = subprocess.check_output(cmd.format(num_processors, T, mesh_path), shell=True)

    # Expected pressure split
    expected_max_reynolds_number = 355
    expected_mean_reynolds_number = 158

    # Search for pressure value in print
    pattern = re.compile(r"Reynolds number=" + number_pattern)
    reynolds_numbers = []
    for match in pattern.finditer(str(d)):
        reynolds_numbers.append(eval(match.group(1)))

    # Check that computed pressure value is same as expected
    tol = 1E-16  # Set tolerance due to parallel computation

    assert abs(int(reynolds_numbers[0]) - expected_max_reynolds_number) < tol
    assert abs(int(reynolds_numbers[1]) - expected_mean_reynolds_number) < tol


if __name__ == '__main__':
    test_run_artery_problem(1)
    test_run_artery_problem(2)
    test_run_atrium_problem(1)
    test_run_atrium_problem(2)
