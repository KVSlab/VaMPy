import re
import subprocess
from os import chdir, path

import pytest

number_pattern = r"(\d+.\d+)"


def get_data_file_path(filename):
    # Construct the path to the mesh file
    current_script_path = path.abspath(__file__)
    current_directory = path.dirname(current_script_path)
    data_file_path = path.join(current_directory, "test_data", "mesh", filename)

    return data_file_path


@pytest.mark.parametrize("num_processors", [2])
def test_artery_problem(num_processors, save_cwd):
    # Simulation parameters
    mesh_path = get_data_file_path("artery.xml")
    dt = 0.951
    T = 10 * dt

    # Navigate to the simulation directory
    chdir("src/vampy/simulation")

    cmd = [
        "mpirun",
        "-np",
        f"{num_processors}",
        "oasismove",
        "NSfracStep",
        "solver=IPCS_ABCN",
        f"T={T}",
        f"dt={dt}",
        "problem=Artery",
        f"mesh_path={mesh_path}",
        "dump_probe_frequency=10",
    ]

    # Run Oasis
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Search for pressure value in print and assert
    output = result.stdout
    pattern = re.compile(r"pressure: " + number_pattern)
    pressures = []
    for match in pattern.finditer(str(output)):
        pressures.append(eval(match.group(1)))

    expected_pressure_1 = 0.38403
    expected_pressure_2 = 0.61597

    tol = 1e-16

    assert abs(pressures[0] - expected_pressure_1) < tol
    assert abs(pressures[1] - expected_pressure_2) < tol

    # Search for flow rate update and assert
    pattern = re.compile(r"computed flow rate: " + number_pattern)
    flow_rates = []
    for match in pattern.finditer(str(output)):
        flow_rates.append(eval(match.group(1)))

    expected_flow_rate_1 = 0.2801
    expected_flow_rate_2 = 0.7219

    # Lower tolerance in parallel due to round off errors
    tol = 1e-3

    assert abs(flow_rates[0] - expected_flow_rate_1) < tol
    assert abs(flow_rates[1] - expected_flow_rate_2) < tol


def test_atrium_problem(save_cwd):
    # Simulation parameters
    mesh_path = get_data_file_path("atrium.xml")
    dt = 0.00951
    T = 10 * dt

    # Navigate to the simulation directory
    chdir("src/vampy/simulation")

    cmd = [
        "oasismove",
        "NSfracStep",
        "solver=IPCS_ABCN",
        f"T={T}",
        f"dt={dt}",
        "problem=Atrium",
        f"mesh_path={mesh_path}",
    ]

    # Run Oasis
    result = subprocess.run(cmd, capture_output=True, text=True)
    # Assert successful simulation
    assert result.returncode == 0

    # Search for velocity value in print and assert
    output = result.stdout
    pattern = re.compile(r"velocity=" + number_pattern)
    velocities = []
    for match in pattern.finditer(str(output)):
        velocities.append(eval(match.group(1)))

    expected_max_velocity = 0.378
    expected_mean_velocity = 0.171

    tol = 1e-16

    assert abs(velocities[0] - expected_max_velocity) < tol
    assert abs(velocities[1] - expected_mean_velocity) < tol


def test_moving_atrium_problem(save_cwd):
    # Simulation parameters
    mesh_path = get_data_file_path("atrium.xml")
    dt = 0.1
    T = 10 * dt

    # Navigate to the simulation directory
    chdir("src/vampy/simulation")

    cmd = [
        "oasismove",
        "NSfracStepMove",
        f"T={T}",
        f"dt={dt}",
        "problem=MovingAtrium",
        f"mesh_path={mesh_path}",
        "save_flow_metrics_frequency=10",
        "track_blood=False",
    ]

    # Run Oasis
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Search for velocity value in print and assert
    output = result.stdout
    pattern = re.compile(r"velocity=" + number_pattern)
    velocities = []
    for match in pattern.finditer(str(output)):
        velocities.append(eval(match.group(1)))

    expected_max_velocity = 0.388
    expected_mean_velocity = 0.163

    tol = 1e-16

    assert abs(velocities[0] - expected_max_velocity) < tol
    assert abs(velocities[1] - expected_mean_velocity) < tol
