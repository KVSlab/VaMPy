import os
from os import path

from dolfin import Mesh

from vampy.automatedPreprocessing.automated_preprocessing import read_command_line, \
    run_pre_processing
from vampy.automatedPreprocessing.preprocessing_common import read_polydata


def test_mesh_model_with_one_inlet():
    model_path = "tests/test_data/tube/tube.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(dict(meshing_method="diameter",
                             smoothing_method="taubin",
                             refine_region=False,
                             coarsening_factor=1.3,
                             visualize=False,
                             compress_mesh=False,
                             outlet_flow_extension_length=1,
                             inlet_flow_extension_length=1
                             ))

    # Run pre processing
    run_pre_processing(**common_input)

    model_name = model_path.replace(".stl", "")

    # Check that output files exist
    output_files = [model_name + suffix for suffix in ['.vtu', '.xml', '_info.json', '_probe_point.json']]
    check_files(output_files)

    num_points = 3473
    num_cells = 19208
    check_mesh(output_files, num_points=num_points, num_cells=num_cells)

    # Remove output files
    for output_file in output_files:
        os.remove(output_file)


def test_mesh_model_with_one_inlet_and_two_outlets():
    model_path = "tests/test_data/artery/artery.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(dict(meshing_method="curvature",
                             smoothing_method="laplace",
                             refine_region=True,
                             region_points=[33, 30, 40],
                             coarsening_factor=1.3,
                             visualize=False,
                             compress_mesh=True,
                             outlet_flow_extension_length=1,
                             inlet_flow_extension_length=1
                             ))

    # Run pre processing
    run_pre_processing(**common_input)

    model_name = model_path.replace(".stl", "")

    # Define output files and check if they exist
    output_files = [model_name + suffix for suffix in ['.vtu', '.xml.gz', '_info.json', '_probe_point.json']]
    check_files(output_files)

    # Check that mesh is not empty with VTK/morphMan and FEniCS and contains correct amount of points and cells
    num_points = 13247
    num_cells = 74861
    check_mesh(output_files, num_points=num_points, num_cells=num_cells)

    # Remove output files
    for output_file in output_files:
        os.remove(output_file)


def test_mesh_model_with_constant_edge_length():
    model_path = "tests/test_data/artery/artery.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(dict(meshing_method="constant",
                             edge_length=0.35,
                             smoothing_method="laplace",
                             refine_region=False,
                             visualize=False,
                             compress_mesh=True,
                             outlet_flow_extension_length=1,
                             inlet_flow_extension_length=1
                             ))

    # Run pre processing
    run_pre_processing(**common_input)

    model_name = model_path.replace(".stl", "")

    # Check that output files exist
    output_files = [model_name + suffix for suffix in ['.vtu', '.xml.gz', '_info.json', '_probe_point.json']]
    check_files(output_files)

    num_points = 6964
    num_cells = 38531
    check_mesh(output_files, num_points=num_points, num_cells=num_cells)

    # Remove output files
    for output_file in output_files:
        os.remove(output_file)


def check_files(output_files):
    # Assert that output files from meshing exist
    for output_file in output_files:
        assert path.exists(output_file)
        assert path.getsize(output_file) > 0


def check_mesh(output_files, num_points, num_cells):
    # Check that mesh is not empty with VTK/morphMan and FEniCS and contains correct amount of points and cells
    mesh_path_vtu = output_files[0]
    mesh_path_xml = output_files[1]
    mesh_vtu = read_polydata(mesh_path_vtu)
    mesh_xml = Mesh(mesh_path_xml)

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells


if __name__ == "__main__":
    test_mesh_model_with_one_inlet()
    test_mesh_model_with_one_inlet_and_two_outlets()
    test_mesh_model_with_constant_edge_length()
