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

    # Check that mesh is created
    mesh_path_vtu = model_path.replace("stl", "vtu")
    mesh_path_xml = model_path.replace("stl", "xml")

    assert path.isfile(mesh_path_vtu)
    assert path.isfile(mesh_path_xml)

    # Check that mesh is not empty with VTK/morphMan and FEniCS and contains correct amount of points and cells
    mesh_vtu = read_polydata(mesh_path_vtu)
    mesh_xml = Mesh(mesh_path_xml)

    num_points = 3473
    num_cells = 19208

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells


def test_mesh_model_with_one_inlet_and_two_outlets():
    model_path = "tests/test_data/artery/artery.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(dict(meshing_method="curvature",
                             smoothing_method="taubin",
                             refine_region=False,
                             coarsening_factor=1.9,
                             visualize=False,
                             compress_mesh=False,
                             outlet_flow_extension_length=1,
                             inlet_flow_extension_length=1
                             ))

    # Run pre processing
    run_pre_processing(**common_input)

    # Check that mesh is created
    mesh_path_vtu = model_path.replace("stl", "vtu")
    mesh_path_xml = model_path.replace("stl", "xml")

    assert path.isfile(mesh_path_vtu)
    assert path.isfile(mesh_path_xml)

    # Check that mesh is not empty with VTK/morphMan and FEniCS and contains correct amount of points and cells
    mesh_vtu = read_polydata(mesh_path_vtu)
    mesh_xml = Mesh(mesh_path_xml)

    num_points = 5733
    num_cells = 31457

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells


if __name__ == "__main__":
    test_mesh_model_with_one_inlet()
    test_mesh_model_with_one_inlet_and_two_outlets()
