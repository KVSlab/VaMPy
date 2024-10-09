import os
from os import path

from dolfin import Mesh

from vampy.automatedPreprocessing.automated_preprocessing import (
    read_command_line,
    run_pre_processing,
)
from vampy.automatedPreprocessing.preprocessing_common import read_polydata


def test_mesh_model_with_one_inlet():
    model_path = "tests/test_data/ventricle/ventricle.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(
        dict(
            meshing_method="diameter",
            smoothing_method="taubin",
            refine_region=False,
            coarsening_factor=1.3,
            visualize=False,
            compress_mesh=False,
            outlet_flow_extension_length=1,
            inlet_flow_extension_length=1,
        )
    )

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

    tear_down(model_path)


def test_mesh_model_with_one_inlet_and_two_outlets():
    model_path = "tests/test_data/artery/artery.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(
        dict(
            meshing_method="curvature",
            smoothing_method="taubin",
            refine_region=False,
            coarsening_factor=1.9,
            visualize=False,
            compress_mesh=False,
            outlet_flow_extension_length=1,
            inlet_flow_extension_length=1,
        )
    )
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

    tear_down(model_path)


def test_mesh_model_with_one_inlet_and_one_outlet():
    model_path = "tests/test_data/vein/vein.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(
        dict(
            meshing_method="diameter",
            refine_region=False,
            coarsening_factor=1.3,
            visualize=False,
            compress_mesh=False,
            outlet_flow_extension_length=0.5,
            inlet_flow_extension_length=0.5,
        )
    )

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

    num_points = 6159
    num_cells = 34337

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells

    tear_down(model_path)


def test_mesh_model_with_geodesic_meshing():
    model_path = "tests/test_data/artery/artery.stl"
    # Get default input parameters
    common_input = read_command_line(model_path)
    region_point = [33.6612, 32.8443, 40.9184]
    common_input.update(
        dict(
            meshing_method="geodesic",
            max_geodesic_distance=2.5,
            edge_length=0.5,
            refine_region=True,
            region_points=region_point,
            visualize=False,
            compress_mesh=False,
            outlet_flow_extension_length=1,
            inlet_flow_extension_length=1,
        )
    )

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

    num_points = 6261
    num_cells = 34849

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells

    tear_down(model_path)

    # Remove additional file when performing refinement
    os.remove(model_path.replace(".stl", "_sac_centerline_0.vtp"))


def tear_down(model_path):
    file_extensions = [".vtu", ".xml", "_info.json", "_probe_point.json"]
    model_name = model_path.split(".")[0]
    for file_extension in file_extensions:
        file_to_remove = model_name + file_extension
        os.remove(file_to_remove)


if __name__ == "__main__":
    test_mesh_model_with_one_inlet()
    test_mesh_model_with_one_inlet_and_one_outlet()
    test_mesh_model_with_one_inlet_and_two_outlets()
    test_mesh_model_with_geodesic_meshing()
