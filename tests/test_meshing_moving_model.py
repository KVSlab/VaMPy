import shutil
from os import listdir, path, remove

from dolfin import Mesh
from morphman import read_polydata

from vampy.automatedPreprocessing.automated_preprocessing import run_pre_processing, read_command_line


def test_meshing_moving_model():
    model_path = "tests/test_data/moving/moving.vtp"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(dict(meshing_method="diameter",
                             coarsening_factor=1.3,
                             smoothing_method="laplace",
                             moving_mesh=True,
                             is_atrium=True,
                             clamp_boundaries=True,
                             refine_region=False,
                             visualize=False,
                             compress_mesh=False,
                             add_boundary_layer=False,
                             outlet_flow_extension_length=1,
                             inlet_flow_extension_length=1
                             ))

    # Run pre processing
    run_pre_processing(**common_input)

    model_name = model_path.replace(".vtp", "")

    # Check that output files exist (standard output files plus artery_points.np)
    output_files = [model_name + suffix for suffix in ['.vtu', '.xml', '_info.json', '_probe_point.json']]
    for output_file in output_files:
        assert path.exists(output_file)
        assert path.getsize(output_file) > 0

    # Check moving domain meshing files
    displacement_file = [model_name + "_points.np"]

    extended_folder = model_name + "_extended"
    assert path.exists(extended_folder)
    assert listdir(extended_folder) != []

    num_points = 3141
    num_cells = 16801
    mesh_path_vtu = output_files[0]
    mesh_path_xml = output_files[1]
    mesh_vtu = read_polydata(mesh_path_vtu)
    mesh_xml = Mesh(mesh_path_xml)

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells

    # Remove output files
    for output_file in output_files + displacement_file:
        remove(output_file)

    shutil.rmtree(extended_folder)
