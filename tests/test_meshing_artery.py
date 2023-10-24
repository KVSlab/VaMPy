from os import path, remove

from dolfin import Mesh
from morphman import read_polydata

from vampy.automatedPreprocessing.automated_preprocessing import run_pre_processing, read_command_line


def test_meshing_artery():
    model_path = "tests/test_data/artery/artery.vtp"
    # Get default input parameters
    common_input = read_command_line(model_path)
    common_input.update(dict(meshing_method="constant",
                             edge_length=0.50,
                             smoothing_method="laplace",
                             refine_region=True,
                             region_points=[33, 30, 40],
                             visualize=False,
                             compress_mesh=False,
                             outlet_flow_extension_length=1,
                             inlet_flow_extension_length=1
                             ))

    # Run pre processing
    run_pre_processing(**common_input)

    model_name = model_path.replace(".vtp", "")

    # Define output files and check if they exist
    output_files = [model_name + suffix for suffix in ['.vtu', '.xml', '_info.json', '_probe_point.json']]
    for output_file in output_files:
        assert path.exists(output_file)
        assert path.getsize(output_file) > 0

    # Check that mesh is not empty with VTK/morphMan and FEniCS and contains correct amount of points and cells
    num_points = 13744
    num_cells = 79700

    mesh_path_vtu = output_files[0]
    mesh_path_xml = output_files[1]
    mesh_vtu = read_polydata(mesh_path_vtu)
    mesh_xml = Mesh(mesh_path_xml)

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells

    # Remove output files
    for output_file in output_files:
        remove(output_file)
