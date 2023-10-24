from os import path, remove

from dolfin import Mesh

from vampy.automatedPreprocessing.automated_preprocessing import read_command_line, run_pre_processing
from vampy.automatedPreprocessing.preprocessing_common import read_polydata


def test_meshing_tube():
    model_path = "tests/test_data/tube/tube.vtp"
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

    model_name = model_path.replace(".vtp", "")

    # Check that output files exist
    output_files = [model_name + suffix for suffix in ['.vtu', '.xml', '_info.json', '_probe_point.json']]
    for output_file in output_files:
        assert path.exists(output_file)
        assert path.getsize(output_file) > 0

    # Check that mesh is not empty with VTK/morphMan and FEniCS and contains correct amount of points and cells
    num_points = 3524
    num_cells = 19465
    mesh_path_vtu = output_files[0]
    mesh_path_xml = output_files[1]
    mesh_vtu = read_polydata(mesh_path_vtu)
    mesh_xml = Mesh(mesh_path_xml)

    assert mesh_vtu.GetNumberOfPoints() == num_points
    assert mesh_xml.num_cells() == num_cells

    # Remove output files
    for output_file in output_files:
        remove(output_file)
