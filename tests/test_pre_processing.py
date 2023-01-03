from os import path

from vampy.automatedPreProcessing.automatedPreProcessing import read_command_line, run_pre_processing
from vampy.automatedPreProcessing.common import read_polydata


def test_pre_processing():
    model_path = "tests/Case_test_tiny_artery/tiny_artery.stl"
    # Get default input parameters
    common_input = read_command_line()
    common_input.update(dict(meshing_method="diameter",
                             smoothing_method="taubin",
                             filename_model=model_path,
                             refine_region=False,
                             coarsening_factor=1.3,
                             viz=False,
                             compress_mesh=False,
                             outlet_flow_extension_length=1,
                             inlet_flow_extension_length=1
                             ))

    # Run pre processing
    run_pre_processing(**common_input)

    # Check that mesh is created
    mesh_path = model_path.replace("stl", "vtu")

    assert path.isfile(mesh_path)

    # Check that mesh is not empty
    mesh = read_polydata(mesh_path)

    assert mesh.GetNumberOfPoints() > 0


if __name__ == "__main__":
    test_pre_processing()
