import sys

sys.path.append("..")
sys.path.append("../automatedPreProcessing")

from automatedPreProcessing.common import *
from automatedPreProcessing.automatedPreProcessing import read_command_line, run_pre_processing


def test_pre_processing():
    model_path = "Case_test_71/Case_test_71.vtp"
    # Get default inputsc
    common_input = read_command_line()
    common_input.update(dict(meshing_method="diameter",
                             filename_model=model_path,
                             aneurysm_present=False,
                             coarsening_factor=1.3,
                             viz=False,
                             compress_mesh=False))

    # Run pre processesing 
    run_pre_processing(**common_input)

    # Check that mesh is created
    mesh_path = model_path.replace("vtp", "vtu")

    assert path.isfile(mesh_path)

    # Check that mesh is not empty
    mesh = ReadPolyData(mesh_path)

    assert mesh.GetNumberOfPoints() > 0


if __name__ == "__main__":
    test_pre_processing()
