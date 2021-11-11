import os
import sys

sys.path.append("..")
sys.path.append("../postprocessing")

from postprocessing.visualize_probes import visualize_probes, read_command_line


def test_visualize_probes():
    # Define input parameters
    folder = "Case_test_71/Probes"
    dt = 0.0951 * 5
    no_of_cycles = 0.05

    # Run probe visualization
    visualize_probes(folder, dt, no_of_cycles, show_figure=False, save_figure=True)

    # Check that figure is saved
    save_path = os.path.join(folder, "Probes.png")
    assert os.path.exists(save_path)

    # Remove figure
    os.remove(save_path)


if __name__ == "__main__":
    test_visualize_probes_2()
