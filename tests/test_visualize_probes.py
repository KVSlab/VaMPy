from os import path, remove

from vampy.automatedPostprocessing.visualize_probes import visualize_probes


def test_visualize_probes():
    # Path to test results and params
    probe_path = "tests/test_data/results/Probes"
    probe_frequency = 100
    file_path = path.join(probe_path, "Probes.png")

    # Run post-processing
    visualize_probes(probe_path, probe_frequency, save_figure=True, show_figure=False)

    # Check that output image exists
    assert path.exists(file_path)
    assert path.getsize(file_path) > 0

    # Remove generatedoutput
    remove(file_path)


if __name__ == "__main__":
    test_visualize_probes()
