from os import path, remove

from vampy.automatedPostprocessing.visualize_probes import visualize_probes


def test_visualize_all_probes():
    # Path to test results and params
    T = 100
    dt = 0.1
    probe_path = "tests/test_data/results/Probes"
    probe_frequency = 100
    u_p_path = path.join(probe_path, "velocity_and_pressure.png")
    ke_path = path.join(probe_path, "kinetic_energy.png")
    spectrum_path = path.join(probe_path, "energy_spectrum.png")
    spectrogram_path = path.join(probe_path, "spectrogram.png")

    # Run post-processing
    visualize_probes(probe_path, probe_frequency, T, dt, probes_to_plot=[], save_figure=True, show_figure=False)

    # Check that output image exists
    for file_path in [u_p_path, ke_path, spectrum_path, spectrogram_path]:
        assert path.exists(file_path)
        assert path.getsize(file_path) > 0

        # Remove generatedoutput
        remove(file_path)


def test_visualize_two_probes():
    # Path to test results and params
    probes = [1, 2]
    T = 100
    dt = 0.1
    probe_path = "tests/test_data/results/Probes"
    probe_frequency = 100

    # Run post-processing
    visualize_probes(probe_path, probe_frequency, T, dt, probes_to_plot=probes, save_figure=True, show_figure=False)

    # Check that output image exists
    for probe in probes:
        results_path = path.join(probe_path, f"probes_{probe}.png")
        assert path.exists(results_path)
        assert path.getsize(results_path) > 0

        # Remove generatedoutput
        remove(results_path)


if __name__ == "__main__":
    test_visualize_all_probes()
    test_visualize_two_probes()
