from __future__ import print_function

from os import path

import matplotlib.pyplot as plt
import numpy as np

from vampy.automatedPostProcessing.postprocessing_common import read_command_line


def visualize_probes(case_path, probe_saving_frequency, show_figure=True, save_figure=False):
    """
    Loads probe points from completed CFD simulation,
    and visualizes velocity (magnitude) and pressure at respective probes
    Assuming results are in mm/s.

    Args:
        probe_saving_frequency (int): Interval between saving probes.
        case_path (str): Path to results from simulation
        show_figure (bool): Shows figure if True
        save_figure (bool): Saves figure if True
    """
    # Plotting variables
    colors = ["red", "blue"]
    max_U = 0
    max_P = 0
    num_probes = 0
    num_rows = 0
    num_cols = 0
    velocitites = []
    pressures = []

    # Load and plot probes
    counter = 0
    while True:
        tstep = probe_saving_frequency * (counter + 1)

        # Load velocity component and pressure probes
        u = path.join(case_path, "u_x_{}.probes".format(tstep))
        v = path.join(case_path, "u_y_{}.probes".format(tstep))
        w = path.join(case_path, "u_z_{}.probes".format(tstep))
        p = path.join(case_path, "p_{}.probes".format(tstep))

        try:
            p_probe = np.load(p, allow_pickle=True)
            u_probe = np.load(u, allow_pickle=True)
            v_probe = np.load(v, allow_pickle=True)
            w_probe = np.load(w, allow_pickle=True)
        except:
            print("-- Finished reading in probes")
            break

        # Create velocity magnitude
        U = np.sqrt(u_probe ** 2 + v_probe ** 2 + w_probe ** 2)
        if counter == 0:
            # Set plot parameters
            num_probes = U.shape[0]
            num_rows = 5
            num_cols = int(np.ceil(num_probes / num_rows))

        # Find velocity and pressure for scaling windows
        max_U = np.max(U) if np.max(U) > max_U else max_U
        p_id = int(0.1 * probe_saving_frequency)  # Avoid scaling from initial pressure
        max_P = np.max(p_probe[:, p_id:]) if np.max(p_probe[:, p_id:]) > max_P else max_P

        # Store values in lists for plotting
        for k in range(num_probes):
            if counter == 0:
                velocitites.append(list(U[k]))
                pressures.append(list(p_probe[k]))
            else:
                velocitites[k] += list(U[k])
                pressures[k] += list(p_probe[k])

        counter += 1

    # Get number of timesteps / store values
    if len(velocitites[0]) > 0:
        n_timesteps = len(velocitites[0])
    else:
        print("-- No data to visualize")
        exit()

    # Generate plots
    for k in range(num_probes):
        # Create subplots
        ax0 = plt.subplot(num_rows, num_cols, k + 1)
        ax1 = ax0.twinx()
        time_interval = np.linspace(0, n_timesteps, n_timesteps)
        ax0.plot(time_interval, velocitites[k], color=colors[0], label="")
        ax1.plot(time_interval, pressures[k], color=colors[1], label="")

        # Set axis limits (y-dir)
        ax0.set_ylim(-1E-2, max_U * 1.5)
        ax1.set_ylim(-1E-2, max_P * 1.5)

        # Set axis labels
        ax0.set_ylabel("Velocity [m/s]", fontsize=12, color=colors[0])
        ax1.set_ylabel("Pressure [Pa]", fontsize=12, color=colors[1])
        ax0.set_xlabel("Time [s]", fontsize=12)

        # Color axis ticks
        ax0.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[0])
        ax1.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[1])

        # Set title to probe number
        ax0.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)

    if save_figure:
        save_path = path.join(case_path, "Probes.png")
        print("-- Saving figure of velocity magnitude and pressure at probes at {}".format(save_path))
        plt.savefig(save_path)

    if show_figure:
        print("-- Plotting velocity magnitude and pressure at probes")
        plt.show()

    # Clear plotting figure
    plt.clf()


if __name__ == '__main__':
    folder, _, _, dt, _, _, probe_freq, _, _, _, _, _, _ = read_command_line()
    visualize_probes(folder, probe_freq, save_figure=True)
