from __future__ import print_function

from os import path

import matplotlib.pyplot as plt
import numpy as np

from vampy.automatedPostprocessing.postprocessing_common import read_command_line


def visualize_probes(case_path, probe_saving_frequency, T, dt, show_figure=True, save_figure=False):
    """
    Loads probe points from completed CFD simulation,
    and visualizes velocity (magnitude) and pressure at respective probes
    Assuming results are in mm/s.

    Args:
        probe_saving_frequency (int): Interval between saving probes.
        case_path (str): Path to results from simulation
        dt (float): Time step of simulation
        T (float): One cardiac cycle, in [ms]
        show_figure (bool): Shows figure if True
        save_figure (bool): Saves figure if True
    """
    # Plotting variables
    colors = ['red', 'blue', 'purple']
    max_U = 0
    max_P = 0
    max_ke = 0
    max_tke = 0
    num_probes = 0
    num_rows = 0
    num_cols = 0
    velocities = []
    velocities_u = []
    velocities_v = []
    velocities_w = []
    pressures = []

    # Load and plot probes
    counter = 0
    while True:
        tstep = probe_saving_frequency * (counter + 1)

        # Load velocity component and pressure probes
        u_path = path.join(case_path, "u_x_{}.probes".format(tstep))
        v_path = path.join(case_path, "u_y_{}.probes".format(tstep))
        w_path = path.join(case_path, "u_z_{}.probes".format(tstep))
        p_path = path.join(case_path, "p_{}.probes".format(tstep))

        try:
            p = np.load(p_path, allow_pickle=True)
            u = np.load(u_path, allow_pickle=True)
            v = np.load(v_path, allow_pickle=True)
            w = np.load(w_path, allow_pickle=True)
        except Exception:
            print("-- Finished reading in probes")
            break

        # Create velocity magnitude
        U = np.linalg.norm([u, v, w], axis=0)
        if counter == 0:
            # Set plot parameters
            num_probes = u.shape[0]
            num_rows = 5
            num_cols = int(np.ceil(num_probes / num_rows))

        # Find velocity and pressure for scaling windows
        max_U = np.max(U) if np.max(U) > max_U else max_U
        p_id = int(0.1 * probe_saving_frequency)  # Avoid scaling from initial pressure
        max_P = np.max(p[:, p_id:]) if np.max(p[:, p_id:]) > max_P else max_P

        # Store values in lists for plotting
        for k in range(num_probes):
            if counter == 0:
                velocities.append(list(U[k]))
                velocities_u.append(list(u[k]))
                velocities_v.append(list(v[k]))
                velocities_w.append(list(w[k]))
                pressures.append(list(p[k]))
            else:
                velocities[k] += list(U[k])
                velocities_u[k] += list(u[k])
                velocities_v[k] += list(v[k])
                velocities_w[k] += list(w[k])
                pressures[k] += list(p[k])

        counter += 1

    # Get number of timesteps / store values
    velocities = np.array(velocities)
    velocities_u = np.array(velocities_u)
    velocities_v = np.array(velocities_v)
    velocities_w = np.array(velocities_w)
    pressures = np.array(pressures)

    # FIXME: Remove
    n_stop = 2853
    velocities = velocities[:, :n_stop]
    velocities_u = velocities_u[:, :n_stop]
    velocities_v = velocities_v[:, :n_stop]
    velocities_w = velocities_w[:, :n_stop]
    pressures = pressures[:, :n_stop]
    #
    if len(velocities[0]) > 0:
        n_timesteps = velocities.shape[1]
    else:
        print("-- No data to visualize")
        exit()

    # Generate probe based flow quantities
    # FIXME: Revert
    saved_points_per_cycle = 951  # int(T / dt)
    #
    n_cycles = int(n_timesteps / saved_points_per_cycle)

    mean_velocities = np.zeros((num_probes, n_timesteps))
    mean_velocities_u = np.zeros((num_probes, n_timesteps))
    mean_velocities_v = np.zeros((num_probes, n_timesteps))
    mean_velocities_w = np.zeros((num_probes, n_timesteps))
    kinetic_energy = np.zeros((num_probes, n_timesteps))
    turbulent_kinetic_energy = np.zeros((num_probes, n_timesteps))
    for k in range(num_probes):
        U = velocities[k]
        u = velocities_u[k]
        v = velocities_v[k]
        w = velocities_w[k]
        u_mean = 0
        u_mean_u = 0
        u_mean_v = 0
        u_mean_w = 0
        for n in range(n_cycles):
            u_mean += U[n * saved_points_per_cycle:(n + 1) * saved_points_per_cycle]
            u_mean_u += u[n * saved_points_per_cycle:(n + 1) * saved_points_per_cycle]
            u_mean_v += v[n * saved_points_per_cycle:(n + 1) * saved_points_per_cycle]
            u_mean_w += w[n * saved_points_per_cycle:(n + 1) * saved_points_per_cycle]

        u_mean_u = np.array(u_mean_u)
        u_mean_v = np.array(u_mean_v)
        u_mean_w = np.array(u_mean_w)

        u_mean /= n_cycles
        u_mean_u /= n_cycles
        u_mean_v /= n_cycles
        u_mean_w /= n_cycles
        mean_velocities[k] = list(u_mean) * n_cycles
        mean_velocities_u[k] = list(u_mean_u) * n_cycles
        mean_velocities_v[k] = list(u_mean_v) * n_cycles
        mean_velocities_w[k] = list(u_mean_w) * n_cycles

        # Calculate total kinetic energy
        kinetic_energy[k] = 0.5 * (u ** 2 + v ** 2 + w ** 3)

        # Calculate fluctuating velocities
        fluctuating_velocities_u = u - mean_velocities_u[k]
        fluctuating_velocities_v = v - mean_velocities_v[k]
        fluctuating_velocities_w = w - mean_velocities_w[k]

        turbulent_kinetic_energy[k] = 0.5 * ((fluctuating_velocities_u ** 2) +
                                             (fluctuating_velocities_v ** 2) +
                                             (fluctuating_velocities_w ** 2))

        max_ke = np.max(kinetic_energy) if np.max(kinetic_energy) > max_ke else max_ke
        max_tke = np.max(turbulent_kinetic_energy) if np.max(turbulent_kinetic_energy) > max_tke else max_tke

    # Generate plots
    for k in range(num_probes):
        # Create subplots
        ax0 = plt.subplot(num_rows, num_cols, k + 1)
        ax1 = ax0.twinx()
        time_interval = np.linspace(0, n_timesteps, n_timesteps)
        ax0.plot(time_interval, velocities[k], color=colors[0], label="")
        ax0.plot(time_interval, mean_velocities[k], color=colors[2], linestyle='--', label="")
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
        save_path = path.join(case_path, "velocity_and_pressure.png")
        print("-- Saving figure of velocity magnitude and pressure at probes at {}".format(save_path))
        plt.savefig(save_path)

    if show_figure:
        print("-- Plotting velocity magnitude and pressure at probes")
        plt.show()

    # Generate plots for kinetic energy
    for k in range(num_probes):
        # Create subplots
        ax0 = plt.subplot(num_rows, num_cols, k + 1)
        ax1 = ax0.twinx()
        time_interval = np.linspace(0, n_timesteps, n_timesteps)
        ax0.plot(time_interval, kinetic_energy[k], color=colors[0], label="")
        ax0.plot(time_interval, turbulent_kinetic_energy[k], color=colors[1], label="")

        # Set axis limits (y-dir)
        ax0.set_ylim(-1E-2, max_ke * 1.1)
        ax1.set_ylim(-1E-2, max_tke * 1.1)

        # Set axis labels
        ax0.set_ylabel("Kinetic energy [m$^2$/s$^2$]", fontsize=12, color=colors[0])
        ax1.set_ylabel("Turbulent kinetic\n energy [m$^2$/s$^2$]", fontsize=12, color=colors[1])
        ax0.set_xlabel("Time [s]", fontsize=12)

        # Color axis ticks
        ax0.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[0])
        ax1.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[1])

        # Set title to probe number
        ax0.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)

    if save_figure:
        save_path = path.join(case_path, "kinetic_energy.png")
        print("-- Saving figure of kinetic energy and turbulent kinetic energy at probes at {}".format(save_path))
        plt.savefig(save_path)

    if show_figure:
        print("-- Plotting kinetic energy and turbulent kinetic energy at probes")
        plt.show()

    # Clear plotting figure
    plt.clf()


def main_probe():
    folder, _, _, dt, _, _, probe_freq, T, _, _, _, _, _ = read_command_line()
    visualize_probes(folder, probe_freq, T, dt, save_figure=True)


if __name__ == '__main__':
    folder, _, _, dt, _, _, probe_freq, T, _, _, _, _, _ = read_command_line()
    visualize_probes(folder, probe_freq, T, dt, save_figure=True)
