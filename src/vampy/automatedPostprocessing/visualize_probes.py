from os import path

import matplotlib.pyplot as plt
import numpy as np

from vampy.automatedPostprocessing.postprocessing_common import read_command_line

# Plotting variables
colors = ['red', 'blue', 'purple']


def visualize_probes(case_path, probe_saving_frequency, T, dt, probes_to_plot, show_figure=True, save_figure=False):
    """
    Loads probe points from completed CFD simulation and visualizes velocity (magnitude)
    and pressure at respective probes. Assuming results are in mm/s.

    Args:
        case_path (str): Path to results from simulation.
        probe_saving_frequency (int): Interval between saving probes.
        T (float): One cardiac cycle, in [ms].
        dt (float): Time step of simulation.
        show_figure (bool): Shows figure if True.
        save_figure (bool): Saves figure if True.
    """

    max_P, max_U, n_cols, n_probes, n_rows, pressures, velocity, velocity_u, velocity_v, velocity_w, n_timesteps \
        = load_probes(case_path, probe_saving_frequency)

    mean_velocity, kinetic_energy, turbulent_kinetic_energy, max_ke, max_tke = compute_mean_velocity_and_kinetic_energy(
        T, dt, n_timesteps, n_probes, velocity, velocity_u, velocity_v, velocity_w)

    kinetic_energy_spectrum, freq, max_kes = compute_energy_spectrum(n_probes, velocity_u, velocity_v, velocity_w)

    if not probes_to_plot:
        # Plot all probes in same plot

        # Plot velocity
        for k in range(n_probes):
            ax = plt.subplot(n_rows, n_cols, k + 1)
            plot_velocity_and_pressure(k, ax, max_P, max_U, mean_velocity, n_timesteps, pressures, velocity)

        if save_figure:
            save_path = path.join(case_path, "velocity_and_pressure.png")
            print("-- Saving figure of velocity magnitude and pressure at probes at {}".format(save_path))
            plt.savefig(save_path)
        if show_figure:
            print("-- Plotting velocity magnitude and pressure at probes")
            plt.show()

        # Plot kinetic energy
        for k in range(n_probes):
            ax = plt.subplot(n_rows, n_cols, k + 1)
            plot_kinetic_energy(k, ax, kinetic_energy, max_ke, max_tke, n_timesteps, turbulent_kinetic_energy)

        if save_figure:
            save_path = path.join(case_path, "kinetic_energy.png")
            print("-- Saving figure of kinetic energy and turbulent kinetic energy at probes at {}".format(save_path))
            plt.savefig(save_path)
        if show_figure:
            print("-- Plotting kinetic energy and turbulent kinetic energy at probes")
            plt.show()

        # Plot energy spectrum
        for k in range(n_probes):
            ax = plt.subplot(n_rows, n_cols, k + 1)
            plot_energy_spectrum(k, ax, freq, kinetic_energy_spectrum, max_kes)

        if save_figure:
            save_path = path.join(case_path, "energy_spectrum.png")
            print("-- Saving figure of energy spectrum at probes at {}".format(save_path))
            plt.savefig(save_path)
        if show_figure:
            print("-- Plotting energy spectrum at probes")
            plt.show()

        # Plot spectrogram
        for k in range(n_probes):
            ax = plt.subplot(n_rows, n_cols, k + 1)
            plot_spectrogram(k, ax, velocity_u)

        if save_figure:
            save_path = path.join(case_path, "spectrogram.png")
            print("-- Saving figure of spectrogram at probes at {}".format(save_path))
            plt.savefig(save_path)
        if show_figure:
            print("-- Plotting spectrogram at probes")
            plt.show()

    else:
        # Plot probes in separate plots
        for k in probes_to_plot:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
            plot_velocity_and_pressure(k, ax1, max_P, max_U, mean_velocity, n_timesteps, pressures, velocity)
            plot_kinetic_energy(k, ax2, kinetic_energy, max_ke, max_tke, n_timesteps, turbulent_kinetic_energy)
            plot_energy_spectrum(k, ax3, freq, kinetic_energy_spectrum, max_kes)
            plot_spectrogram(k, ax4, velocity_u)

            if save_figure:
                save_path = path.join(case_path, f"probes_{k}.png")
                print(f"-- Saving figure for probes at {save_path}")
                plt.savefig(save_path)
            if show_figure:
                print(f"-- Plotting quantities for probe {k}")
                plt.show()


def plot_spectrogram(k, ax, velocity):
    """
    Plots spectrogram at probe points.

    Args:
        velocity(numpy.ndarray): Velocity data.
    """
    # Filter out negative values
    ax.specgram(velocity[k], NFFT=256, Fs=2, noverlap=128, cmap="jet")

    # Set axis labels
    ax.set_ylabel("Frequency [Hz]", fontsize=12)
    ax.set_xlabel("Time [ms]", fontsize=12)

    # Set title to probe number
    ax.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)


def plot_energy_spectrum(k, ax, freq, kinetic_energy_spectrum, max_kes):
    """
    Plots energy spectrum at probe points.

    Args:
        case_path (str): Path to results from simulation.
        freq (numpy.ndarray): Frequency data.
        kinetic_energy_spectrum (numpy.ndarray): Kinetic energy spectrum data.
        max_kes (float): Maximum value in the kinetic energy spectrum.
        n_cols (int): Number of columns in the plot grid.
        n_probes (int): Number of probe points.
        n_rows (int): Number of rows in the plot grid.
        save_figure (bool): Saves figure if True.
        show_figure (bool): Shows figure if True.
    """
    # Create subplots
    ax.plot(freq, kinetic_energy_spectrum[k], color=colors[1], label="")
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Set axis limits
    ax.set_ylim(None, max_kes * 1.1)
    ax.set_xlim(min(freq), max(freq))

    # Set axis labels
    ax.set_ylabel("E(k)", fontsize=12, color=colors[0])
    ax.set_xlabel("k", fontsize=12)

    # Color axis ticks
    ax.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[0])

    # Set title to probe number
    ax.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)


def plot_kinetic_energy(k, ax, kinetic_energy, max_ke, max_tke, n_timesteps, turbulent_kinetic_energy):
    """
    Plots kinetic energy and turbulent kinetic energy at probe points.

    Args:
        kinetic_energy (numpy.ndarray): Kinetic energy data.
        max_ke (float): Maximum kinetic energy value.
        max_tke (float): Maximum turbulent kinetic energy value.
        n_timesteps (int): Number of time steps.
        turbulent_kinetic_energy (numpy.ndarray): Turbulent kinetic energy data.
    """
    # Generate plots for kinetic energy
    # Create subplots
    ax_twinx = ax.twinx()
    time_interval = np.linspace(0, n_timesteps, n_timesteps)
    ax.plot(time_interval, kinetic_energy[k], color=colors[0], label="")
    ax.plot(time_interval, turbulent_kinetic_energy[k], color=colors[1], label="")

    # Set axis limits
    ax.set_ylim(-1E-2, max_ke * 1.1)
    ax_twinx.set_ylim(-1E-2, max_tke * 1.1)
    ax.set_xlim(min(time_interval), max(time_interval))

    # Set axis labels
    ax.set_ylabel("Kinetic energy [m$^2$/s$^2$]", fontsize=12, color=colors[0])
    ax_twinx.set_ylabel("Turbulent kinetic\n energy [m$^2$/s$^2$]", fontsize=12, color=colors[1])
    ax.set_xlabel("Time [s]", fontsize=12)

    # Color axis ticks
    ax.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[0])
    ax_twinx.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[1])

    # Set title to probe number
    ax.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)


def plot_velocity_and_pressure(k, ax, max_P, max_U, mean_velocity, n_timesteps, pressures, velocity):
    """
    Plots velocity magnitude and pressure at probe points.

    Args:
        case_path (str): Path to results from simulation.
        max_P (float): Maximum pressure value.
        max_U (float): Maximum velocity value.
        mean_velocity (numpy.ndarray): Mean velocity data.
        n_timesteps (int): Number of time steps.
        n_cols (int): Number of columns in the plot grid.
        n_probes (int): Number of probe points.
        n_rows (int): Number of rows in the plot grid.
        pressures (numpy.ndarray): Pressure data.
        save_figure (bool): Saves figure if True.
        show_figure (bool): Shows figure if True.
        velocity (numpy.ndarray): Velocity data.
    """
    # Generate plots
    ax_twinx = ax.twinx()
    time_interval = np.linspace(0, n_timesteps, n_timesteps)
    ax.plot(time_interval, velocity[k], color=colors[0], label="")
    ax.plot(time_interval, mean_velocity[k], color=colors[2], linestyle='--', label="")
    ax_twinx.plot(time_interval, pressures[k], color=colors[1], label="")

    # Set axis limits
    ax.set_ylim(-1E-2, max_U * 1.5)
    ax_twinx.set_ylim(-1E-2, max_P * 1.5)
    ax.set_xlim(min(time_interval), max(time_interval))

    # Set axis labels
    ax.set_ylabel("Velocity [m/s]", fontsize=12, color=colors[0])
    ax_twinx.set_ylabel("Pressure [Pa]", fontsize=12, color=colors[1])
    ax.set_xlabel("Time [s]", fontsize=12)

    # Color axis ticks
    ax.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[0])
    ax_twinx.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[1])

    # Set title to probe number
    ax.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)


def compute_mean_velocity_and_kinetic_energy(T, dt, n_timesteps, n_probes, velocity, velocity_u, velocity_v,
                                             velocity_w):
    """
    Computes mean velocity and kinetic energy at probe points.

    Args:
        T (float): One cardiac cycle, in [ms].
        dt (float): Time step of simulation.
        n_timesteps (int): Number of time steps.
        n_probes (int): Number of probe points.
        velocity (numpy.ndarray): Velocity data.
        velocity_u (numpy.ndarray): Velocity component 'u' data.
        velocity_v (numpy.ndarray): Velocity component 'v' data.
        velocity_w (numpy.ndarray): Velocity component 'w' data.

    Returns:
        numpy.ndarray: Kinetic energy data.
        float: Maximum kinetic energy value.
        float: Maximum turbulent kinetic energy value.
        numpy.ndarray: Mean velocity data.
        numpy.ndarray: Turbulent kinetic energy data.
    """
    max_ke = 0
    max_tke = 0
    # FIXME: Revert
    saved_points_per_cycle = int(T / dt)
    #saved_points_per_cycle = 951
    n_cycles = int(n_timesteps / saved_points_per_cycle)
    mean_velocity = np.zeros((n_probes, n_timesteps))
    mean_velocity_u = np.zeros((n_probes, n_timesteps))
    mean_velocity_v = np.zeros((n_probes, n_timesteps))
    mean_velocity_w = np.zeros((n_probes, n_timesteps))
    kinetic_energy = np.zeros((n_probes, n_timesteps))
    turbulent_kinetic_energy = np.zeros((n_probes, n_timesteps))
    for k in range(n_probes):
        U = velocity[k]
        u = velocity_u[k]
        v = velocity_v[k]
        w = velocity_w[k]
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
        mean_velocity[k] = list(u_mean) * n_cycles
        mean_velocity_u[k] = list(u_mean_u) * n_cycles
        mean_velocity_v[k] = list(u_mean_v) * n_cycles
        mean_velocity_w[k] = list(u_mean_w) * n_cycles

        # Calculate total kinetic energy
        kinetic_energy[k] = 0.5 * (u ** 2 + v ** 2 + w ** 3)

        # Calculate fluctuating velocity
        fluctuating_velocity_u = u - mean_velocity_u[k]
        fluctuating_velocity_v = v - mean_velocity_v[k]
        fluctuating_velocity_w = w - mean_velocity_w[k]

        turbulent_kinetic_energy[k] = 0.5 * ((fluctuating_velocity_u ** 2) +
                                             (fluctuating_velocity_v ** 2) +
                                             (fluctuating_velocity_w ** 2))

        max_ke = np.max(kinetic_energy) if np.max(kinetic_energy) > max_ke else max_ke
        max_tke = np.max(turbulent_kinetic_energy) if np.max(turbulent_kinetic_energy) > max_tke else max_tke
    return mean_velocity, kinetic_energy, turbulent_kinetic_energy, max_ke, max_tke


def load_probes(case_path, probe_saving_frequency):
    """
    Loads probe data from a completed CFD simulation.

    Args:
        case_path (str): Path to results from simulation.
        probe_saving_frequency (int): Interval between saving probes.

    Returns:
        float: Maximum pressure value.
        float: Maximum velocity value.
        int: Number of columns in the plot grid.
        int: Number of probe points.
        int: Number of rows in the plot grid.
        numpy.ndarray: Pressure data.
        numpy.ndarray: Velocity data.
        numpy.ndarray: Velocity component 'u' data.
        numpy.ndarray: Velocity component 'v' data.
        numpy.ndarray: Velocity component 'w' data.
        int: Number of time steps.
    """
    max_U = 0
    max_P = 0
    n_probes = 0
    n_rows = 0
    n_cols = 0
    velocity = []
    velocity_u = []
    velocity_v = []
    velocity_w = []
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
            n_probes = u.shape[0]
            n_rows = 5
            n_cols = int(np.ceil(n_probes / n_rows))

        # Find velocity and pressure for scaling windows
        max_U = np.max(U) if np.max(U) > max_U else max_U
        p_id = int(0.1 * probe_saving_frequency)  # Avoid scaling from initial pressure
        max_P = np.max(p[:, p_id:]) if np.max(p[:, p_id:]) > max_P else max_P

        # Store values in lists for plotting
        for k in range(n_probes):
            if counter == 0:
                velocity.append(list(U[k]))
                velocity_u.append(list(u[k]))
                velocity_v.append(list(v[k]))
                velocity_w.append(list(w[k]))
                pressures.append(list(p[k]))
            else:
                velocity[k] += list(U[k])
                velocity_u[k] += list(u[k])
                velocity_v[k] += list(v[k])
                velocity_w[k] += list(w[k])
                pressures[k] += list(p[k])

        counter += 1
    # Get number of timesteps / store values
    velocity = np.array(velocity)
    velocity_u = np.array(velocity_u)
    velocity_v = np.array(velocity_v)
    velocity_w = np.array(velocity_w)
    pressures = np.array(pressures)
    # FIXME: Remove
    # n_stop = 2853
    # velocity = velocity[:, :n_stop]
    # velocity_u = velocity_u[:, :n_stop]
    # velocity_v = velocity_v[:, :n_stop]
    # velocity_w = velocity_w[:, :n_stop]
    # pressures = pressures[:, :n_stop]

    # Check if data is available
    if len(velocity[0]) > 0:
        n_timesteps = velocity.shape[1]
    else:
        print("-- No data to visualize. Exiting..")
        exit()

    return max_P, max_U, n_cols, n_probes, n_rows, pressures, velocity, velocity_u, velocity_v, velocity_w, n_timesteps


def compute_energy_spectrum(n_probes, velocity_u, velocity_v, velocity_w):
    """
    Computes the energy spectrum for each probe based on velocity components.

    Args:
        n_probes (int): Number of probe points.
        velocity_u (numpy.ndarray): Velocity component 'u' data for each probe.
        velocity_v (numpy.ndarray): Velocity component 'v' data for each probe.
        velocity_w (numpy.ndarray): Velocity component 'w' data for each probe.

    Returns:
        numpy.ndarray: Kinetic energy spectrum for each probe.
        numpy.ndarray: Frequency values.
        float: Maximum value in the kinetic energy spectrum.
    """
    max_kes = 0

    # Define FFT function
    def get_fft(x):
        n = len(x)
        fft = np.fft.fft(x)
        freq = np.fft.fftfreq(n)
        return fft, freq

    # Create empty list to store kinetic energy spectrum
    kinetic_energy_spectrum = []

    # Iterate over each probe
    for k in range(n_probes):
        # FFT of the velocity signals
        fft_u, freq_u = get_fft(velocity_u[k])
        fft_v, freq_v = get_fft(velocity_v[k])
        fft_w, freq_w = get_fft(velocity_w[k])

        # Compute the energy spectrum (proportional to square of absolute fft values)
        E_u = np.abs(fft_u) ** 2
        E_v = np.abs(fft_v) ** 2
        E_w = np.abs(fft_w) ** 2

        # Summing up energy spectrum from all 3 velocity components
        E_total = E_u + E_v + E_w

        # Store positive values of total kinetic energy spectrum
        positive_indices = freq_u > 0
        freq = freq_u[positive_indices]
        kinetic_energy_spectrum.append(E_total[positive_indices])
        max_kes = np.max(kinetic_energy_spectrum) if np.max(kinetic_energy_spectrum) > max_kes else max_kes

    # Convert list to numpy array
    kinetic_energy_spectrum = np.array(kinetic_energy_spectrum)

    return kinetic_energy_spectrum, freq, max_kes


def main_probe():
    # Read command line arguments
    folder, _, _, dt, _, _, probe_freq, T, _, _, _, _, _, probes_to_plot = read_command_line()
    visualize_probes(folder, probe_freq, T, dt, probes_to_plot, save_figure=True)


if __name__ == '__main__':
    main_probe()
