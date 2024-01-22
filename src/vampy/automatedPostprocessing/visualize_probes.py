from os import path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from vampy.automatedPostprocessing.postprocessing_common import read_command_line

# Plotting variables
colors = ['red', 'blue', 'purple']


def visualize_probes(case_path, probe_saving_frequency, T, dt, probes_to_plot, show_figure=False, save_figure=False):
    """
    Loads probe points from completed CFD simulation and visualizes (1) velocity (magnitude), mean velocity,
    and pressure, (2) kinetic energy, turbulent kinetic energy, and turbulence intensity, (3) kinetic energy spectrum,
    and (4) spectrogram of velocity (u-component) at respective probes. Assuming results are in m/s.

    Args:
        case_path: Path to results from simulation.
        probe_saving_frequency: Interval between saving probes.
        T: One cardiac cycle, in [ms].
        dt: Time step of simulation.
        probes_to_plot: List of integers corresponding to single probe points
        show_figure: Shows figure if True.
        save_figure: Saves figure if True.
    """
    max_P, max_U, n_cols, n_probes, n_rows, pressures, velocity, velocity_u, velocity_v, velocity_w, n_timesteps \
        = load_probes(case_path, probe_saving_frequency)

    mean_velocity, kinetic_energy, turbulent_kinetic_energy, turbulence_intensity, max_ke, max_tke \
        = compute_mean_velocity_and_kinetic_energy(T, dt, n_timesteps, n_probes, velocity, velocity_u, velocity_v,
                                                   velocity_w)
    kinetic_energy_spectrum, freq, max_kes = compute_energy_spectrum(n_probes, velocity_u, velocity_v, velocity_w, dt)

    if not probes_to_plot:
        # Plot all probes in same plot
        plot_all_probes(n_probes, n_rows, n_cols, case_path, max_P, max_U, mean_velocity, n_timesteps, pressures,
                        velocity, kinetic_energy, max_ke, max_tke, turbulent_kinetic_energy, turbulence_intensity, freq,
                        kinetic_energy_spectrum, max_kes, velocity_u, save_figure, show_figure)
    else:
        # Plot probes in separate plots
        for k in probes_to_plot:
            PLOT = False
            if PLOT:
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
                plot_velocity_and_pressure(k, ax1, max_P, max_U, mean_velocity, n_timesteps, pressures, velocity)
                plot_kinetic_energy(k, ax2, kinetic_energy, max_ke, max_tke, n_timesteps, turbulent_kinetic_energy,
                                    turbulence_intensity)
                plot_energy_spectrum(k, ax3, freq, kinetic_energy_spectrum, max_kes)
                plot_spectrogram(k, ax4, velocity_u)

                print(f"-- Plotting quantities for probe {k}")
                save_and_show_plot(case_path, f"probes_{k}.png", save_figure, show_figure)
            else:
                id_half = int(5000 / 2)
                peak_filling = max(mean_velocity[k][:id_half])
                peak_emptying = max(mean_velocity[k][id_half:id_half * 2])

                return peak_filling, peak_emptying


def save_and_show_plot(case_path, filename, save, show):
    """
    Saves and/or shows a figure based on the given flags.

    Args:
        case_path: Path to the results from the simulation.
        filename: Name of the file to save the plot.
        save: If True, the figure will be saved.
        show: If True, the figure will be displayed.
    """
    if save:
        save_path = path.join(case_path, filename)
        print(f"-- Saving figure at {save_path}")
        plt.savefig(save_path)
    if show:
        print("-- Showing figure")
        plt.show()


def plot_all_probes(n_probes, n_rows, n_cols, case_path, max_P, max_U, mean_velocity, n_timesteps, pressures,
                    velocity, kinetic_energy, max_ke, max_tke, turbulent_kinetic_energy, turbulence_intensity, freq,
                    kinetic_energy_spectrum, max_kes, velocity_u, save_figure, show_figure):
    """
    Plots data for all probes.

    Args:
        All the input parameters for the respective plotting functions.
    """
    # Plot velocity and pressure
    for k in range(n_probes):
        ax = plt.subplot(n_rows, n_cols, k + 1)
        plot_velocity_and_pressure(k, ax, max_P, max_U, mean_velocity, n_timesteps, pressures, velocity)
    save_and_show_plot(case_path, "velocity_and_pressure.png", save_figure, show_figure)

    # Plot kinetic energy
    for k in range(n_probes):
        ax = plt.subplot(n_rows, n_cols, k + 1)
        plot_kinetic_energy(k, ax, kinetic_energy, max_ke, max_tke, n_timesteps, turbulent_kinetic_energy,
                            turbulence_intensity)
    save_and_show_plot(case_path, "kinetic_energy.png", save_figure, show_figure)

    # Plot energy spectrum
    for k in range(n_probes):
        ax = plt.subplot(n_rows, n_cols, k + 1)
        plot_energy_spectrum(k, ax, freq, kinetic_energy_spectrum, max_kes)
    save_and_show_plot(case_path, "energy_spectrum.png", save_figure, show_figure)

    # Plot spectrogram
    for k in range(n_probes):
        ax = plt.subplot(n_rows, n_cols, k + 1)
        plot_spectrogram(k, ax, velocity_u)
    save_and_show_plot(case_path, "spectrogram.png", save_figure, show_figure)


def plot_spectrogram(k, ax, velocity, color_map='jet', font_size=12):
    """
    Plots spectrogram at probe points.

    Args:
        k (int): Current probe index.
        ax (matplotlib.axis.Axes): The matplotlib axis to plot on.
        velocity(numpy.ndarray): Velocity data.
        color_map (str): Color map for the spectrogram plot.
        font_size (int): Font size for the labels.
    """
    # Filter out negative values
    ax.specgram(velocity[k], NFFT=256, Fs=2, noverlap=128, cmap=color_map)

    # Set axis labels
    ax.set_ylabel("Frequency [Hz]", fontsize=font_size)
    ax.set_xlabel("Time [ms]", fontsize=font_size)

    # Set title to probe number
    ax.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)


def plot_energy_spectrum(k, ax, freq, kinetic_energy_spectrum, max_kes):
    """
    Plots energy spectrum at probe points.

    Args:
        k (int): Current probe index.
        ax (matplotlib.axis.Axes): The matplotlib axis to plot on.
        freq (numpy.ndarray): Frequency data.
        kinetic_energy_spectrum (numpy.ndarray): Kinetic energy spectrum data.
        max_kes (float): Maximum value in the kinetic energy spectrum.
    """

    # Create subplots
    def E(k):
        """
        Computes the energy spectrum E(k) of turbulent flow, using Kolmogorov's -5/3 law.

        Parameters:
        k (numpy.array): List of wave numbers in Hz.

        Returns:
        E_k (numpy.array): Energy spectrum E(k) at the given wave numbers.
        """
        return k ** (-5 / 3)

    energy_spectrum = E(freq)
    ax.plot(freq, energy_spectrum, color="black", linestyle="--", linewidth=3, label="")
    ax.plot(freq, kinetic_energy_spectrum[k], color=colors[1], label="")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.annotate('$k^{-5/3}$', xy=(np.mean(freq) / 4, np.mean(kinetic_energy_spectrum[k])), xycoords='data',
                textcoords='data', color='black', fontsize=15)

    # Set axis limits
    ax.set_ylim(None, max_kes * 1.1)
    ax.set_xlim(min(freq), max(freq))

    # Set axis labels
    ax.set_ylabel("E(k)", fontsize=12)
    ax.set_xlabel("k", fontsize=12)

    # Set title to probe number
    ax.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)


def plot_kinetic_energy(k, ax, kinetic_energy, max_ke, max_tke, n_timesteps, turbulent_kinetic_energy,
                        turbulence_intensity):
    """
    Plots kinetic energy and turbulent kinetic energy at probe points.

    Args:
        k (int): Current probe index.
        ax (matplotlib.axis.Axes): The matplotlib axis to plot on.
        kinetic_energy (numpy.ndarray): Kinetic energy data.
        max_ke (float): Maximum kinetic energy value.
        max_tke (float): Maximum turbulent kinetic energy value.
        n_timesteps (int): Number of time steps.
        turbulent_kinetic_energy (numpy.ndarray): Turbulent kinetic energy data.
    """
    # Generate plots for kinetic energy
    ax_twinx = ax.twinx()
    time_interval = np.linspace(0, n_timesteps, n_timesteps)
    ax.plot(time_interval, kinetic_energy[k], color=colors[0], label="")
    ax.plot(time_interval, turbulent_kinetic_energy[k], color=colors[1], label="")
    ax_twinx.plot(time_interval, turbulence_intensity[k], color=colors[2], label="")

    # Set axis limits
    ax.set_ylim(-1E-2, max_ke * 1.1)
    ax.set_xlim(min(time_interval), max(time_interval))

    # Set axis labels
    ax.set_ylabel("Energy [m$^2$/s$^2$]", fontsize=12)
    ax_twinx.set_ylabel("Turbulence intensity [-]", fontsize=12, color=colors[2])
    ax.set_xlabel("Time [s]", fontsize=12)

    # Color axis ticks
    ax.tick_params(axis='y', which='major', labelsize=12)
    ax_twinx.tick_params(axis='y', which='major', labelsize=12, labelcolor=colors[2])

    # Set title to probe number
    ax.set_title('Probe {}'.format(k + 1), y=1.0, pad=-14)


def plot_velocity_and_pressure(k, ax, max_P, max_U, mean_velocity, n_timesteps, pressures, velocity):
    """
    Plots velocity magnitude and pressure at probe points.

    Args:
        k (int): Current probe index.
        ax (matplotlib.axis.Axes): The matplotlib axis to plot on.
        max_P (float): Maximum pressure value.
        max_U (float): Maximum velocity value.
        mean_velocity (numpy.ndarray): Mean velocity data.
        n_timesteps (int): Number of time steps.
        pressures (numpy.ndarray): Pressure data.
        velocity (numpy.ndarray): Velocity data.
    """
    # Generate plots
    id_half = int(5000 / 2)
    peak_filling = max(mean_velocity[k][:id_half])
    peak_emptying = max(mean_velocity[k][id_half:id_half * 2])
    print(f"Peak emptying is {peak_emptying * 100:.02f} cm/s")
    print(f"Peak filling is {peak_filling * 100:.02f} cm/s")

    ax_twinx = ax.twinx()
    time_interval = np.linspace(0, n_timesteps, n_timesteps)
    ax.plot(time_interval, velocity[k], color=colors[0], label="")
    ax.plot(time_interval, mean_velocity[k], color=colors[2], linestyle='--', label="")
    # ax_twinx.plot(time_interval, pressures[k], color=colors[1], label="")

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
    saved_points_per_cycle = int(T / dt)
    # FIXME: Revert
    # saved_points_per_cycle = 951
    n_cycles = int(n_timesteps / saved_points_per_cycle)

    mean_velocity = np.zeros((n_probes, n_timesteps))
    kinetic_energy = np.zeros((n_probes, n_timesteps))
    turbulent_kinetic_energy = np.zeros((n_probes, n_timesteps))
    turbulence_intensity = np.zeros((n_probes, n_timesteps))
    for k in range(n_probes):
        U = velocity[k]
        u = velocity_u[k]
        v = velocity_v[k]
        w = velocity_w[k]

        # Calculate mean velocities
        u_mean = np.mean(U.reshape(-1, saved_points_per_cycle), axis=0)
        u_mean_u = np.mean(u.reshape(-1, saved_points_per_cycle), axis=0)
        u_mean_v = np.mean(v.reshape(-1, saved_points_per_cycle), axis=0)
        u_mean_w = np.mean(w.reshape(-1, saved_points_per_cycle), axis=0)

        # Replicate mean velocities to match original array shape
        mean_velocity[k] = np.tile(u_mean, n_cycles)
        u_mean_u = np.tile(u_mean_u, n_cycles)
        u_mean_v = np.tile(u_mean_v, n_cycles)
        u_mean_w = np.tile(u_mean_w, n_cycles)

        # Calculate total kinetic energy
        kinetic_energy[k] = 0.5 * (u ** 2 + v ** 2 + w ** 2)

        # Calculate fluctuating velocity
        fluctuating_velocity_u = u - u_mean_u
        fluctuating_velocity_v = v - u_mean_v
        fluctuating_velocity_w = w - u_mean_w

        turbulent_kinetic_energy[k] = 0.5 * ((fluctuating_velocity_u ** 2) +
                                             (fluctuating_velocity_v ** 2) +
                                             (fluctuating_velocity_w ** 2))
        u_prime = np.sqrt(2 / 3 * turbulent_kinetic_energy[k])
        U_bar = np.sqrt(u_mean_u ** 2 + u_mean_v ** 2 + u_mean_w ** 2)
        turbulence_intensity[k] = u_prime / U_bar

    max_ke = np.max(kinetic_energy)
    max_tke = np.max(turbulent_kinetic_energy)

    return mean_velocity, kinetic_energy, turbulent_kinetic_energy, turbulence_intensity, max_ke, max_tke


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
        # 25000 => Ignore firs5
        tstep = probe_saving_frequency * (counter + 1) + 25000

        probe_data = load_probe_data(case_path, tstep)
        if probe_data is None:
            break

        p, U, u, v, w = probe_data

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
    # # # FIXME: Remove
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


def load_probe_data(case_path, tstep):
    """
    Loads probe data for a given timestep.

    Args:
        case_path (str): Path to results from simulation.
        tstep (int): Current timestep.

    Returns:
        numpy.ndarray: Pressure data.
        numpy.ndarray: Velocity data.
    """
    # Load velocity component and pressure probes
    u_path = path.join(case_path, f"u_x_{tstep}.probes")
    v_path = path.join(case_path, f"u_y_{tstep}.probes")
    w_path = path.join(case_path, f"u_z_{tstep}.probes")
    p_path = path.join(case_path, f"p_{tstep}.probes")

    try:
        p = np.load(p_path, allow_pickle=True)
        u = np.load(u_path, allow_pickle=True)
        v = np.load(v_path, allow_pickle=True)
        w = np.load(w_path, allow_pickle=True)
    except Exception:
        print("-- Finished reading in probes")
        return None

    # Create velocity magnitude
    U = np.linalg.norm([u, v, w], axis=0)

    return p, U, u, v, w


def compute_energy_spectrum(n_probes, velocity_u, velocity_v, velocity_w, dt):
    """
    Computes the energy spectrum for each probe based on velocity components.

    Args:
        n_probes (int): Number of probe points.
        velocity_u (numpy.ndarray): Velocity component 'u' data for each probe.
        velocity_v (numpy.ndarray): Velocity component 'v' data for each probe.
        velocity_w (numpy.ndarray): Velocity component 'w' data for each probe.
        dt (float): Sampling rate at probes

    Returns:
        numpy.ndarray: Kinetic energy spectrum for each probe.
        numpy.ndarray: Frequency values.
        float: Maximum value in the kinetic energy spectrum.
    """

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
        # Compute the Fourier Transform of the velocity fields
        Fu = np.fft.fft(velocity_u[k])
        Fv = np.fft.fft(velocity_v[k])
        Fw = np.fft.fft(velocity_w[k])

        # Compute the energy density
        Eu = np.abs(Fu) ** 2
        Ev = np.abs(Fv) ** 2
        Ew = np.abs(Fw) ** 2
        E_total = 0.5 * (Eu + Ev + Ew)

        # Compute the frequencies for each time point
        N = velocity_u[k].size
        freqs = np.fft.fftfreq(N, dt)

        # Store positive values of total kinetic energy spectrum
        freq = freqs[:N // 2]
        E_total = E_total[:N // 2]

        kinetic_energy_spectrum.append(E_total)

    # Convert list to numpy array
    kinetic_energy_spectrum = np.array(kinetic_energy_spectrum)

    # Compute the maximum value
    max_kes = np.max(kinetic_energy_spectrum)

    return kinetic_energy_spectrum, freq, max_kes


def main_probe():
    # Read command line arguments
    folder, _, _, dt, _, _, probe_freq, T, _, _, _, _, _ = read_command_line()
    # Load probes

    condition = "af"
    data = {
        'case_id': [],
        'velocity_filling': [],
        'velocity_emptying': []
    }
    case_to_id = pd.read_csv(f"probe_ids_{condition}.csv")
    cases = ["0001", "0003", "0004", "0005", "0006", "0007", "0008", "0009", "0019", "0020", "0021", "0022", "0023",
             "0024", "0025", "0026", "0027", "0028", "0029", "0030", "0031", "0032", "0033", "0034", "0035", "0074",
             "0075", "0076", "0077", "0078", "0079", "0080", "0081", "1029", "1030", "1031", "1032", "1033", "1034",
             "1035", "1036", "1037", "1038", "1039", "2022"]
    folders = []
    for folder, case in zip(folders, cases):
        probes_to_plot = case_to_id[case]
        try:
            filling, emptying = visualize_probes(folder, probe_freq, T, dt, probes_to_plot, show_figure=True,
                                                 save_figure=True)
            data['case_id'] = case
            data['velocity_filling'] = filling
            data['velocity_emptying'] = emptying
        except:
            print(f"-- Failed for case {case}")

    df = pd.DataFrame(data)
    df.to_csv(f"laa_velocities_{condition}.csv")


if __name__ == '__main__':
    main_probe()
