import errno
import json
from os import path, walk

import paramiko


def exists(sftp, path):
    """
    os.path.exists for paramiko's SCP object
    """
    try:
        sftp.stat(path)
    except IOError as e:
        if e.errno == errno.ENOENT:
            return False
    else:
        return True


def run_simulation(config_path, local_dir, case_name):
    """
    Run simulation of case on a remote ssh server with
    given input configuration.

    Args:
        config_path (str): Path to configuration file
        local_dir (str): Path to case folder
        case_name (str): Case name
    """
    client = paramiko.SSHClient()
    client.load_system_host_keys()

    config = json.load(open(config_path))

    try:
        hostname = config['hostname']
        username = config['username']
        password = config['password']
        remote_folder = config['remoteDataDir']
        local_folder = config['localDataDir']
        simulation_folder = config['simDataDir']
        post_processing_folder = config['postProcDataDir']
        job_script = config['job_script']
        flow_rate = config['flow_rate']
        problem_file = config['problem_file']
    except KeyError:
        raise ValueError('Invalid configuration file')

    # Use case folder if local folder is blank
    local_dir = local_dir if local_folder == '' else local_folder

    # Get path to home folder on remote
    client.connect(hostname, username=username, password=password)
    stdin, stdout, stderr = client.exec_command('echo $HOME')
    home = str(stdout.read().strip().decode("utf-8"))

    sftp = client.open_sftp()
    remote_folder = path.join(home, remote_folder)

    # Upload run script
    sftp.chdir(remote_folder)
    sftp.put(path.join(local_dir, job_script), job_script)

    # Upload simulation scripts
    sftp.chdir(home)
    for dirpath, _, filenames in walk(path.join(simulation_folder, "probe")):
        try:
            sftp.listdir(path.join(home, "probe"))
        except IOError:
            sftp.mkdir(path.join(home, "probe"))

        for filename in filenames:
            sftp.put(path.join(dirpath, filename), path.join(home, "probe", filename))

    if not exists(sftp, "probe"):
        sftp.put(path.join(simulation_folder, "probe"), "probe")
    if not exists(sftp, "Womersley.py"):
        sftp.put(path.join(simulation_folder, "Womersley.py"), "Womersley.py")
    if not exists(sftp, "Probe.py"):
        sftp.put(path.join(simulation_folder, "Probe.py"), "Probe.py")
    if not exists(sftp, flow_rate):
        sftp.put(path.join(simulation_folder, flow_rate), flow_rate)
    if not exists(sftp, problem_file):
        sftp.put(path.join(simulation_folder, problem_file), problem_file)

    # Upload post-processing scripts
    if not exists(sftp, "compute_flow_and_simulation_metrics.py"):
        sftp.put(path.join(post_processing_folder, "compute_flow_and_simulation_metrics.py"),
                 "compute_flow_and_simulation_metrics.py")
    if not exists(sftp, "compute_hemodynamic_indices.py"):
        sftp.put(path.join(post_processing_folder, "compute_hemodynamic_indices.py"), "compute_hemodynamic_indices.py")
    if not exists(sftp, "visualize_probes.py"):
        sftp.put(path.join(post_processing_folder, "visualize_probes.py"), "visualize_probes.py")
    if not exists(sftp, "postprocessing_common.py"):
        sftp.put(path.join(post_processing_folder, "postprocessing_common.py"), "postprocessing_common.py")

    # Upload mesh, probe points and info
    sftp.chdir(path.join(remote_folder, "oasis/mesh"))
    if not exists(sftp, case_name + ".xml.gz"):
        sftp.put(path.join(local_dir, case_name + ".xml.gz"), case_name + ".xml.gz")
    if not exists(sftp, case_name + "_probe_point"):
        sftp.put(path.join(local_dir, case_name + "_probe_point"), case_name + "_probe_point")
    if not exists(sftp, case_name + "_info.json"):
        sftp.put(path.join(local_dir, case_name + "_info.json"), case_name + "_info.json")

    if not exists(sftp, path.join(remote_folder, "results_{}".format(case_name))):
        sftp.mkdir(path.join(remote_folder, "results_{}".format(case_name)))

    sftp.close()

    # Run script
    stdin, stdout, stderr = client.exec_command('sbatch {}'.format(path.join(remote_folder, job_script)))

    for msg in stdout:
        print(msg)

    for msg in stderr:
        print(msg)

    client.close()
