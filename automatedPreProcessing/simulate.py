from __future__ import print_function

import errno
import json
import os

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
        remote_folder = config['remoteFolder']
    except KeyError:
        raise ValueError('Invalid configuration file')

    # Get path to home folder on remote
    client.connect(hostname, username=username, password=password)
    stdin, stdout, stderr = client.exec_command('echo $HOME')
    home = str(stdout.read().strip())

    sftp = client.open_sftp()
    remote_folder = os.path.join(home, remote_folder)

    # Uplad run script
    sftp.chdir(os.path.join(remote_folder))
    if not exists(sftp, case_name + ".sh"):
        sftp.put(os.path.join(local_dir, case_name + ".sh"), case_name + ".sh")

    # Upload mesh 
    sftp.chdir(os.path.join(remote_folder, "mesh"))
    if not exists(sftp, case_name + ".xml.gz"):
        sftp.put(os.path.join(local_dir, case_name + ".xml.gz"), case_name + ".xml.gz")

    # Upload probe points and info
    sftp.chdir(os.path.join(remote_folder, "input"))
    if not exists(sftp, case_name + "_probe_point"):
        sftp.put(os.path.join(local_dir, case_name + "_probe_point"), case_name + "_probe_point")
    if not exists(sftp, case_name + ".txt"):
        sftp.put(os.path.join(local_dir, case_name + ".txt"), case_name + ".txt")

    if not exists(sftp, os.path.join(remote_folder, "results", case_name)):
        sftp.mkdir(os.path.join(remote_folder, "results", case_name))

    sftp.close()

    # Run script
    script_path = os.path.join(remote_folder, case_name + ".sh")
    stdin, stdout, stderr = client.exec_command(os.path.join(remote_folder, 'run.sh {}'.format(script_path)))

    for msg in stdout:
        print(msg)

    for msg in stderr:
        print(msg)

    client.close()
