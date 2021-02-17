from __future__ import print_function

import errno
import json
import os

import paramiko


def exists(sftp, path):
    """os.path.exists for paramiko's SCP object
    """
    try:
        sftp.stat(path)
    except IOError as e:
        if e.errno == errno.ENOENT:
            return False
    else:
        return True


def run_simulation(config_path, localDir, caseName):
    client = paramiko.SSHClient()
    client.load_system_host_keys()

    config = json.load(open(config_path))

    try:
        hostname = config['hostname']
        username = config['username']
        password = config['password']
        remoteFolder = config['remoteFolder']
        script = config['script']
    except KeyError:
        raise ValueError('Invalid configuration file')

    # Get path to home folder on remote
    client.connect(hostname, username=username, password=password)
    stdin, stdout, stderr = client.exec_command('echo $HOME')
    home = str(stdout.read().strip())

    sftp = client.open_sftp()
    remoteFolder = os.path.join(home, remoteFolder)

    # Uplad run script
    sftp.chdir(os.path.join(remoteFolder))
    if not exists(sftp, caseName + ".sh"):
        sftp.put(os.path.join(localDir, caseName + ".sh"), caseName + ".sh")

    # Upload mesh 
    sftp.chdir(os.path.join(remoteFolder, "mesh"))
    if not exists(sftp, caseName + ".xml.gz"):
        sftp.put(os.path.join(localDir, caseName + ".xml.gz"), caseName + ".xml.gz")

    # Upload probe points and info
    sftp.chdir(os.path.join(remoteFolder, "input"))
    if not exists(sftp, caseName + "_probe_point"):
        sftp.put(os.path.join(localDir, caseName + "_probe_point"), caseName + "_probe_point")
    if not exists(sftp, caseName + ".txt"):
        sftp.put(os.path.join(localDir, caseName + ".txt"), caseName + ".txt")

    if not exists(sftp, os.path.join(remoteFolder, "results", caseName)):
        sftp.mkdir(os.path.join(remoteFolder, "results", caseName))

    sftp.close()

    # Run script
    script_path = os.path.join(remoteFolder, caseName + ".sh")

    stdin, stdout, stderr = client.exec_command(os.path.join(remoteFolder, 'run.sh {}'.format(script_path)))

    for msg in stdout:
        print(msg)

    for msg in stderr:
        print(msg)

    client.close()
