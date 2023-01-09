#!/usr/bin/env python
import os
import subprocess
import sys


def main_oasis():
    # Parse mesh path if present
    for i, arg in enumerate(sys.argv[1:]):
        if "mesh_path" in arg:
            mesh_path = arg.split("=")[-1]
            mesh_path = os.path.join(os.getcwd(), mesh_path)
            sys.argv[i + 1] = "mesh_path={}".format(mesh_path)

    # Run oasis
    command = 'cd src/vampy/simulation; oasis NSfracStep {}'.format(' '.join(sys.argv[1:]))
    subprocess.run(command, shell=True)


if __name__ == '__main__':
    main_oasis()
