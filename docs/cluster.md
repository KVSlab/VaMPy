# Running a simulation on an HPC cluster

Here, we present an automated pipeline for performing CFD simulations followed by post-processing of the results on an
HPC cluster. For this tutorial, we have used the supercomputer
[Saga](https://documentation.sigma2.no/hpc_machines/saga.html) as our starting point.

The automated procedure for running simulations on a cluster is part of the preprocessing
script `automatedPreProcessing.py` called through the  `vampy-mesh` command. The script includes the `--config-path`
flag (`-cp` for short), which should be followed by the path to a configuration file for the remote simulation. As a
template, we have included the configuration file
`ssh_config.json`, which the user will need to edit. The template contains the following fields, and is set up for
the `Artery.py`
problem:

``` json
{
  "hostname": "YOUR_HOSTNAME",
  "username": "YOUR_USERNAME",
  "password": "YOUR_PASSWORD",
  
  "remote_vampy_folder": "VaMPy",
  "local_mesh_folder": "tests/Case_test_artery",
  "job_script": "artery_job.sh"
}
```

Here, we assume that there is an available *VaMPy* environment installed on the cluster, and that the *VaMPy* repository
has been cloned to the root directory of the cluster. The capitalized values are to be edited by the user, and possibly
the `remote_vampy_folder` depending on the location of the remote *VaMPy* folder. The
`local_mesh_folder` points to the path where the mesh, mesh information, and probe points are stored, and we assume that the
file defined by the
`job_script` key is located there as well. The job script currently located within the `tests/Case_test_artery` folder
is specific for the
`Artery.py` problem run on Saga, but is adaptable to other simulations and clusters. The user will also need to
edit `artery_job.sh` with their username information, cluster project, and VaMPy environment, all highlighted in
the example job script. Then, to perform pre-processing, CFD simulation, and post-processing through a single script,
run the following command:

``` console
$ vampy-mesh -m diameter -i tests/Case_test_artery/artery.vtp -c 1.3 -cp src/vampy/automatedPreProcessing/ssh_config.json -viz False
```

If the script is successful, it should output:

``` console
Submitted batch job XXXXXXX
```

where `XXXXXXX` is the job ID.
