.. title::  Running a CFD simulation on a HPC cluster

.. _cluster:

=========================================
Running a CFD simulation on a HPC cluster
=========================================
.. highlight:: console

Here, we present an automated pipeline for performing CFD simulations followed by post-processing of the results on a HPC cluster.
For this tutorial, we have used the supercomputer `Saga <https://documentation.sigma2.no/hpc_machines/saga.html>`_ as our starting point.

The automated procedure for running simulations on a cluster is part of the preprocessing script `automatedPreProcessing.py`.
The script includes the ``--simulation-config`` flag (``-sc`` for short), which should be followed by the path to a configuration file for the remote simulation.
As a template, we have included the configuration file ``ssh_config.json``, which the user will need to edit. The template contains the following fields, and is set up for the ``Artery.py`` problem:

.. code-block:: json

    {
      "hostname": "YOUR_HOSTNAME",
      "username": "YOUR_USERNAME",
      "password": "YOUR_PASSWORD",

      "remoteDataDir": "Oasis",
      "localDataDir": "test/Case_test_artery",
      "simDataDir": "simulation",
      "postProcDataDir": "automatedPostProcessing",

      "job_script": "artery_job.sh",
      "problem_file": "Artery.py",
      "flow_rate": "ICA_values"
    }

Here, we assume that the `Oasis` framework has been cloned or copied to the home directory on the cluster.
The capitalized values are to be edited by the user, and possibly the ``remoteDataDir`` depending on the location of the remote `Oasis` folder.
The ``localDataDir`` points to the path where the mesh, mesh information, and probe points are stored, and we assume that the file defined by the ``job_script`` key is located there as well.
The job script currently located within the ``test/Case_test_artery`` folder is specific for the ``Artery.py`` problem run on Saga, but is adaptable to other simulations and clusters.
The user will also need to edit ``artery_job.sh`` with their username information, cluster project, and FEniCS module location, all highlighted in the example job script.
Then, to perform pre-processing, CFD simulation, and post-processing through a single script, run the following command::

   $ python automatedPreProcessing/automatedPreProcessing.py -m diameter -i test/Case_test_artery/artery.vtp -c 1.3 -sc automatedPreProcessing/ssh_config.json -vz False

If the script is successful, it should output::

    Submitted batch job XXXXXXX

where `XXXXXXX` is the job ID.