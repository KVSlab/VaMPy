# automatedPreProcessing
Describes a pipeline for automated pre processing of surfaces using vmtk.


## Running remote jobs
Create a configuration file following the format below. 
```
{
    "hostname": "host.name.com",
    "username": "username",
    "password": "password",

    /* Copy data in localDataDir to this directory on remote server*/
    "remoteDataDir": "REMOTE_FOLDER",

    /* All files in this directory are copied to remoteDataDir */
    "localDataDir": "MESH_FOLDER",

    /* After the files are copied this script is run. On Saga
      this script should create a job script and add to the queue */
    "script": "example/example_script.sh"
}
```

To run a simulation on a remote ssh server start the script with ´--simulationConfig path/to/config.json´.
For example, to run the example configuration file enter the following in the terminal:

```bash
$ python automatedPreProcessing/automatedPreProcessing.py --simulationConfig automatedPreProcessing/ssh_config.json
```
