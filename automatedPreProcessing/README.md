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
    "remoteDataDir": "TEST_DIR",

    /* All files in this directory are copied to remoteDataDir */
    "localDataDir": "example",

    /* After the files are copied this script is run. On abel
      this script should create a job script and add to the queue */
    "script": "example/example_script.sh"
}
```

To run simulation on a remote ssh server start the script with ´--simulationConfig path/to/config.json´.

```bash
$ python automatedPreProcessing.py --simulationConfig example/ssh_config.json
```
