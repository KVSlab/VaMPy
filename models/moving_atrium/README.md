# Left atrium model for moving domain simulations

The provided left atrium model (`model.vtp`) is manually created in MeshMixer by [Henrik A. Kjeldsberg](https://github.com/hkjeldsberg/) , and can be meshed
and simulated using VaMPy given that `OasisMove` is installed.

In addition, the folder `model_moved` contains several displaced surface models, where `model_moved_00.vtp`
corresponds to the input model, `model.vtp`. These surface models are used to create the displacement matrix describing
the motion of the geometry.

To mesh the model with dynamic/moving meshing, run the following command:

```
vampy-mesh -i models/moving_atrium/model.vtp -at -dm -cl -m constant -el 1.3 -bl False -fli 1 -flo 1
```

Then, to run a moving simulation, navigate to `/simulation` and run the command

```
oasismove problem=MovingAtrium mesh_path=[PATH_TO_VAMPY]/models/moving_atrium/model.xml.gz 
```


