# Installation
## Installing VaMPy using `conda`

We recommend installing `VaMPy` and its dependencies through `conda`.  
Start by cloning into the repository:

``` console
$ git clone https://github.com/KVSLab/VaMPy.git
$ cd VaMPy
```

Then, using the ``environment.yml`` file in the root of the repository, you can call:

``` console
$ conda env update --file environment.yml --name your_environment
```

Next, can now activate your environment by running::

``` console
$ conda activate your_environment
```

Finally, install `VaMPy` inside your `conda` environment using `pip`: 

``` console
$ python3 -m pip install .
```

Now you are all set, and can start using the Vascular Modeling Pypeline.

## Editable installations

If you want to make changes to any of the underlying packages, you should remove them from the `environment.yml` file,
and install them from source, as described for `Oasis` in the next section.

### Installing Oasis

To install an editable version of [Oasis](https://github.com/mikaem/Oasis) on your machine, run the following commands
inside your conda environment:

``` console
$ git clone https://github.com/mikaem/Oasis
$ pip install cppimport
$ pip install --editable Oasis
```

The ``--editable`` flag installs the project in editable mode meaning that any changes to the original package will be
reflected directly in your environment.

## Running VaMPy using Docker

A `Dockerfile` is supplied in the root directory of the repository, which can build a docker-image with all
dependencies installed. The Docker-image can be built with the following command:

``` console
$ docker build -t name_of_image -f docker/Dockerfile .
```

A Docker-container can then be started with the following command:

``` console
$ docker run -ti --network=host -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $(pwd):/root/shared -w /root/shared --rm --shm-size=512m name_of_image
```

To run the VaMPy GUI, you need to call:

``` console
$ xhost +local:root
```

on your system before running the scripts.

### Note on Docker

Remember to call:

``` console
xhost -local:root
```

on the host system when you are done running the Docker container.
    