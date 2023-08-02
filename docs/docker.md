(install:docker)=

# Installing with Docker

Installing `VaMPy` through Docker may be a good choice if you want to work in an isolated environment that doesn't
interfere with your system's existing libraries or dependencies. Our Docker image contains `VaMPy` with its
dependencies, and can be run on any operating system that supports Docker. The image is built to support the AMD64
architecture, making it compatible with most modern CPUs.

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) must be installed on your system.

## Step 1: Pull the Docker image

You can pull the Docker image for the latest `VaMPy` package from its
official [GitHub container registry](https://github.com/KVSlab/VaMPy/pkgs/container/vampy) using the following command:

``` console
$ docker pull ghcr.io/kvslab/vampy:latest
```

## Step 2: Run the Docker container

After pulling the Docker image, you can run the container using the following command:

``` console
$ docker run -it ghcr.io/kvslab/vampy:latest
```

## Step 3: Verify the installation

You can verify that `VaMPy` is installed correctly by downloading the test dependencies, and running the tests using the
following commands:

``` console
$ python3 -m pip install .[test]
$ python3 -m pytest tests 
```

## Building our own docker image

Instead of pulling the Docker image for `VaMPy` from GitHub, you can build it yourself using the provided `Dockerfile`.
The Dockerfile can be found in the `docker` folder located in the root of the project repository.

To build the Docker image, open a terminal window and navigate to the `docker` folder of the project. Then, execute the
following command:

``` console
docker build -t vampy .
```

This command will use the instructions in the `Dockerfile` to create a Docker image with the name `vampy`. Once the
Docker image is built, you can run a container from it using the `docker run` command.
