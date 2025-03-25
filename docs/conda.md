(install:conda)=

# Installing with `conda`

Using `conda` for installing `VaMPy` allows you to easily manage dependencies and keep them isolated from your system's
existing libraries. We recommend this option for developers who want to contribute to the software, but requires
installation of all the dependencies.

## Prerequisites

- The `conda` package manager must be installed on your computer. It can be installed
  through [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

(install:linux)=

## Installation on Linux or macOS

The easiest way to install `VaMPy` and its dependencies is via `conda-forge`. Run the following command in your terminal:

```
conda create -n your_environment -c conda-forge vampy
```

Once the installation is complete, activate the newly created environment:

```
conda activate your_environment
```

### Alternative: Manual Installation from Source

If you prefer to install `VaMPy` from source, follow these steps:

### Step 1:  Clone the `VaMPy` repository

Start by downloading and navigating to the root directory of `VaMPy` with the following command in your terminal:

``` console
$ git clone https://github.com/KVSLab/VaMPy.git
$ cd VaMPy
```

### Step 2:  Create a `conda` environment

Once you have installed `conda`, create a new environment for `VaMPy` with all its dependencies using the following
command in your terminal:

``` console
$ conda env update --file environment.yml --name your_environment
```

### Step 3: Activate the `conda` environment

After the configuration of the `conda` environment is finished, activate the newly created environment by running the
following command:

``` console
$ conda activate your_environment
```

### Step 4: Install `VaMPy` inside the `conda` environment using `pip`

Finally, you can install the `VaMPy` package inside your environment using `pip`:

``` console
$ python3 -m pip install .
```

### Step 5: Verify the installation

You can verify that `VaMPy` is installed correctly by downloading the test dependencies, and running the tests using the
following commands:

``` console
$ python3 -m pip install .[test]
$ python3 -m pytest tests 
```

## Installation on Windows

We recommend Windows users to use [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install)
and follow the [Linux](install:linux) instructions, or use [Docker](install:docker).

Alternatively, Windows users may install the `FEniCS` dependency from source, by following
the [FEniCS Reference Manual](https://fenics.readthedocs.io/en/latest/installation.html). Then, download the remaining
dependencies through `conda` by removing the `fenics` dependency inside `environment.yml` and follow the steps of
the [Linux/macOS](install:linux) installation instructions.

## Editable installation of VaMPy

If you want to make changes to any of the scripts included in `VaMPy`, you can install an editable version on your
machine by supplying the `--editable` flag:

```
$ python3 -m pip install --editable .
```

The `--editable` flag installs the project in editable mode meaning that any changes to the original package will be
reflected directly in your environment.

