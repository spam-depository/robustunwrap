# Robust Phase Unwrapping

## Compilation and Installation

Ensure you have CMake version 3.17 or higher and a reasonable C/C++ compiler.
Clone this repository to your computer.

### Python Setup for the Python Interface

First, create a virtual environment because for building the Python interface
and more specifically a Python package, we need `setuptools` and a working
Python interpreter. Required is Python >= 3.9 with the `venv` package installed.

Create a new virtual environment directly inside the repository root directory:

```shell
python3 -m venv venv
```

Activate the virtual environment:

```shell
source venv/bin/activate
```

Update everything and ensure you have `setuptools` and `wheel` installed.
You probably also want `numpy` but who knows:

```shell
pip install --upgrade pip
pip install setuptools wheel
```

#### Alternative: Specify an existing Python virtual environment

If you have already an existing Python environment, there is no need to create a new one.
You can specify it in the CMake call by using the addition option

```shell
-DPython3_ROOT_DIR="path/to/another/venv"
```

### Build the Python Interface from Commandline

Ensure the virtual Python environment is activated!
Create a new directory `build` and run CMake setting the flag for the Python
Interface to ON:

```shell
mkdir build
cd build
cmake -DPYTHON_INT=1 ..
make help
```

You see various targets to build. You can call `make` with any target, but
you probably want either `PyPackageBuild` which just _builds_ the Python package
for the interface or `PyPackageInstall` which additionally installs it into 
the current virtual environment.

### Using the Python Interface

After either manually or automatically installing the robustunwrap Python package,
you can load it like any other package in your Python code:

```python
import robustunwrap
```

### Building the Matlab Interface

This should work out of the box if you have a Matlab installation.
The only thing required is that you specify the root folder of Matlab if it
is not in the `PATH` during the run of Cmake

```shell
cmake -DMATLAB_INT=1 -DMatlab_ROOT_DIR="path/to/matlab/R2022a" ..
```