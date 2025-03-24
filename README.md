# Canoe: Comprehensive Atmosphere N' Ocean Engine

[![build](https://github.com/chengcli/canoe/actions/workflows/main.yml/badge.svg)](https://github.com/chengcli/canoe/actions/workflows/main.yml)
[![build](https://github.com/chengcli/canoe/actions/workflows/mac.yml/badge.svg)](https://github.com/chengcli/canoe/actions/workflows/mac.yml)
[![DocStatus](https://readthedocs.org/projects/pycanoe/badge/?version=latest)](https://pycanoe.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![codecov](https://codecov.io/gh/chengcli/canoe/branch/main/graph/badge.svg?token=hKnnv79a09)](https://codecov.io/gh/chengcli/canoe)

## Install system libraries and toolchain
Canoe can be installed on either a Linux distribution or on MacOS. Open a Linux or Mac terminal,
you can clone this repo using the following command:
```
git clone https://github.com/chengcli/canoe
```
This will copy all source files into your local computer. You will need to install a few
system libraries before installing canoe. All following instructions are executed under
the `canoe/` directory, which is referred to as the `root`.

### MacOS Installation Guide
We assume that [homebrew](https://brew.sh/) is already installed on your Mac and we will
use `brew` to install required system libraries. The system libraries are listed in
`Brewfile` at the `root`. To install them all, execute
```
brew bundle
```
### Ubuntu Linux Installation Guide
On a Ubuntu linux system, use `apt` to install
```
sudo apt install $(cat packages_debian.txt)
```

### Redhat Linux Installation Guide
On a Redhat linux system, use `yum` to install
```
sudo yum -y install $(cat packages_centos.txt)
```

### Multi-core execution
If multi-core parallelization is needed, these extra pacakges should be install
- mpich parallel library

Ubuntu linux:
```
sudo apt install libmpich-dev
```
Redhat linux:
```
sudo yum install mpich-devel
source ~/.bash_profile
```
- pnetcdf output

Redhat linux does not support pnetcdf natively. So it should be downloaded and install.
The default installation directory is $HOME/opt/
```
cd external
./fetch_pnetcdf.sh
./install_pnetcdf.sh
cd ..
```

### Notes on using openmpi
Some system, especially conda, uses openmpi by default.
It is know that openmpi sometimes causes trouble in the simulation and mpich works
better with canoe. If you have to use openmpi and the run fails immediately after
execution, try to run with a single core first and then multi-core.

## Install python libraries
The minimum python version is 3.8.
All needed python libraries are collected in `requirements.txt`. We suggest using a
python [virtual environment](https://docs.python.org/3/library/venv.html) to install
these packages. If you are already using a virtual enviroment, install python packages
by
```
pip3 install -r requirements.txt
```
Otherwise, to create a python virtual environment:
```
python -m venv pyenv
```
This command will create an environment named `pyenv` in your current directory. Then, you
can use the previous command to install the python packages.

## Install pre-commit
Register your `pre-commit` hooks using
```
pre-commit install
```
The contributor's guide explains the meaning of `pre-commit`.

## Environment variables
These following environment variables are important for the system to find the
appropriate MPI:

```
export PATH=$PATH:/usr/lib64/mpich/bin
export LD_LIBRARY_PATH=/usr/lib64/mpich/lib:$LD_LIBRARY_PATH
export MPICC=/usr/lib64/mpich/bin/mpicc
export MPICXX=/usr/lib64/mpich/bin/mpicxx
```

## Compiler
gcc9 or clang is supported and tested.

## How to build and test
After you completed the installation steps, you can build the canoe library.
The easiest way is to build it in-place, meaning that the build (binary files) are
located under `root`. To do so, make a new directory named `build`
```
mkdir build
```
All build files will be generated and placed under this directory. It is completely safe
to delete the whole directory if you want another build. `cd` to build and `cmake`

```
cd build
cmake ..
```
This command tells the cmake command to look for `CMakeFiles.txt` in the parent directory,
and start configuring the compile environment. Then compile the code by
```
make -j4
```
This comman will use 4 cores to compile the code in parallel. Once complete, all executable
files will be placed in `build/bin`.

## Optional packages
- The [Reference Forward Model](http://eodg.atm.ox.ac.uk/RFM/) (RFM) is provided optionally as
a tool to generate opacity tables. The source code of this package is not publically available.
Please contact [Anu Dudhia](mailto:anu.dudhia@physics.ox.ac.uk) ar [Cheng Li](mailto:chengcli@umich.edu) to obtain access.
The build process turns off RFM by default. To turn on building RFM, use
```
cmake .. -DRFM=ON
```
- The [DIScrete Ordinate Radiative Transfer](https://doi.org/10.1016/j.jqsrt.2011.03.019) (DISORT)
is provided optionally as a plane-parallel radiative transfer solver.
The original source code was in Fortran77. Tim Downling translated it to C in 2011.
The C-source code, version 2.1.3, is hosted at [libradtran.org](http://libradtran.org/doku.php).
The build process turns off DISORT by default. To turn on building DISORT, use
```
cmake .. -DDISORT=ON
```
- The [Parallel Kernel-Independent Fast Multipole Method](https://github.com/dmalhotra/pvfmm) (PVFMM)
is provided optionally as a potential solver for N-Body simulation.
The upstream source code is hosted at [pvfmm](https://github.com/dmalhotra/pvfmm).
We made a customized fork hosted at [chengcli:pvfmm](https://github.com/chengcli/pvfmm/).
The build process turns off PVFMM by default. To turn on building PVFMM, use
```
cmake .. -DPVFMM=ON
```

## Installation and Setup Instructions on UCLA Hoffman2

Follow the steps below if you are working on Hoffman2.

### 1. Configure Your Shell Environment

Edit your `~/.bashrc` file to ensure all shell processes use the appropriate modules:

```bash
module load gcc/11.3.0 > /dev/null 2>&1
module load python/3.12.9 > /dev/null 2>&1
module load mpich/4.3 > /dev/null 2>&1
```

### 2. Request an Interactive Session

```bash
qrsh -l h_data=1G,h_rt=2:00:00,arch=intel-E5* -pe shared 8
```

### 3. Load Required Environment Modules

```bash
module load cmake/3.30.0
module load wget
module load make
module load intel/2020.4

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/lib
export PATH=$PATH:$HOME/opt/bin
export MPICC=/u/local/mpi/mpich/4.3.0/gcc/11.3.0/bin/mpicc
export MPICXX=/u/local/mpi/mpich/4.3.0/gcc/11.3.0/bin/mpicxx
export CPATH=$CPATH:$HOME/opt/include
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HOME/opt/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$HOME/opt/include
```

As a less-cumbersome alternative, I've included canoe/dependencies file that can run everything at once from the command line (provided that you've already set up your python virtual environment; see step 4 below):

```bash
source ./dependencies
```

### 4. Set Up Python Virtual Environment

```bash
cd canoe
python -m venv ~/.virtualenvs/3.12.9/canoe
source ~/.virtualenvs/3.12.9/canoe/bin/activate
python -m pip install --upgrade pip
pip3 install -r requirements.txt
```

### 5. Build External Dependencies

```bash
cd external
./fetch_cantera.sh
./fetch_eigen.sh
./fetch_pnetcdf.sh

./install_eigen.sh
./install_cantera.sh
./install_pnetcdf.sh
```

### 6. Reinstall PyTorch (CPU Version)

```bash
pip3 install --force-reinstall torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

### 7. Compile the Model

```bash
mkdir build; cd build
cmake -DTASK=bryan -DCMAKE_INSTALL_PREFIX=$HOME/canoe/INSTDIR -DEIGEN3_INCLUDE_DIR=$HOME/opt/include/include/eigen3 ..
make -j 4
```
### 📁 Setup a New Experiment

To set up a new experiment in CANOE, follow the steps below. As an example, we'll walk through how to register a new experiment named `2025-Earth-RCE`.

#### 1. Prepare Experiment Directory
First, assemble the minimum required files under:
```
canoe/examples/2025-Earth-RCE/
```

#### 2. Register the Experiment in CANOE
You will need to update three files to properly register the new experiment:

---

#### 🔧 `canoe/examples/CMakeLists.txt`

This file adds `2025-Earth-RCE` to the compilation path. Add the following snippet:

```cmake
if (${TASK} STREQUAL "earth")
  add_subdirectory(2025-Earth-RCE)
endif()
```

---

#### 🔧 `canoe/examples/2025-Earth-RCE/CMakeLists.txt`

You can copy this file from an existing example:

```bash
cp canoe/examples/2019-Li-snap/CMakeLists.txt canoe/examples/2025-Earth-RCE/CMakeLists.txt
```

Then, edit it to reflect the new experiment setup. The file typically includes:

```cmake
# 1. Compile "earth" problem
setup_problem(earth)

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${input})
endforeach()
```

---

#### 🔧 `canoe/cmake/examples/earth.cmake`

This file configures the hydrodynamics and other runtime options for the `earth` problem. Add the following configuration block:

```cmake
# configuration for earth hydrodynamics

# athena variables
set(NUMBER_GHOST_CELLS 3)
set(EQUATION_OF_STATE ideal_moist)
set(NON_BAROTROPIC_EOS 1)
set(RSOLVER lmars)

# canoe variables
set(MPI ON)
set(PNETCDF ON)
```

---

Once these changes are made, re-run `cmake` to regenerate the build system and proceed with compiling the code as usual.

### Run the model

#### 1. Request an Interactive Session and Load Dependencies

```bash
qrsh -l h_data=1G,h_rt=2:00:00,arch=intel-E5* -pe shared 8
source ./dependencies
```

#### 2. Run the executable

```bash
mpirun -n 8 ./earth.release -i earth.inp -d /u/scratch/f/fspauldi/canoe/earth > log.earth &
```

#### 3. Process the output
If the run was successful, a set of .nc files are deposited in the /run directory. To combine them in a single .nc file, 

```bash
python3 combine.py
```
