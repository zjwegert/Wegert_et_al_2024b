# Level set-based inverse homogenisation of three-dimensional piezoelectric materials

This repository contains the algorithms and results for the following manuscript:

> Zachary J. Wegert, Anthony P. Roberts, and Vivien J. Challis (2024). "Level set-based inverse homogenisation of three-dimensional piezoelectric materials". [arXiv:2410.03148](https://arxiv.org/abs/2410.03148).

If you use or modify this package, we ask that you cite the above manuscript and others listed under [GridapTopOpt](https://github.com/zjwegert/GridapTopOpt.jl).

## Results
The results for the above paper are available in `/res`. This directory contains a subpackage for processing data and outputting figures (see `/res/scripts`). Raw data is contained in `/res/data`.

## Usage
This is a distributed code base and scripts require modification to be run in serial. Please read the documentation of [GridapDistributed](https://github.com/gridap/GridapDistributed.jl) prior to attempting to run the scripts.

### Usage on a PC
All scripts can be run on a local machine. Ensure that you adjust the number of CPUs and number of elements to make computations viable on your machine. Scripts can be called from the command line using the optional arguments specified by the scripts in `/scripts`. E.g.,
```
mpiexecjl -n 16 julia --project scripts/3d_max_Bulk_constrain_Vol_and_dh.jl -OPTIONS ...
```

### Usage on a HPC
Please carefully read the following:
- You should use the MPI implementation provided by the HPC system. See [these instructions](https://juliaparallel.org/MPI.jl/stable/configuration/#using_system_mpi) for details on how to detect the MPI backend.
- You should build PETSc using the system MPI implementation. See the following [instructions](https://github.com/gridap/GridapPETSc.jl?tab=readme-ov-file#installation) for details.
In `/jobs` we include job templates for PBS Pro along with installation scripts for PETSc and P4est. Note that job and installation scripts are HPC dependent, so please take this into account when running or installing software.

## Software versions
At the time of writing, we use Julia 1.9.4 and the package versions are set under `[compat]` in the `Project.toml` file. In addition, we used PETSc v3.19.5. We also installed P4est v2.8.5, although this is not used in this work.

If you wish to use newer versions of these packages, please remove the equality next to the version number under the `[compat]` heading of the `Project.toml` file. Note that although we expect that the implemetation will work for future versions of the project dependencies, this is not guaranteed. Please contact us if you notice issues with one of the dependencies.
