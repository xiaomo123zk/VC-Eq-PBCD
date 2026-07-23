# VC-Eq-PBCD and VC-PBCD for User Equilibrium Traffic Assignment

This repository implements the VC-Eq-PBCD algorithm and the VC-PBCD algorithm for solving user equilibrium traffic assignment problems.

## Overview

This work proposes two algorithms for solving user equilibrium traffic assignment:

- VC-Eq-PBCD: a parallel and efficient method for user equilibrium traffic assignment.
- VC-PBCD: a related baseline algorithm for comparison and evaluation.

The implementation is written in C++ and uses OpenMP for parallel computation.

## Requirements

The code is designed for a C++11 compiler with OpenMP support.

### Recommended environment

- Windows operating system
- C++11-compatible compiler
  - MinGW-w64 / g++
  - Visual Studio with OpenMP enabled
- OpenMP library installed and enabled

### Build requirements

Compile with OpenMP enabled. For example, with g++:

```bash
g++ -O2 -std=c++11 -fopenmp VC-Eq-PBCD.cpp -o VC-Eq-PBCD.exe
g++ -O2 -std=c++11 -fopenmp VC-PBCD.cpp -o VC-PBCD.exe
```

If you are using Visual Studio, enable OpenMP in the project settings (for example, use the `/openmp` compiler option).

## Project structure

- code/
  - VC-Eq-PBCD.cpp: main implementation of the VC-Eq-PBCD algorithm
  - VC-PBCD.cpp: main implementation of the VC-PBCD algorithm
  - Eq-PBCD/: comparison algorithm code
  - iGP/: comparison algorithm code
  - iGreedy/: comparison algorithm code
  - iTAPAS/: comparison algorithm code
- data/
  - Network data files for different test cases, including Birmingham, Anaheim, Sioux Falls, Chicago Sketch, and Philadelphia
- paper/
  - Related paper materials
- result/
  - Result figures and generated output files

## Data and network examples

The code is configured to use the network data under the data folder. The source code contains hard-coded paths such as:

- F:\VC-Eq-PBCD\data\Birm
- F:\VC-Eq-PBCD\data\Anaheim
- F:\VC-Eq-PBCD\data\Phi

If you run this project from a different directory, make sure the file paths in the source code are updated accordingly.

## Birmingham network example

The repository includes a Birmingham network example. The source code provides a sample command-line invocation for the VC-Eq-PBCD solver:

```bash
VC-Eq-PBCD.exe 3 24 1600 3 8 2 1 1 1 1 1 1 1
```

Meaning of the arguments:

1. Network selection (`3` for Birmingham)
2. Number of OpenMP threads
3. Parallel parameter / control parameter
4. Demand or scenario selection parameter
5. Additional simulation or configuration parameter
6-13. Additional algorithm parameters

For the VC-PBCD implementation, a simpler example is:

```bash
VC-PBCD.exe 3 24
```

## Notes on execution

- The programs generate output logs and result values to the console.
- The code is intended for research and algorithm comparison purposes.
- The repository also includes comparison algorithms and related experimental results for reference.

## Related algorithms included

This repository also provides implementations or references for the following algorithms:

- Eq-PBCD
- iGP
- iGreedy
- iTAPAS

These can be used for benchmarking and comparison with the proposed VC-Eq-PBCD and VC-PBCD methods.

## Citation and reference

Please refer to the accompanying paper materials in the paper/ directory for the theoretical background, algorithm details, and experimental results.
