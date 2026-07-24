# VC-Eq-PBCD and VC-PBCD

This directory contains the C++ implementations used for user-equilibrium traffic assignment experiments.

## Requirements

- Windows
- VScode
- C++11 compiler
- OpenMP support, for example MinGW-w64 `g++` with `-fopenmp`

## Main Files

- `VC-Eq-PBCD.cpp`: VC-Eq-PBCD implementation.
- `VC-PBCD.cpp`: VC-PBCD baseline implementation.
- `Eq-PBCD`, `iGP`, `iGreedy`, `iTAPAS`: comparison algorithm entries or supporting code.

The source files currently use data files under:

```text
F:\VC-Eq-PBCD\data
```

If the project is moved, update the hard-coded data prefixes in the source files.

## Build

From `F:\VC-Eq-PBCD\code`:

```powershell
g++ -o VC-Eq-PBCD -std=c++11 -g .\VC-Eq-PBCD.cpp -fopenmp
g++ -o VC-PBCD -std=c++11 -g .\VC-PBCD.cpp -fopenmp
```

## VC-Eq-PBCD Usage

Example:

```powershell
.\VC-Eq-PBCD.exe 3 24 1600 3 8 2 1 1 1 1 1 1 1
```

Arguments:

1. Network selector: `1=Anaheim`, `2=Chicago-Sketch`, `3=Birmingham`, `4=Philadelphia`.
2. Number of OpenMP threads.
3. Initial `elseparalevel`.
4. Path-count threshold `num`.
5. OD coloring file suffix `sim`; for example `8` reads `*_od_VertexColor_equ_sim08.csv`.
6. `parameter[4]`
7. `parameter[5]`
8. `parameter[6]`
9. `parameter[7]`
10. `parameter[8]`
11. `parameter[9]`
12. `parameter[10]`
13. `parameter[11]`

Anaheim example:

```powershell
.\VC-Eq-PBCD.exe 1 4 100 3 8 2 1 2 1 1 1 1 1 > Ana-VC-Eq-PBCD.log
```

## VC-PBCD Usage

Example:

```powershell
.\VC-PBCD.exe 1 8 > Ana-VC-PBCD.log
```

Arguments:

1. Network selector: `1=Anaheim`, `2=Chicago-Sketch`, `3=Birmingham`, `4=Philadelphia`.
2. Number of OpenMP threads.


## Output

Both programs print iteration progress, objective value, UE gap, runtime, and `LinkFlow[0]` to the console. Redirect output to a log file when running experiments.
