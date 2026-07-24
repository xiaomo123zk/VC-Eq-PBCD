# VC-Eq-PBCD and VC-PBCD

This directory contains the C++ implementations used for user-equilibrium traffic assignment experiments.

## Requirements

- Windows
- A C++11 compiler
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
6. `parameter[4]`: divisor applied to `elseparalevel` after `UEGap < 1e-4`.
7. `parameter[5]`: divisor applied to `elseparalevel` after `UEGap < 1e-5`.
8. `parameter[6]`: divisor applied to `elseparalevel` after `UEGap < 1e-6`.
9. `parameter[7]`: divisor applied to `elseparalevel` after `UEGap < 1e-7`.
10. `parameter[8]`: divisor applied to `elseparalevel` after `UEGap < 1e-8`.
11. `parameter[9]`: divisor applied to `elseparalevel` after `UEGap < 1e-9`.
12. `parameter[10]`: divisor applied to `elseparalevel` after `UEGap < 1e-10`.
13. `parameter[11]`: divisor applied to `elseparalevel` after `UEGap < 1e-11`.

After each outer iteration, the code first sets `elseparalevel = level / 2` once `UEGap < 1`. It then checks the thresholds above in order. If multiple thresholds are satisfied in the same iteration, the divisors are applied cumulatively. A value of `1` leaves `elseparalevel` unchanged at that threshold; values greater than `1` reduce the batch size used for restricted updates of OD pairs that still violate the OD gap.

Anaheim example:

```powershell
.\VC-Eq-PBCD.exe 1 4 100 3 8 2 1 2 1 1 1 1 1 > Ana-VC-Eq-PBCD.log
```

## VC-PBCD Usage

Example:

```powershell
.\VC-PBCD.exe 3 8
```

Arguments:

1. Network selector: `1=Anaheim`, `2=Chicago-Sketch`, `3=Birmingham`, `4=Philadelphia`.
2. Number of OpenMP threads.

The OD coloring suffix is hard-coded in `VC-PBCD.cpp` by network:

- Anaheim: `sim09`
- Chicago-Sketch: `sim07`
- Birmingham: `sim08`
- Philadelphia: `sim09`

## Output

Both programs print iteration progress, objective value, UE gap, runtime, and `LinkFlow[0]` to the console. Redirect output to a log file when running experiments.
