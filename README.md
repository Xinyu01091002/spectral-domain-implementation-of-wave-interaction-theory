# MF12 Wave Reconstruction: MATLAB, Python, and C++

This repository implements and verifies nonlinear multidirectional wave reconstruction based on the third-order Madsen-Fuhrman 2012 (MF12) framework.

It now contains:

- a MATLAB reference implementation
- a Python spectral port
- a C++ spectral port
- a shared cross-language verification and benchmarking workflow

MATLAB remains the numerical ground truth. Python and C++ are checked against MATLAB on shared deterministic cases.

## Acknowledgement

This repository builds directly on the MATLAB implementation released by David R. Fuhrman through DTU Data:

- David R. Fuhrman, "Matlab implementation of the third-order, multi-directional, irregular wave theory of Madsen & Fuhrman (2012, J. Fluid Mech. 698, 304-334)", DTU Data, posted 2024-02-08:
  https://data.dtu.dk/articles/code/Matlab_implementation_of_the_third-order_multi-directional_irregular_wave_theory_of_Madsen_Fuhrman_2012_J_Fluid_Mech_698_304-334_/22060124

The coefficient computation in this repository follows that MATLAB implementation, and the retained MATLAB workflows here are used as the ground truth for verification of the Python and C++ ports.

## What This Repo Does

The core MF12 workflows are:

- direct reconstruction: build MF12 coefficients and evaluate the field directly in physical space
- spectral reconstruction: build spectral coefficients and reconstruct the field efficiently on a regular grid

The main outputs are:

- free-surface elevation `eta`
- free-surface velocity potential `phi`

In practice, this repo is used for:

- research verification of MF12 direct vs spectral implementations
- focused wave-group and crossing-sea studies
- benchmark comparisons of accuracy and speed
- cross-language port validation against the MATLAB reference

## Repository Layout

```text
.
|- matlab/                     MATLAB reference implementation and research scripts
|- python/                     Python spectral port, tests, and plotting helpers
|- cpp/                        C++ spectral port and build files
|- cross_language_comparison/  shared cases, comparison metadata, and workflow notes
|- docs/                       supporting notes and archived benchmark logs
|- outputs/                    generated results and figures
|- TODO.md
|- LICENSE
`- README.md
```

More specifically:

- [matlab](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab) contains the MF12 reference implementation, examples, tests, reproduction scripts, and verification plots.
- [python](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/python) contains the Python spectral-only package, CLI, tests, and plotting/benchmark helpers.
- [cpp](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cpp) contains the C++ spectral port and its CMake-based build.
- [cross_language_comparison](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cross_language_comparison) contains shared cases and the MATLAB/Python/C++ comparison workflow.

## Current Status

- MATLAB direct and spectral workflows are the reference implementation.
- Python matches MATLAB on the shared spectral cases to machine precision.
- C++ matches MATLAB on the shared spectral cases to machine precision.
- the MATLAB coefficient and reconstruction workflow in this repository follows and is verified against the DTU MATLAB release by David R. Fuhrman
- Shared regression cases currently include:
  - `minimal_small`
  - `wavegroup_regression`
  - `benchmark_medium`
  - `benchmark_dense_300`
  - `benchmark_dense_600`

## Quick Start

### MATLAB

From the repository root in MATLAB:

```matlab
run('matlab/setup_paths.m');
run('matlab/examples/run_direct_minimal.m');
run('matlab/examples/run_spectral_minimal.m');
run('matlab/tests/smoke_test_minimal.m');
```

Useful MATLAB entry points:

- [run_direct_minimal.m](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/examples/run_direct_minimal.m)
- [run_spectral_minimal.m](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/examples/run_spectral_minimal.m)
- [smoke_test_minimal.m](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/tests/smoke_test_minimal.m)
- [plot_phi3_wavegroup_lines.m](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/verification/plot_phi3_wavegroup_lines.m)

### Python

Recommended setup from the repository root:

```powershell
python -m pip install -e python
```

Or without installation:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
```

Useful commands:

```powershell
python -m unittest discover python/tests
python -m mf12_python.cli verify cross_language_comparison/cases/minimal_small --repeats 1
python python/scripts/run_benchmarks.py
python python/scripts/compare_python_cpp_matlab.py cross_language_comparison/cases/minimal_small cross_language_comparison/cases/wavegroup_regression cross_language_comparison/cases/benchmark_medium
```

### C++

Build from the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build cpp/build
```

Useful commands:

```powershell
cpp/build/mf12_cpp.exe inspect cross_language_comparison/cases/minimal_small
cpp/build/mf12_cpp.exe validate cross_language_comparison/cases/minimal_small
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
```

## Core MATLAB APIs

The preferred MF12 function entry points are:

- `mf12_direct_coefficients`
- `mf12_direct_surface`
- `mf12_spectral_coefficients`
- `mf12_spectral_surface`

Typical usage:

```matlab
coeffs_direct = mf12_direct_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta_direct, phi_direct] = mf12_direct_surface(order, coeffs_direct, X, Y, t);

coeffs_spec = mf12_spectral_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta_spec, phi_spec] = mf12_spectral_surface(coeffs_spec, Lx, Ly, Nx, Ny, t);
```

## Cross-Language Workflow

The shared workflow uses MATLAB as ground truth and compares Python and C++ against the same archived reference outputs.

### 1. Export or refresh shared MATLAB cases

In MATLAB:

```matlab
run('matlab/repro/export_cross_language_cases.m');
```

This writes shared deterministic cases under [cross_language_comparison/cases](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cross_language_comparison/cases).

### 2. Verify Python or C++ against MATLAB

Python:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python -m mf12_python.cli verify cross_language_comparison/cases/wavegroup_regression --repeats 1
```

C++:

```powershell
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/wavegroup_regression outputs/cross_language_comparison/verify_cpp/wavegroup_regression
```

### 3. Generate comparison summaries and figures

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python python/scripts/compare_python_cpp_matlab.py cross_language_comparison/cases/minimal_small cross_language_comparison/cases/wavegroup_regression cross_language_comparison/cases/benchmark_medium
python python/scripts/plot_component_lines_python_cpp_matlab.py
python python/scripts/plot_third_order_python_cpp_matlab.py
```

Typical outputs are written under [outputs/cross_language_comparison](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/outputs/cross_language_comparison).

## Shared Case Format

Each shared case folder contains:

- `case.json`: scalar inputs, metadata, file paths, and tolerances
- `inputs/*.csv`: input arrays such as `a`, `b`, `kx`, and `ky`
- `reference/matlab/*.csv`: MATLAB reference `eta`, `phi`, `x`, and `y`
- `reference/matlab/result.json`: MATLAB runtime metadata

For the wave-group line-comparison figures, MATLAB component references are also archived under:

- `reference/matlab_components/`

## Recommended Starting Points

- new MATLAB user:
  - [run_spectral_minimal.m](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/examples/run_spectral_minimal.m)
- direct-vs-spectral reference check:
  - [smoke_test_minimal.m](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/tests/smoke_test_minimal.m)
- wave-group verification figure:
  - [plot_phi3_wavegroup_lines.m](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/verification/plot_phi3_wavegroup_lines.m)
- Python validation:
  - [cli.py](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/python/src/mf12_python/cli.py)
- C++ validation:
  - [main.cpp](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cpp/src/main.cpp)

## Outputs and Generated Files

- Generated figures, verification results, and benchmark summaries go under [outputs](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/outputs).
- C++ build artifacts now live under [cpp/build](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cpp/build).
- `outputs/` and `cpp/build/` are generated folders and should not be versioned.

## Language-Specific Notes

See these for more focused usage details:

- [matlab/README.md](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/README.md)
- [python/README.md](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/python/README.md)
- [cpp/README.md](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cpp/README.md)
- [cross_language_comparison/README.md](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cross_language_comparison/README.md)

## Citation and License

- License: MIT, see [LICENSE](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/LICENSE)
- Citation metadata: [CITATION.cff](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/CITATION.cff)
