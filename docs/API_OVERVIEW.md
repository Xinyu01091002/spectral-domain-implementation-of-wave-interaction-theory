# API Overview

This repository now exposes three actively used implementation surfaces:

- MATLAB direct and spectral reconstruction, treated as the numerical reference
- Python spectral reconstruction and verification tools
- C++ spectral reconstruction and verification tools

The shared cross-language workflow uses archived MATLAB reference outputs under `cross_language_comparison/cases/` to validate the Python and C++ ports on deterministic cases.

## Recommended Mental Model

- use MATLAB when you need the reference implementation, direct reconstruction, or to refresh archived shared cases
- use Python when you want a lightweight spectral port, automated verification, or comparison and plotting scripts
- use C++ when you want the compiled spectral port, validation tooling, or performance benchmarking

## MATLAB Reference APIs

MATLAB remains the primary scientific reference for the repository.

### Direct reconstruction

```matlab
[X, Y] = meshgrid(x, y);
coeffs = mf12_direct_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta, phi] = mf12_direct_surface(order, coeffs, X, Y, t);
```

Use this path for direct physical-space evaluation and reference checks against the spectral implementation.

### Spectral reconstruction

```matlab
coeffs = mf12_spectral_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta, phi] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);
```

Use this path for FFT-based reconstruction and for exporting the reference outputs used by the Python and C++ ports.

### MATLAB entry points

- `matlab/setup_paths.m`: adds the bundled MF12 source tree to the MATLAB path
- `matlab/irregularWavesMF12/Source/mf12_direct_coefficients.m`: preferred public name for direct coefficient generation
- `matlab/irregularWavesMF12/Source/mf12_direct_surface.m`: preferred public name for direct reconstruction
- `matlab/irregularWavesMF12/Source/mf12_spectral_coefficients.m`: preferred public name for spectral coefficient generation
- `matlab/irregularWavesMF12/Source/mf12_spectral_surface.m`: preferred public name for spectral reconstruction
- `matlab/examples/run_direct_minimal.m`: smallest direct example
- `matlab/examples/run_spectral_minimal.m`: smallest spectral example
- `matlab/tests/smoke_test_minimal.m`: fast regression check
- `matlab/tests/regression_wavegroup_phi3.m`: representative directional wave-group regression
- `matlab/repro/export_cross_language_cases.m`: refreshes shared archived cases for Python and C++

## Python Spectral Interface

The Python implementation is a spectral-only port packaged under `python/src/mf12_python/`.

### Installation and test entry points

From the repository root:

```powershell
python -m pip install -e python
python -m unittest discover python/tests
```

Without installation:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
```

### Python CLI

The primary Python interface is the module CLI:

```powershell
python -m mf12_python.cli <command> ...
```

Supported commands:

- `verify <case_dir> [--repeats N] [--warmup] [--output-dir DIR]`: run the Python spectral path and compare to archived MATLAB references
- `benchmark <case_dirs...> [--repeats N] [--warmup] [--output-dir DIR]`: benchmark one or more shared cases and write CSV/JSON summaries

Typical examples:

```powershell
python -m mf12_python.cli verify cross_language_comparison/cases/minimal_small --repeats 1
python -m mf12_python.cli benchmark cross_language_comparison/cases/minimal_small cross_language_comparison/cases/benchmark_medium --repeats 3 --warmup
```

### Python supporting scripts

- `python/scripts/run_benchmarks.py`: convenience benchmark runner
- `python/scripts/compare_python_cpp_matlab.py`: aggregate MATLAB/Python/C++ comparison summaries
- `python/scripts/plot_component_lines_python_cpp_matlab.py`: comparison plots for line samples
- `python/scripts/plot_third_order_python_cpp_matlab.py`: third-order comparison plots

Generated Python outputs default under `outputs/cross_language_comparison/`.

## C++ Spectral Interface

The C++ implementation is a compiled spectral-only port with a CLI under `cpp/src/main.cpp`.

### Build entry points

From the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build cpp/build
```

### C++ CLI

The primary executable is:

```powershell
cpp/build/mf12_cpp.exe <command> ...
```

Supported commands:

- `inspect <case_dir>`: load a shared case and print a summary
- `validate <case_dir>`: validate required arrays and metadata for a shared case
- `verify <case_dir> [output_dir] [order]`: run the C++ spectral path and compare to archived MATLAB references
- `benchmark <case_dir> [repeats] [warmup]`: benchmark the spectral path without requiring MATLAB references
- `dump-coeffs <case_dir> <output_dir> [order]`: dump coefficient arrays for debugging and port-comparison work

Typical examples:

```powershell
cpp/build/mf12_cpp.exe inspect cross_language_comparison/cases/minimal_small
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
cpp/build/mf12_cpp.exe benchmark cross_language_comparison/cases/benchmark_medium 3 1
```

Optional build acceleration is available through FFTW and OpenMP, configured in `cpp/CMakeLists.txt`.

## Shared Cross-Language Case Interface

The contract between MATLAB, Python, and C++ lives under `cross_language_comparison/cases/`.

Each case directory typically contains:

- `case.json`: scalar inputs, metadata, tolerances, and file locations
- `inputs/*.csv`: input arrays such as `a`, `b`, `kx`, and `ky`
- `reference/matlab/*.csv`: archived MATLAB reference outputs
- `reference/matlab/result.json`: MATLAB runtime metadata

This layout is what the Python and C++ CLIs consume for verification and benchmarking.

## Public vs Internal Interfaces

Preferred public interfaces in this repository are:

- MATLAB: `mf12_direct_*` and `mf12_spectral_*`
- Python: `python -m mf12_python.cli ...`
- C++: `cpp/build/mf12_cpp(.exe) ...`
- shared validation data: `cross_language_comparison/cases/*`

Specialized implementation files inside the bundled MATLAB source tree, the Python package internals, and the C++ library internals should generally be treated as implementation details unless you are doing porting or debugging work.
