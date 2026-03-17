# Cross-Language Workflow

This directory documents the shared MATLAB, Python, and C++ workflow for the spectral MF12 path.

MATLAB remains the numerical ground truth. Python and C++ consume the same archived shared cases and compare their spectral outputs against the MATLAB references.

## Current Status

- `minimal_small` matches MATLAB at machine precision
- `wavegroup_regression` matches MATLAB at machine precision for both `eta` and `phi`
- `benchmark_medium` matches MATLAB at machine precision for both `eta` and `phi`
- `benchmark_dense_300` and `benchmark_dense_600` are available as larger retained-component benchmark cases

## Folder Structure

- `cross_language_comparison/`: shared portable cases plus this workflow note
- `cross_language_comparison/cases/`: language-neutral validation and benchmark cases
- `matlab/repro/`: MATLAB exporters that materialize shared cases and MATLAB-only reference artifacts
- `python/`: Python package, tests, and helper scripts for the spectral MF12 path
- `cpp/`: C++ spectral CLI that reads the same shared case format

## Standard Workflow

1. Generate or refresh shared cases from MATLAB:

```matlab
run('matlab/repro/export_cross_language_cases.m');
```

2. Verify one case from Python:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python -m mf12_python.cli verify cross_language_comparison/cases/minimal_small --repeats 1
```

3. Build the C++ CLI:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build cpp/build
```

4. Verify one case from C++:

```powershell
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
```

5. Run Python benchmarks across the shared cases:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python python/scripts/run_benchmarks.py
```

6. Generate a MATLAB-vs-Python-vs-C++ summary across the shared cases:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python python/scripts/compare_python_cpp_matlab.py cross_language_comparison/cases/minimal_small cross_language_comparison/cases/wavegroup_regression cross_language_comparison/cases/benchmark_medium
```

7. Run a retained-component and grid-size speed sweep when needed:

```powershell
python cross_language_comparison/run_speed_sweep.py --warmup
```

This script:

- derives temporary benchmark cases from `benchmark_medium`
- truncates retained components by descending `sqrt(a^2 + b^2)`
- sweeps square reconstruction grids with `Nx = Ny`
- times every available implementation among MATLAB, Python, and C++
- writes CSV and JSON summaries plus figures under `outputs/cross_language_comparison/speed_sweep/`

## Case And Result Format

Each shared case folder contains:

- `case.json`: scalar inputs, file references, metadata, and tolerances
- `inputs/*.csv`: `a`, `b`, `kx`, and `ky`
- `reference/matlab/*.csv`: MATLAB `eta`, `phi`, `x`, and `y`
- `reference/matlab/result.json`: MATLAB runtime metadata

Python and C++ verification outputs both write the same result bundle shape:

- `eta.csv`
- `phi.csv`
- `x.csv`
- `y.csv`
- `result.json`

The comparison metrics reported today are:

- `eta_max_abs_err`, `phi_max_abs_err`
- `eta_rms_err`, `phi_rms_err`
- `eta_relative_l2_err`, `phi_relative_l2_err`
- `speedup_vs_matlab_total`
- `speedup_vs_matlab_reconstruction`

## Porting Note

One important implementation trap should stay documented for future porting work:

- `Lambda3` is easy to mistranscribe because only its first grouped contribution is multiplied by the leading `h^2 / (4 * beta)` prefactor
- the later `Fnpm`, `Fnpp`, `Fmpp`, `Gnpm`, `Gnpp`, and `Gmpp` terms are separate additive terms with their own `1 / beta` scaling
- if those later terms are accidentally pulled under the same outer prefactor, the third-order `G_*` coefficients become much too large
- that bug can make `eta` wrong while `phi` still looks correct, because the retained MF12 `phi` reconstruction depends primarily on the `F_*` and `mu_*` side, while `eta` uses the `G_*` and `Lambda3` side directly

If a future port shows `phi` matching MATLAB closely but third-order `eta` being badly wrong, check the `Lambda3` transcription first.
