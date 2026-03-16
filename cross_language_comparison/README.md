# Cross-Language Port Workflow

This repository now includes a Python-first cross-language scaffold for the MF12 spectral workflow while keeping MATLAB as the numerical ground truth.

Current status:

- `minimal_small` matches MATLAB at machine precision.
- `wavegroup_regression` and `benchmark_medium` now also match MATLAB at machine precision for both `eta` and `phi`.

## Folder Structure

- `cross_language_comparison/`: shared portable cases plus this workflow note
- `cross_language_comparison/cases/`: language-neutral validation and benchmark cases
- `python/`: Python package, tests, and helper scripts for the spectral-only MF12 path
- `cpp/`: C++ scaffold that reads the same shared case format
- `matlab/repro/`: MATLAB exporters that materialize shared cases and MATLAB-only reference artifacts

## Workflow

1. Generate or refresh shared cases from MATLAB:

```matlab
run('matlab/repro/export_cross_language_cases.m');
```

2. Verify one case from Python:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python -m mf12_python.cli verify cross_language_comparison/cases/minimal_small --repeats 1
```

3. Run Python benchmarks across all shared cases:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python python/scripts/run_benchmarks.py
```

4. Build the C++ case-inspection scaffold when a C++ compiler is available:

```powershell
cmake -S cpp -B cpp/build
cmake --build cpp/build
cpp/build/mf12_cpp inspect cross_language_comparison/cases/minimal_small
```

5. Run the current C++ order-2 verification path on the smallest shared case:

```powershell
cpp/build/mf12_cpp verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
```

6. Generate a Python-vs-C++-vs-MATLAB summary across the shared cases:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python python/scripts/compare_python_cpp_matlab.py cross_language_comparison/cases/minimal_small cross_language_comparison/cases/wavegroup_regression cross_language_comparison/cases/benchmark_medium
```

## Case And Result Format

Each case folder contains:

- `case.json`: scalar inputs, file references, and tolerances
- `inputs/*.csv`: `a`, `b`, `kx`, and `ky`
- `reference/matlab/*.csv`: MATLAB `eta`, `phi`, `x`, and `y`
- `reference/matlab/result.json`: MATLAB runtime metadata

Python verification outputs write:

- `eta.csv`
- `phi.csv`
- `x.csv`
- `y.csv`
- `result.json`

The C++ scaffold currently reads the same case format and now supports:

- a fully matching `minimal_small` run
- an in-progress third-order path for the larger cases

As the solver grows, it should keep writing the same result bundle shape as Python:

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

One important implementation trap showed up during the Python port and should be preserved for future C++ work:

- `Lambda3` is easy to mistranscribe because only its first grouped contribution is multiplied by the leading `h^2 / (4 * beta)` prefactor.
- The later `Fnpm`, `Fnpp`, `Fmpp`, `Gnpm`, `Gnpp`, and `Gmpp` terms are separate additive terms with their own `1 / beta` scaling.
- If those later terms are accidentally pulled under the same outer prefactor, the third-order `G_*` coefficients become much too large.
- That bug can make `eta` wrong while `phi` still looks correct, because the retained MF12 `phi` reconstruction depends primarily on the `F_*` / `mu_*` side, while the `eta` reconstruction uses the `G_*` / `Lambda3` side directly.

In practice, if a future port shows:

- `phi` matching MATLAB closely
- but third-order `eta`, especially the `n+m+p` superharmonic content, being badly wrong

then `Lambda3` transcription should be checked before anything else.
