# Python Port

This directory contains the Python spectral MF12 port, including constant-z kinematics reconstruction.

Layout:

- `src/mf12_python/`: implementation and CLI
- `tests/`: Python unit tests
- `scripts/`: plotting and benchmark helpers
- `pyproject.toml`: package metadata for editable installs

Recommended setup from the repository root:

```powershell
python -m pip install -e python
```

If you do not want to install it, use:

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

Generated outputs default to `outputs/cross_language_comparison/`.

The shared-case verification path now checks:

- surface fields: `eta`, `phi`
- constant-z kinematics: `u`, `v`, `w`, `p`, `phi_vol`, `uV`, `vV`, `a_x`, `a_y`
