# Python Port

This directory contains the Python spectral-only MF12 port.

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
```

Generated outputs default to `outputs/cross_language_comparison/`.
