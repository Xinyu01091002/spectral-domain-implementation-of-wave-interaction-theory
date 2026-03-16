# MATLAB Workspace

This directory is the MATLAB home for the repository.

Contents:

- `setup_paths.m`: adds the bundled MF12 source directory to the MATLAB path
- `irregularWavesMF12/Source/`: core MF12 direct and spectral implementation
- `examples/`: minimal runnable entry points
- `tests/`: fast regression and smoke checks
- `repro/`: scripts that export shared cases or reproduce representative paper outputs
- `verification/`: manuscript-facing comparisons, figures, and performance studies

Recommended starting points:

- `run('matlab/setup_paths.m')`
- `run('matlab/examples/run_spectral_minimal.m')`
- `run('matlab/tests/smoke_test_minimal.m')`
