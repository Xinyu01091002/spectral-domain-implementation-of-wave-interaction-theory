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
- `run('matlab/verification/plot_phi3_direct_vs_spectral.m')`
- `run('matlab/verification/plot_kinematics_direct_vs_spectral.m')`

Current MATLAB status:

- `mf12_direct_coefficients` + `mf12_direct_surface` remain the physical-space reference workflow.
- `mf12_spectral_coefficients` + `mf12_spectral_surface` reproduce the validated surface fields to machine precision on the retained MATLAB verification cases.
- `mf12_spectral_kinematics` reproduces `kinematicsMF12` on constant-`z` planes to machine precision for `u`, `v`, `w`, `p`, `phi`, `uV`, `vV`, `a_x`, and `a_y`.
