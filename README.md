# Spectral Domain Implementation for Crossing Sea Studies

This repository contains MATLAB scripts for generating and comparing nonlinear wave surfaces in crossing-sea conditions, with a focus on MF12-based spectral methods and validation workflows.

## Project Layout

```text
.
|- compare_crossing_sea_methods.m        # Main comparison workflow (3rd-order methods)
|- generate_crossing_sea.m               # Crossing-sea case generation
|- check_correctness.m                   # Validation script for implementation consistency
|- check_speed.m                         # Speed/performance checks
|- compare_eta33_vwa_ow3d.m              # Eta33 comparison workflow
|- compare_eta33_vwa_ow3d_stream.m       # Streamlined eta33 comparison
|- test_superharmonic_fix.m              # Superharmonic fix test
|- analytic2D.m / compute_kxky.m         # Utility functions
|- my2nd_directional_generator.m         # Directional wave generator utility
|- irregularWavesMF12/                   # MF12 library (source + examples)
|- figures/                              # Generated figures (ignored by git)
`- processed_eta33/                      # Generated processed data (ignored by git)
```

## Requirements

- MATLAB (recommended R2021a or newer)
- No additional toolbox is strictly required for the included scripts

## Quick Start

From MATLAB in the repository root:

```matlab
addpath(genpath(fullfile(pwd, 'irregularWavesMF12')));
rehash toolboxcache;

% Generate a crossing-sea state
generate_crossing_sea;

% Run main method comparison
compare_crossing_sea_methods;
```

## Notes on Tracked vs Generated Files

- Source code (`*.m`) is tracked.
- Logs and generated outputs are ignored via `.gitignore`, including:
  - `compare_log.txt`
  - `run_output.txt`
  - `figures/`
  - `processed_eta33/`
  - root-level generated comparison PNG files

If you want to keep selected figures in GitHub, move them to a dedicated folder (for example `docs/images/`) and remove that folder from `.gitignore`.

## Suggested First Commit

```bash
git init
git add .
git commit -m "Initial commit: spectral-domain MATLAB project"
```
