# Spectral Domain Implementation for Crossing Sea Studies

This repository contains MATLAB scripts for generating and comparing nonlinear wave surfaces in crossing-sea conditions, with a focus on MF12-based spectral methods and validation workflows.

## Project Layout

```text
.
|- verification/                         # Validation/benchmark/plot scripts
|  |- compare_crossing_sea_methods.m
|  |- benchmark_mf12_direct_vs_spectral.m
|  |- benchmark_mf12_speed_memory.m
|  |- plot_mf12_theoretical_complexity_memory.m
|  |- plot_mf12_theoretical_three_methods.m
|  |- plot_phi3_direct_vs_spectral.m
|  |- plot_phi3_wavegroup_lines.m
|  `- plot_eta_wavegroup_lines.m
|- outputs/                              # Generated CSV/MAT/PNG/log outputs
|  |- figures/
|  `- processed_eta33/
|- generate_crossing_sea.m               # Crossing-sea case generation
|- analytic2D.m / compute_kxky.m         # Utility functions
|- my2nd_directional_generator.m         # Directional wave generator utility
|- irregularWavesMF12/                   # MF12 library (source + examples)
`- README.md
```

## Folder Convention

- Put validation, benchmark, and plotting scripts in `verification/`.
- Put generated outputs (`.csv/.mat/.png/.log`) in `outputs/`.
- Scripts in `verification/` are written to resolve paths from script location, so they can be run directly without switching MATLAB working directory.

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
run('verification/compare_crossing_sea_methods.m');
```

## Benchmark Commands

From repository root:

```matlab
% 1) Direct MF12 summation vs spectral reconstruction
run('verification/benchmark_mf12_direct_vs_spectral.m');

% 2) Multi-method benchmark (wrapper/streaming/spectral variants)
run('verification/benchmark_mf12_speed_memory.m');

% 3) Theory figures (no hardware timing dependence)
run('verification/plot_mf12_theoretical_complexity_memory.m');
run('verification/plot_mf12_theoretical_three_methods.m');

% 4) Third-order phi comparison figures
run('verification/plot_phi3_direct_vs_spectral.m');
run('verification/plot_phi3_wavegroup_lines.m');

% 5) Matching eta decomposition figures
run('verification/plot_eta_wavegroup_lines.m');
```

## Harmonic Decomposition Figures

- `verification/plot_phi3_wavegroup_lines.m` now generates two fixed paper-ready PNG files in `outputs/`:
  - `mf12_phi3_wavegroup_lines_comparison_pub.png`
  - `mf12_phi3_wavegroup_lines_error_pub.png`
- `verification/plot_eta_wavegroup_lines.m` generates in `outputs/`:
  - `mf12_eta_wavegroup_lines_comparison_pub.png`
  - `mf12_eta_wavegroup_lines_error_pub.png`

Both scripts use crossing-sea wave-group components and provide:
- 2x4 layout (top: centerline, bottom: diagonal)
- Columns: first harmonic, second superharmonic, second subharmonic, third superharmonic
- x-axis normalized by `\lambda_p`

## Recent Fixes

- Fixed third-order potential phase bug in:
  - `irregularWavesMF12/Source/surfaceMF12_new.m`
  - `irregularWavesMF12/Source/surfaceMF12.m`
- The `mu_2npm` contribution now uses `sin(theta_2npm)` (consistent phase variable).
- This removes the stable percent-level phase discrepancy previously observed in `phi` 3rd-order comparisons.

## Notes on Tracked vs Generated Files

- Source code (`*.m`) is tracked.
- Logs and generated outputs are ignored via `.gitignore`, including:
  - `outputs/compare_log.txt`
  - `outputs/run_output.txt`
  - `outputs/figures/`
  - `outputs/processed_eta33/`
  - generated PNG/CSV/MAT/log under `outputs/`

If you want to keep selected figures in GitHub, move them to a dedicated folder (for example `docs/images/`) and remove that folder from `.gitignore`.

## Suggested First Commit

```bash
git init
git add .
git commit -m "Initial commit: spectral-domain MATLAB project"
```
