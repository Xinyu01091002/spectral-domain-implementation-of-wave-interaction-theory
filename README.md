# Spectral Domain Implementation for Crossing Sea Studies

This repository contains MATLAB scripts for generating and comparing nonlinear wave surfaces in crossing-sea conditions, with a focus on MF12-based spectral methods and validation workflows.

## Project Layout

```text
.
|- compare_crossing_sea_methods.m        # Main comparison workflow (3rd-order methods)
|- generate_crossing_sea.m               # Crossing-sea case generation
|- benchmark_mf12_direct_vs_spectral.m   # Direct vs spectral speed/error benchmark
|- benchmark_mf12_speed_memory.m         # Multi-method speed/memory benchmark
|- benchmark_mf12_plot.m                 # Plot benchmark results from latest MAT file
|- plot_phi3_direct_vs_spectral.m        # 2D phi(3rd) direct vs spectral comparison
|- plot_phi3_wavegroup_lines.m           # Wave-group phi(3rd) center/diagonal line comparison
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

## Benchmark Commands

From repository root:

```matlab
% 1) Direct MF12 summation vs spectral reconstruction
benchmark_mf12_direct_vs_spectral;

% 2) Multi-method benchmark (wrapper/streaming/spectral variants)
benchmark_mf12_speed_memory;
benchmark_mf12_plot;

% 3) Third-order phi comparison figures
plot_phi3_direct_vs_spectral;
plot_phi3_wavegroup_lines;
```

## Recent Fixes

- Fixed third-order potential phase bug in:
  - `irregularWavesMF12/Source/surfaceMF12_new.m`
  - `irregularWavesMF12/Source/surfaceMF12.m`
- The `mu_2npm` contribution now uses `sin(theta_2npm)` (consistent phase variable).
- This removes the stable percent-level phase discrepancy previously observed in `phi` 3rd-order comparisons.

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
