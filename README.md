# Spectral Domain Implementation for Crossing Sea Studies

This repository contains MATLAB scripts for generating and comparing nonlinear wave surfaces in crossing-sea conditions, with a focus on MF12-based spectral methods and validation workflows.

## Project Layout

```text
.
|- compare_crossing_sea_methods.m        # Main comparison workflow (3rd-order methods)
|- generate_crossing_sea.m               # Crossing-sea case generation
|- benchmark_mf12_direct_vs_spectral.m   # Direct vs spectral speed/error benchmark
|- benchmark_mf12_speed_memory.m         # Multi-method speed/memory benchmark
|- plot_mf12_theoretical_complexity_memory.m  # End-to-end theory: direct vs spectral (time+memory)
|- plot_mf12_theoretical_three_methods.m      # Theory: direct vs spectral vs streaming
|- plot_phi3_direct_vs_spectral.m        # 2D phi(3rd) direct vs spectral comparison
|- plot_phi3_wavegroup_lines.m           # Crossing-sea wave-group phi_s harmonic decomposition plots
|- plot_eta_wavegroup_lines.m            # Crossing-sea wave-group eta harmonic decomposition plots
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

% 3) Theory figures (no hardware timing dependence)
plot_mf12_theoretical_complexity_memory;
plot_mf12_theoretical_three_methods;

% 4) Third-order phi comparison figures
plot_phi3_direct_vs_spectral;
plot_phi3_wavegroup_lines;

% 5) Matching eta decomposition figures
plot_eta_wavegroup_lines;
```

## Harmonic Decomposition Figures

- `plot_phi3_wavegroup_lines.m` now generates two fixed paper-ready PNG files:
  - `mf12_phi3_wavegroup_lines_comparison_pub.png`
  - `mf12_phi3_wavegroup_lines_error_pub.png`
- `plot_eta_wavegroup_lines.m` generates:
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
