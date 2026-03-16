# MF12 Spectral Reconstruction for Directional Wave Groups

This repository contains MATLAB code for nonlinear multidirectional wave reconstruction based on the third-order Madsen-Fuhrman 2012 (MF12) framework. The public interface is organized around two workflows:

- direct reconstruction: `mf12_direct_coefficients` + `mf12_direct_surface`
- spectral reconstruction: `mf12_spectral_coefficients` + `mf12_spectral_surface`

The repository is still research-oriented, but the entry points below are intended to let a new user understand the codebase quickly and run a representative case from a clean MATLAB session.

## Recent Updates

Recent cleanup and restructuring focused on making the repository easier to understand and easier to run as a public codebase:

- standardized the preferred MF12 function names in `matlab/irregularWavesMF12/Source` around
  `mf12_direct_coefficients`, `mf12_direct_surface`,
  `mf12_spectral_coefficients`, and `mf12_spectral_surface`
- removed legacy wrapper-style interfaces and unused bundled examples that were no longer part of the active workflow
- simplified `matlab/setup_paths.m` so it only adds the MF12 source directory required by the current scripts
- added clearer MATLAB entry points under `matlab/examples/`, `matlab/tests/`, `matlab/repro/`, and `matlab/verification/`
- added a spectral-only 3D section visualization for crossing directional wave groups in
  `matlab/verification/plot_phi3_wavegroup_surface_sections.m`

## Requirements

- MATLAB R2021a or newer is recommended
- no nonstandard MATLAB toolboxes are currently required by the included entry points

## Quick Start

From the repository root in MATLAB:

```matlab
run('matlab/setup_paths.m');
```

Run one minimal example for each supported workflow:

```matlab
run('matlab/examples/run_direct_minimal.m');
run('matlab/examples/run_spectral_minimal.m');
```

Run the release-readiness smoke test:

```matlab
run('matlab/tests/smoke_test_minimal.m');
```

Run one representative manuscript-facing reproduction script:

```matlab
run('matlab/repro/reproduce_manuscript_minimal.m');
```

## Where To Start

- new user: run `matlab/examples/run_direct_minimal.m` or `matlab/examples/run_spectral_minimal.m`
- reviewer or reproducer: start with `matlab/repro/reproduce_manuscript_minimal.m`
- contributor: run `matlab/tests/smoke_test_minimal.m` before and after changes
- manuscript-facing investigation: look under `matlab/verification/`

## Main Workflows

### Direct reconstruction

The direct path builds MF12 coefficients in physical form and evaluates the reconstruction directly on a spatial grid:

```matlab
[X, Y] = meshgrid(x, y);
coeffs = mf12_direct_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta, phi] = mf12_direct_surface(order, coeffs, X, Y, t);
```

This is the reference path used for direct-vs-spectral comparisons.

### Spectral reconstruction

The spectral path builds the superharmonic-oriented coefficient structure and reconstructs the field through FFT-based accumulation:

```matlab
coeffs = mf12_spectral_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta, phi] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);
```

This is the primary high-performance path used in the benchmark and large-`N` verification scripts.

## Repository Layout

```text
.
|- matlab/                MATLAB implementation, examples, tests, and verification scripts
|- python/                Python spectral port, tests, and helper scripts
|- c/                     reserved home for the future C or C++ implementation
|- cross_language_comparison/  shared cases, comparison workflow notes, and language-neutral metadata
|- docs/                  supporting repository notes and archived reference logs
|- outputs/               generated outputs, ignored by git
|- TODO.md                repository cleanup checklist
|- LICENSE
`- README.md
```

## Recommended Entry Points

- `matlab/setup_paths.m`: add the MF12 source directory to the MATLAB path
- `matlab/irregularWavesMF12/Source/mf12_direct_coefficients.m`, `matlab/irregularWavesMF12/Source/mf12_direct_surface.m`: preferred direct-workflow entry points
- `matlab/irregularWavesMF12/Source/mf12_spectral_coefficients.m`, `matlab/irregularWavesMF12/Source/mf12_spectral_surface.m`: preferred spectral-workflow entry points
- `matlab/examples/run_direct_minimal.m`: smallest direct reconstruction example
- `matlab/examples/run_spectral_minimal.m`: smallest spectral reconstruction example
- `matlab/tests/smoke_test_minimal.m`: fast direct-vs-spectral consistency check
- `matlab/tests/regression_wavegroup_phi3.m`: representative directional wave-group regression check
- `matlab/repro/reproduce_manuscript_minimal.m`: minimal manuscript reproduction

## Verification and Benchmark Scripts

The scripts under `matlab/verification/` are manuscript-facing and remain closer to the research workflow than to a polished package API. The most useful ones for understanding current usage are:

- `plot_phi3_wavegroup_lines.m`
- `plot_phi3_wavegroup_surface_sections.m`
- `plot_eta_wavegroup_lines.m`
- `compare_crossing_sea_methods.m`
- `benchmark_direct_vs_spectral_vsN.m`
- `benchmark_mf12_speed_memory.m`
- `test_new_spectral_realistic_sea_highN.m`

These scripts are the best references for how the two reconstruction pipelines are used in realistic directional and crossing-sea cases.

In particular:

- `plot_phi3_wavegroup_lines.m` compares direct and spectral reconstructions for directional or crossing wave-group sections
- `plot_phi3_wavegroup_surface_sections.m` provides a spectral-only 3D visualization with centerline/diagonal section overlays for focused crossing wave groups

## Relationship to Bundled MF12 Code

The directory `matlab/irregularWavesMF12/` contains the bundled MF12 source implementation used by this repository. This repository adds public-facing examples, smoke tests, manuscript reproduction entry points, and verification scripts centered on direct versus spectral reconstruction.

The `mf12_*` function names in `matlab/irregularWavesMF12/Source` are the preferred public and implementation names in this repository.

Inside `matlab/irregularWavesMF12/Source`, the preferred clear-name aliases are:

- `mf12_direct_coefficients`
- `mf12_spectral_coefficients`
- `mf12_direct_surface`
- `mf12_spectral_surface`

If the repository is prepared for a more formal release later, this boundary should likely be made more explicit, for example by moving the bundled code under a dedicated `external/` or `third_party/` location.

## Reproducibility Notes

- generated figures and logs should go under `outputs/`
- a sample benchmark log is archived in `docs/benchmarks/`
- the current reproduction entry point intentionally covers one representative manuscript case rather than every figure in the paper

## Cross-Language Scaffold

A Python-first spectral-only cross-language scaffold now lives under `python/`, with shared MATLAB-ground-truth cases under `cross_language_comparison/cases/` and the exporter at `matlab/repro/export_cross_language_cases.m`.

See `cross_language_comparison/README.md` for the shared case format, verification flow, and benchmark entry points.

## Citation and License

- license: MIT, see `LICENSE`
- software citation metadata: `CITATION.cff`

Until the citation metadata is finalized, update the placeholder author and repository information in `CITATION.cff` before public release.
