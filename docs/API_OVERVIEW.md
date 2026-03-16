# API Overview

This repository currently has two public reconstruction workflows that should be treated as the primary supported interfaces.

## Primary workflows

### Direct reconstruction

```matlab
[X, Y] = meshgrid(x, y);
coeffs = mf12_direct_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta, phi] = mf12_direct_surface(order, coeffs, X, Y, t);
```

Use this path when you want a direct physical-space reference calculation.

### Spectral reconstruction

```matlab
coeffs = mf12_spectral_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta, phi] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);
```

Use this path for FFT-based reconstruction and large-`N` performance studies.

## Supporting entry points

- `matlab/setup_paths.m`: adds the bundled MF12 source tree to the MATLAB path
- `matlab/irregularWavesMF12/Source/mf12_direct_coefficients.m`: preferred public name for direct coefficient generation
- `matlab/irregularWavesMF12/Source/mf12_direct_surface.m`: preferred public name for direct reconstruction
- `matlab/irregularWavesMF12/Source/mf12_spectral_coefficients.m`: preferred public name for spectral coefficient generation
- `matlab/irregularWavesMF12/Source/mf12_spectral_surface.m`: preferred public name for spectral reconstruction
- `matlab/examples/run_direct_minimal.m`: smallest direct example
- `matlab/examples/run_spectral_minimal.m`: smallest spectral example
- `matlab/tests/smoke_test_minimal.m`: fast release-readiness check
- `matlab/tests/regression_wavegroup_phi3.m`: more representative directional wave-group comparison
- `matlab/repro/reproduce_manuscript_minimal.m`: minimal manuscript reproduction

## Source-level implementation names

Within `matlab/irregularWavesMF12/Source`, the clearer source-level names are:

- `mf12_direct_coefficients`
- `mf12_spectral_coefficients`
- `mf12_direct_surface`
- `mf12_spectral_surface`

The preferred implementation files now live under these clear names.

## Internal specialized files

The bundled MF12 source tree still contains a small number of specialized internal files such as:

- `surfaceMF12_integrated_opt.m`

These are not part of the recommended public interface.
