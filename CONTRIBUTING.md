# Contributing

This repository is still organized around an active research workflow, so contributions should stay conservative and reproducibility-focused.

## Expectations

- preserve the two main method pipelines already used by the verification scripts:
  - direct: `mf12_direct_coefficients` + `mf12_direct_surface`
  - spectral: `mf12_spectral_coefficients` + `mf12_spectral_surface`
- prefer the clearer public `mf12_*` names from `matlab/irregularWavesMF12/Source` in new user-facing examples and documentation
- avoid changing scientific behavior unless the change is documented and validated
- keep generated figures, logs, and temporary outputs out of tracked source directories

## Preferred change types

- documentation and entry-point cleanup
- tests and reproducibility improvements
- performance improvements that do not change numerical results
- bug fixes with a reproducible verification case

## Before submitting changes

- run the smallest relevant entry point you touched
- run `matlab/tests/smoke_test_minimal.m` for general release-readiness
- if you changed a reconstruction path, also run `matlab/tests/regression_wavegroup_phi3.m`

## Notes

- update `CITATION.cff` before public release
- do not move `manuscript.tex` or the bundled `matlab/irregularWavesMF12/` tree without checking downstream references first
- delete files only after confirming they are not used by `README`, `matlab/examples`, `matlab/tests`, `matlab/repro`, or the retained `matlab/verification` scripts
