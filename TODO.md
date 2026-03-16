# Repository Cleanup TODO

This checklist tracks the work needed to turn the current research workspace into a public-facing repository that is understandable, reproducible, and maintainable.

## Phase 1: Public entry points and documentation

- [x] Identify the two stable method pipelines used by the verification scripts:
  - direct: `mf12_direct_coefficients` + `mf12_direct_surface`
  - spectral: `mf12_spectral_coefficients` + `mf12_spectral_surface`
- [x] Rewrite the top-level README around these two workflows.
- [x] Add a simple path-setup entry point for new users.
- [x] Add minimal runnable examples for the direct and spectral pipelines.
- [x] Add software citation metadata (`CITATION.cff`).

## Phase 2: Repository structure cleanup

- [x] Separate public examples, tests, benchmarks, and manuscript reproduction more clearly.
- [x] Move or replace ambiguous generated artifacts currently stored in source locations.
- [x] Remove unreferenced top-level legacy scripts and isolated benchmark helpers.
- [ ] Decide whether the manuscript source should move under `paper/` or remain at the repository root.
- [ ] Decide whether the embedded `matlab/irregularWavesMF12/` code should remain in-place or move under `external/`.

## Phase 3: Naming and API surface

- [x] Review public-facing file names and mark unclear names for rename or deprecation.
- [x] Document which functions are primary public entry points versus legacy or benchmark variants.
- [x] Add thin wrappers if needed so users do not need to discover internal function variants by name alone.
- [x] Promote the clearer `mf12_*` names to the preferred implementation names in `matlab/irregularWavesMF12/Source`.

## Phase 4: Tests and reproducibility

- [x] Keep one fast smoke test that runs in seconds.
- [x] Add one regression test for direct-vs-spectral consistency.
- [x] Keep one representative manuscript reproduction entry point.
- [ ] Document expected output locations and runtime expectations.

## Phase 5: Release polish

- [x] Add contribution guidance for future collaborators.
- [ ] Update `CITATION.cff` so it reflects the current multi-language repository rather than the earlier MATLAB-only placeholder metadata.
- [ ] Refresh `CONTRIBUTING.md` for the current `matlab/`, `python/`, `cpp/`, and `cross_language_comparison/` workflow.
- [ ] Add release notes or a lightweight changelog.
- [ ] Verify the repository can be cloned and run from a clean MATLAB session.

## Phase 6: Cross-language implementations and performance benchmarking

- [ ] Decide the scope of the non-MATLAB ports:
  - C++ reference implementation of the core spectral reconstruction path
  - Python implementation for accessibility and analysis workflows
  - clarify whether the first target is spectral-only or both direct and spectral workflows
- [x] Define a shared benchmark case format so MATLAB, C++, and Python run the same inputs and grids.
- [x] Add fixed validation cases with archived expected outputs so cross-language comparisons are numerical, not only visual.
- [x] Implement the core MF12 spectral reconstruction in C++.
- [x] Implement the same reconstruction path in Python.
- [ ] Decide whether Python should begin as pure NumPy or include an accelerated path (for example `numba`, `pybind11`, or a C extension).
- [x] Add a cross-language verification script that compares MATLAB, C++, and Python outputs for `eta` and selected `phi` diagnostics.
- [ ] Add a reproducible benchmark harness that records:
  - coefficient-generation time
  - reconstruction time
  - total runtime
  - grid size and number of wave components
  - hardware and compiler/interpreter metadata
- [ ] Add a three-way speed comparison workflow for MATLAB, Python, and C++ using the shared cross-language cases.
- [ ] Decide and document which timing comparisons are primary:
  - MATLAB direct vs MATLAB spectral
  - MATLAB spectral vs Python spectral
  - MATLAB spectral vs C++ spectral
- [ ] Generate one archived reference benchmark summary and figure that report accuracy and speed for all three implementations on the shared cases.
- [ ] Archive benchmark outputs under `outputs/` and keep one reference log under `docs/benchmarks/`.
- [ ] Document build and run instructions for the C++ and Python implementations in the README once the first working versions exist.

## Cross-language implementation notes

- [x] Record the `Lambda3` transcription pitfall discovered during the Python port:
  - only the first grouped contribution carries the outer `h^2 / (4*beta)` prefactor
  - mis-grouping the later additive terms can break third-order `eta` while leaving `phi` apparently correct
  - check this first when porting the retained `G_np2m`, `G_2npm`, and `G_npmpp` branches to C++
- [x] Record the C++ third-order phase-correction pitfall discovered during the C++ port:
  - once third-order corrections update `omega`, the second self-harmonic phase array must also be refreshed as `omega2 = 2 * omega`
  - leaving `omega2` stale can produce small but persistent order-3 field mismatches even when the retained third-order coefficients already match Python/MATLAB
  - if future compiled ports match coefficient dumps but not the final order-3 fields, compare the pre-IFFT spectra and check stale derived phase arrays before re-debugging `Lambda3`
