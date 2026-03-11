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
- [ ] Decide whether the embedded `irregularWavesMF12/` code should remain in-place or move under `external/`.

## Phase 3: Naming and API surface

- [x] Review public-facing file names and mark unclear names for rename or deprecation.
- [x] Document which functions are primary public entry points versus legacy or benchmark variants.
- [x] Add thin wrappers if needed so users do not need to discover internal function variants by name alone.
- [x] Promote the clearer `mf12_*` names to the preferred implementation names in `irregularWavesMF12/Source`.

## Phase 4: Tests and reproducibility

- [x] Keep one fast smoke test that runs in seconds.
- [x] Add one regression test for direct-vs-spectral consistency.
- [x] Keep one representative manuscript reproduction entry point.
- [ ] Document expected output locations and runtime expectations.

## Phase 5: Release polish

- [x] Add contribution guidance for future collaborators.
- [ ] Add release notes or a lightweight changelog.
- [ ] Verify the repository can be cloned and run from a clean MATLAB session.
