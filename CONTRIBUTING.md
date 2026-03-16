# Contributing

This repository is still research-oriented, so changes should stay conservative, reproducible, and easy to validate across languages.

## Working Principles

- preserve the scientific behavior of the retained MF12 workflows unless a change is intentional, documented, and validated
- treat MATLAB as the numerical ground truth for the shared cross-language cases
- preserve attribution to the DTU MATLAB implementation released by David R. Fuhrman, since the coefficient workflow in this repository follows and is verified against that implementation
- keep user-facing examples and documentation aligned with the current folder layout:
  - `matlab/`
  - `python/`
  - `cpp/`
  - `cross_language_comparison/`
- keep generated outputs, figures, logs, and build artifacts out of tracked source directories

## Main Workflows

The primary workflows in this repository are:

- MATLAB direct reconstruction:
  - `mf12_direct_coefficients`
  - `mf12_direct_surface`
- MATLAB spectral reconstruction:
  - `mf12_spectral_coefficients`
  - `mf12_spectral_surface`
- Python spectral reconstruction and verification:
  - `python/src/mf12_python/`
- C++ spectral reconstruction and verification:
  - `cpp/src/`

For new examples and docs, prefer the clearer `mf12_*` names from [matlab/irregularWavesMF12/Source](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/irregularWavesMF12/Source).

## Good Contribution Types

- documentation cleanup and onboarding improvements
- tests and reproducibility improvements
- bug fixes with a clear validation case
- performance improvements that do not change numerical results
- cross-language verification improvements

## Before You Submit Changes

Run the smallest relevant checks for the area you touched.

If you changed MATLAB core logic, run:

```matlab
run('matlab/setup_paths.m');
run('matlab/tests/smoke_test_minimal.m');
run('matlab/tests/regression_wavegroup_phi3.m');
```

If you changed the shared case export or MATLAB reference outputs, also run:

```matlab
run('matlab/repro/export_cross_language_cases.m');
```

If you changed Python code, run from the repository root:

```powershell
python -m pip install -e python
python -m unittest discover python/tests
python -m mf12_python.cli verify cross_language_comparison/cases/minimal_small --repeats 1
```

If you changed C++ code, run from the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build cpp/build
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
```

If you changed shared comparison logic, figures, or case definitions, run at least one cross-language comparison:

```powershell
$env:PYTHONPATH = (Resolve-Path 'python/src').Path
python python/scripts/compare_python_cpp_matlab.py cross_language_comparison/cases/minimal_small cross_language_comparison/cases/wavegroup_regression cross_language_comparison/cases/benchmark_medium
```

## File and Folder Hygiene

- keep generated files under `outputs/`
- keep C++ build artifacts under `cpp/build/`
- do not commit temporary debug dumps unless they are intentionally promoted to a retained artifact
- when deleting files, first confirm they are not referenced by:
  - [README.md](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/README.md)
  - [matlab/examples](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/examples)
  - [matlab/tests](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/tests)
  - [matlab/repro](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/matlab/repro)
  - [python/scripts](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/python/scripts)
  - [cross_language_comparison](C:/Research/spectral%20domain%20implementation%20of%20wave%20interaction%20theory/cross_language_comparison)

## Notes for Future Porting Work

- when documenting or presenting the implementation, acknowledge the DTU MATLAB release by David R. Fuhrman:
  https://data.dtu.dk/articles/code/Matlab_implementation_of_the_third-order_multi-directional_irregular_wave_theory_of_Madsen_Fuhrman_2012_J_Fluid_Mech_698_304-334_/22060124
- if third-order `phi` looks correct but third-order `eta` is wrong, check the `Lambda3` transcription first
- if third-order coefficients match but final order-3 fields still differ in a compiled port, check stale derived phase arrays such as `omega2 = 2 * omega`
- when comparing diagonal lines, make sure all languages use the same archived sampling path from the MATLAB reference
