# C++ Port

This directory now contains the initial C++ scaffold for the MF12 spectral workflow.

Current contents:

- `CMakeLists.txt`: small build entry point for the C++ CLI
- `include/mf12_cpp/`: headers for shared-case and CSV loading
- `src/`: CLI and implementation files

The first milestone remains spectral-only and should reuse:

- `cross_language_comparison/cases/` for shared benchmark and validation inputs
- the same result schema used by the Python verification flow
- the MATLAB outputs under each shared case as the numerical ground truth

Current CLI scope:

- `inspect <case_dir>`: load a shared case and print a summary
- `validate <case_dir>`: check that required arrays and reference files are present and dimensionally consistent
- `verify <case_dir> [output_dir]`: run the current C++ spectral path and compare against MATLAB reference outputs

Current implementation status:

- order `<= 3` is verified against the shared MATLAB reference cases
- superharmonic-only spectral path
- `minimal_small`, `wavegroup_regression`, and `benchmark_medium` pass against MATLAB reference outputs
- MATLAB, Python, and C++ now agree to machine precision on the shared spectral-only cases

Recent porting pitfall that was worth recording:

- after the third-order frequency correction updates `omega`, the second self-harmonic phase must also be refreshed through `omega2 = 2 * omega`
- leaving `omega2` stale can make the pre-IFFT spectra differ even when the third-order coefficients themselves already match
- this showed up as a small but visible order-3 field mismatch in `wavegroup_regression` while `benchmark_medium` happened to remain much less sensitive

Build commands from the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build cpp/build
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
```
