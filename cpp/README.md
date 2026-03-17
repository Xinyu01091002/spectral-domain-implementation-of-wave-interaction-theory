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
- `benchmark <case_dir> [repeats] [warmup]`: time the C++ spectral path without requiring MATLAB reference files

Current implementation status:

- order `<= 3` is verified against the shared MATLAB reference cases
- superharmonic-only spectral path
- `minimal_small`, `wavegroup_regression`, and `benchmark_medium` pass against MATLAB reference outputs
- MATLAB, Python, and C++ now agree to machine precision on the shared spectral-only cases
- C++ reconstruction now prefers FFTW3 when available and otherwise falls back to the in-repo radix-2 inverse FFT path for power-of-two grids, with a direct inverse-DFT fallback for other sizes
- third-order coefficient profiling now separates `np2m`, `2npm`, and `npmpp`, and the `npmpp` hotspot can use optional OpenMP parallelization

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

Optional acceleration flags:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_OPENMP=ON
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_FFTW=ON
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_OPENMP=ON -DMF12_ENABLE_FFTW=ON
```

FFTW notes:

- the default `MF12_FFTW_ROOT` now points to `C:/tools/fftw3-mingw-dll`
- this layout expects `include/fftw3.h`, `lib/libfftw3.dll.a`, and `bin/libfftw3-3.dll`
- the project first tries `find_package(FFTW3)` and then falls back to this simple MinGW DLL layout
- when the DLL is found, CMake copies it next to `mf12_cpp.exe` after linking so the CLI can run without editing `PATH`

Optimization log:

- switched C++ reconstruction from the original direct inverse-DFT path to FFT-based reconstruction:
  - first with the in-repo radix-2 inverse FFT
  - later with FFTW3 on Windows/MinGW through a DLL import-library layout
  - result: reconstruction cost dropped noticeably, especially for larger grids
- added optional OpenMP parallelization only around the `npmpp` hotspot in the third-order coefficient build:
  - result: meaningful speedup on larger retained-component cases
  - this should be treated as part of the current C++ benchmark configuration, not as a single-thread baseline
- added profiling for:
  - linear coefficient time
  - second-order coefficient time
  - third-order coefficient time
  - `np2m`, `2npm`, and `npmpp` sub-breakdowns
  - result: `npmpp` was confirmed as the dominant hotspot
- tried more aggressive caching/reordering inside third-order coefficient generation:
  - result: some variants preserved `eta` but degraded `phi`
  - these were intentionally not kept because matching MATLAB to machine-precision scale is the hard constraint
- reduced non-essential `npmpp` bookkeeping during normal runs by avoiding debug index storage outside coefficient-dump workflows:
  - result: lower memory traffic and a small benchmark improvement without changing numerical output

Current interpretation:

- the largest remaining cost is still third-order coefficient generation, especially `npmpp`
- FFT optimization helps, but coefficient generation remains the main scaling bottleneck for large component counts
- C++ performance comparisons should always be reported together with whether OpenMP and FFTW were enabled

Possible future optimization directions:

- evaluate a benchmark-only streaming or chunked `npmpp` accumulation path to reduce peak memory use for very large cases
- continue optimizing `npmpp` data access patterns and temporary-array pressure without changing the numerical formula structure
- measure explicit thread-count scaling (`1, 2, 4, 8, ...`) so speedup claims are separated from algorithmic improvements
- if needed, investigate whether reconstruction can benefit further from FFTW planning strategy changes beyond the current simple integration
