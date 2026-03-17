# C++ Port

This directory contains the C++ spectral MF12 implementation and CLI.

The C++ port uses the same shared-case format as the Python tools and treats the archived MATLAB outputs under `cross_language_comparison/cases/` as the numerical ground truth.

## Layout

- `CMakeLists.txt`: build entry point for the C++ CLI
- `include/mf12_cpp/`: headers for case I/O, comparison, and spectral reconstruction
- `src/`: CLI, implementation, and result-writing code
- `build/`: generated build artifacts

## Supported CLI Commands

The CLI entry point is `cpp/build/mf12_cpp.exe` on Windows and `cpp/build/mf12_cpp` on Linux.

- `inspect <case_dir>`: load a shared case and print a summary
- `validate <case_dir>`: check that required arrays and reference files are present and dimensionally consistent
- `verify <case_dir> [output_dir] [order]`: run the C++ spectral path and compare against MATLAB reference outputs
- `benchmark <case_dir> [repeats] [warmup]`: time the C++ spectral path without requiring MATLAB reference files
- `dump-coeffs <case_dir> <output_dir> [order]`: dump coefficient arrays for debugging and cross-language comparison work

## Current Implementation Status

- the current port is spectral-only
- order `<= 3` is verified against the shared MATLAB reference cases
- `minimal_small`, `wavegroup_regression`, and `benchmark_medium` pass against MATLAB reference outputs
- MATLAB, Python, and C++ agree to machine precision on the shared spectral regression cases
- reconstruction prefers FFTW3 when available and otherwise falls back to the in-repo radix-2 inverse FFT path for power-of-two grids, with a direct inverse-DFT fallback for other sizes
- third-order coefficient profiling separates `np2m`, `2npm`, and `npmpp`, and the `npmpp` hotspot can use optional OpenMP parallelization

## Prerequisites

- a C++17-capable compiler such as `g++` 9+ or a recent MinGW-w64 / WinLibs / MSYS2 toolchain
- CMake 3.18 or newer
- Ninja is recommended for the commands below
- FFTW3 is strongly recommended for better FFT performance
- OpenMP comes from the compiler/runtime, not from this repository

FFTW behavior:

- when FFTW3 is available and `MF12_ENABLE_FFTW=ON`, the CLI links to FFTW
- when FFTW is not found, the code falls back to the in-repo radix-2 inverse FFT for power-of-two grids and a direct inverse-DFT fallback otherwise

## Installation / Build

### Ubuntu / Debian-like Linux

Install dependencies:

```bash
sudo apt update
sudo apt install -y build-essential cmake ninja-build libfftw3-dev pkg-config
```

Verify the toolchain:

```bash
cmake --version
g++ --version
pkg-config --modversion fftw3
```

Configure and build from the repository root:

```bash
cmake -S cpp -B cpp/build -G Ninja -DMF12_ENABLE_OPENMP=ON -DMF12_ENABLE_FFTW=ON
cmake --build cpp/build
```

Run a quick verification command:

```bash
./cpp/build/mf12_cpp verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp_linux/minimal_small
```

### Windows

Recommended toolchain:

- use MinGW-w64, WinLibs, or MSYS2 with `g++`
- install CMake and Ninja and ensure they are available on `PATH`
- install FFTW3 system-wide, or place it in a local folder and pass `-DMF12_FFTW_ROOT=C:/path/to/fftw`
- OpenMP support is provided by the chosen compiler/runtime

Configure and build from the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_OPENMP=ON -DMF12_ENABLE_FFTW=ON
cmake --build cpp/build
```

If FFTW is installed in a non-default location:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_OPENMP=ON -DMF12_ENABLE_FFTW=ON -DMF12_FFTW_ROOT=C:/tools/fftw3-mingw-dll
cmake --build cpp/build
```

Run a quick verification command:

```powershell
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
```

## Build And Run

From the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build cpp/build
cpp/build/mf12_cpp.exe verify cross_language_comparison/cases/minimal_small outputs/cross_language_comparison/verify_cpp/minimal_small
```

Useful commands:

```powershell
cpp/build/mf12_cpp.exe inspect cross_language_comparison/cases/minimal_small
cpp/build/mf12_cpp.exe validate cross_language_comparison/cases/minimal_small
cpp/build/mf12_cpp.exe benchmark cross_language_comparison/cases/benchmark_medium 3 1
cpp/build/mf12_cpp.exe dump-coeffs cross_language_comparison/cases/wavegroup_regression outputs/cross_language_comparison/dump_coeffs/wavegroup_regression 3
```

## Optional Acceleration

Enable OpenMP and FFTW from the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_OPENMP=ON
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_FFTW=ON
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DMF12_ENABLE_OPENMP=ON -DMF12_ENABLE_FFTW=ON
```

FFTW notes:

- the default `MF12_FFTW_ROOT` points to `C:/tools/fftw3-mingw-dll`
- that layout expects `include/fftw3.h`, `lib/libfftw3.dll.a`, and `bin/libfftw3-3.dll`
- the project first tries `find_package(FFTW3)` and then falls back to this simple MinGW DLL layout
- when the DLL is found, CMake copies it next to `mf12_cpp.exe` after linking so the CLI can run without editing `PATH`

Parallel control:

- when OpenMP is enabled, the runtime thread count can be controlled without recompiling
- on Windows PowerShell:

```powershell
$env:OMP_NUM_THREADS = "8"
cpp/build/mf12_cpp.exe benchmark cross_language_comparison/cases/benchmark_dense_300 1 1
```

- on Linux or macOS shells:

```bash
OMP_NUM_THREADS=8 ./cpp/build/mf12_cpp benchmark cross_language_comparison/cases/benchmark_dense_300 1 1
```

## Troubleshooting / Compatibility

- very old environments, including Ubuntu 16.04-era `cmake` and `g++`, may be too old for this C++17 build
- if OpenMP is not detected, the build can still succeed, but the hotspot loops run without compiler-managed parallelism and performance may drop
- if FFTW is not found, the project does not fail by default; it falls back to slower internal reconstruction paths
- the shared-case inputs are portable across platforms, but binaries are not; rebuild on each target OS

## Porting Notes

One porting pitfall worth preserving:

- after the third-order frequency correction updates `omega`, the second self-harmonic phase must also be refreshed through `omega2 = 2 * omega`
- leaving `omega2` stale can make the pre-IFFT spectra differ even when the third-order coefficients themselves already match
- this showed up as a small but visible order-3 field mismatch in `wavegroup_regression` while `benchmark_medium` was less sensitive

Performance interpretation:

- the largest remaining cost is still third-order coefficient generation, especially `npmpp`
- FFT optimization helps, but coefficient generation remains the main scaling bottleneck for large component counts
- performance comparisons should always be reported together with whether OpenMP and FFTW were enabled
