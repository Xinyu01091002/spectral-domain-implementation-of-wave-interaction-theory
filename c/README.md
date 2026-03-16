# C Port Placeholder

This directory is reserved for the future C or C++ implementation of the MF12 spectral workflow.

When that port starts, it should reuse:

- `cross_language_comparison/cases/` for shared benchmark and validation inputs
- the same result schema used by the Python verification flow
- the MATLAB outputs under each shared case as the numerical ground truth

The current recommendation is to begin with a small standalone CLI for the spectral-only path before introducing a library interface.
