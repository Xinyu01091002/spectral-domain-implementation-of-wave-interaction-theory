# Documentation Guide

This `docs/` folder is for reference material that supports the repository but is not part of the first-run setup flow.

Use the nearest README for practical onboarding:

- repository overview and quick starts: `README.md`
- MATLAB-specific usage: `matlab/README.md`
- Python-specific usage: `python/README.md`
- C++ build and CLI usage: `cpp/README.md`
- shared MATLAB/Python/C++ workflow: `cross_language_comparison/README.md`

Use `docs/` for deeper or more archival material:

- `API_OVERVIEW.md`: cross-cutting API and interface map across MATLAB, Python, and C++
- `benchmarks/`: archived benchmark logs and reproducibility records

## What Belongs In `docs/`

Good fits for this folder:

- architecture notes
- numerical-method caveats and compatibility notes
- validation methodology
- benchmark interpretation and archived logs
- porting notes that are useful for maintainers but too detailed for quick-start READMEs

Usually keep these out of `docs/` and in the nearest README instead:

- installation steps
- build commands
- first-run examples
- common CLI commands
- "start here" guidance for new users

## Practical Rule

If a new user needs the information in their first 10 minutes, keep it in a README close to the code they are using.

If the information is helpful but secondary, historical, or cross-cutting, `docs/` is the right home.
