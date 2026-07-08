# Repository Guidelines

## Project Overview

`CPP/` is the complete and authoritative implementation for calculating URFs from streamlines. Focus development, fixes, validation, and build guidance on this directory.

Other top-level files are reference material only:

- `matlabPrototype.m` is a historical/prototype reference.
- `RUST/` is not the maintained implementation for current work.
- `test_data/` contains sample inputs, including `npsat_urf.opt` and trajectory data.

Do not spend effort on the Rust or MATLAB code unless the user explicitly asks for it.

## Build and Run

Run CMake commands from `CPP/`. The C++ implementation depends on Eigen headers only. If Eigen is not installed under `/usr/include/eigen3`, provide `EIGEN_PATH`:

```bash
cd CPP
cmake -S . -B build -DEIGEN_PATH=/path/to/eigen3
cmake --build build
```

The C++ binary expects a process id argument for normal execution and supports `-v` for version output.

## Coding Conventions

- Prefer minimal, local changes that preserve the existing C++ data flow and file formats.
- Preserve existing option names and output column names unless the task is explicitly about changing the interface.
- Keep the code compatible with GCC 5-era systems and no later than C++11.
- The C++ code uses headers for shared structures/helpers and Eigen sparse matrix types; follow that pattern for nearby changes.
- Avoid large formatting-only changes in legacy C++ files unless requested.

## Validation

When changing C++ code, run the CMake build if Eigen is available:

```bash
cd CPP
cmake -S . -B build -DEIGEN_PATH=/path/to/eigen3
cmake --build build
```

If Eigen is not available in the environment, state that clearly in the final response.

For behavior changes, use files in `test_data/` as sample inputs and document any generated output files.

## Git and Generated Files

- Do not commit build directories, generated `.dat` outputs, or local dependency checkouts.
- Check `git status --short` before editing and before finishing.
- Do not overwrite user changes; if a file has unrelated edits, work around them or ask for direction.
