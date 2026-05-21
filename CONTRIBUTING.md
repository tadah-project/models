# Contributing to Tadah!MODELS

The full multi-repo contributor guide for the Tadah! project lives on
the umbrella repository:
<https://git.ecdf.ed.ac.uk/tadah/tadah.mlip/-/blob/develop/CONTRIBUTING.md>

This file covers what is specific to Tadah!MODELS. Bug reports and
feature requests on this module use the GitLab issue tracker on this
project; cross-cutting changes (work that touches more than one
Tadah! module) should be discussed on the umbrella project first.

## Build & test in isolation

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
ctest --test-dir build/tests --output-on-failure
```

In-source builds are blocked deliberately; always pick a separate
build directory.

### Dependencies

- C++17 compiler.
- [Tadah!CORE](https://git.ecdf.ed.ac.uk/tadah/core) — fetched
  automatically.
- OpenMP (optional; enabled with `-DTADAH_ENABLE_OPENMP=ON`).
- LAPACK / BLAS at link time (provided by the build environment).

## Layout

| Path | Purpose |
|---|---|
| `include/tadah/models/` | Public headers; one per public class. |
| `src/`                  | Implementation translation units. |
| `tests/`                | Unit tests for every descriptor, kernel, and basis function. |

## Adding a new descriptor

The general recipe is

1. Implement a class derived from `D2_Base` or `DM_Base` in
   `include/tadah/models/` and `src/`.
2. Register the new class in `d2_all.h` / `dm_all.h` so the factory
   can dispatch on the configuration name.
3. Add a test under `tests/`.
4. Document the new class in
   [Tadah.MLIP/docs/source/descriptors.rst](https://git.ecdf.ed.ac.uk/tadah/tadah.mlip/-/blob/develop/docs/source/descriptors.rst).
