# Tadah!MODELS

Physics and machine-learning components of the
[Tadah!](https://tadah.readthedocs.io/) MLIP toolchain. Tadah!MODELS
implements:

- **Regressors** — `M_KRR` (kernel ridge regression) and `M_BLR`
  (Bayesian linear regression).
- **Kernels** — `Kern_Linear`, `Kern_LQ`, `Kern_Polynomial`,
  `Kern_Quadratic`, `Kern_RBF`, `Kern_Sigmoid`.
- **Basis functions** — `BF_Linear`, `BF_Polynomial2`.
- **Two-body descriptors** — Behler-Parrinello (`D2_BP`), B-spline
  (`D2_Blip`), tabulated EAM (`D2_EAM`), Mie (`D2_MIE`), ZBL screening
  (`D2_ZBL`), Lennard-Jones (`D2_LJ`), composite (`D2_mJoin`).
- **Many-body descriptors** — embedded atom (`DM_EAM`), modified EAM
  (`DM_EAD`, `DM_mEAD`), B-spline (`DM_Blip`), polynomial
  (`DM_Poly`), composite (`DM_mJoin`).
- **Cutoff functions** — `Cut_Cos`, `Cut_CosS`, `Cut_Poly`,
  `Cut_PolyS`, `Cut_PT`, `Cut_Tanh`, `Cut_SinInv`.
- **Embedding functions** — `F_RLR`, `F_SQRT`, `F_Blip`, `F_mEnv`.

These classes are referenced by name in user-facing configuration
files; for the user view see
[descriptors](https://tadah.readthedocs.io/en/latest/descriptors.html),
[cutoffs](https://tadah.readthedocs.io/en/latest/cutoffs.html), and
[regressors](https://tadah.readthedocs.io/en/latest/regressors.html)
in the umbrella documentation.

## When you use this directly

Tadah!MODELS is not meant to be consumed standalone — Tadah!MLIP, MD,
and HPO all pull it in via `FetchContent`. Build it directly only if
you are developing new descriptors, kernels, or basis functions.

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

## License

GPL-3.0; see [LICENSE](LICENSE).

## Reporting issues

Bugs and feature requests:
[Tadah!MODELS issues](https://git.ecdf.ed.ac.uk/tadah/models/-/issues).
User-facing questions belong on
[Tadah!MLIP](https://git.ecdf.ed.ac.uk/tadah/tadah.mlip).
