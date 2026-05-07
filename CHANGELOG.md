# Changelog

All notable changes to **Tadah.MODELS** (descriptors, cutoffs, basis functions, regressors) are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

Nothing yet.

## [1.3.0-beta.1] — TBD

The largest module-level surface change of this release: a new family of rotationally invariant power-spectrum descriptors (`DM_Poly` and friends), new `F_Blip` / `DM_mBlip` / `DM_mBlipEnv` embedding-function descriptors, in-place HPO `reinit()` plumbing, and a corrected `Cut_PT` normalisation.

### Added
- `DM_Poly<ChannelFn>` — rotationally invariant power-spectrum descriptor with blip radial basis. Registered variants:
  - `DM_PS_SQ`, `DM_PS_SQRT`, `DM_PS_RLR`, `DM_PS_Blip` (channel = power spectrum).
  - `DM_OW_*` and `DM_GOW_*` variants (`AggMode` template parameter — orbital-wise / general orbital-wise aggregation).
  - `DM_GPS_SQ`, `DM_GPS_SQRT`, `DM_GPS_RLR`, `DM_GPS_Blip` (`RadialBasis` template parameter — `GaussianRadial` instead of `BlipRadial`).
- `radial_policies.h` — `BlipRadial` and `GaussianRadial` policies that drive the `RadialBasis` template parameter.
- `F_Blip` embedding function and `DM_mBlip` descriptor.
- `DM_mBlipEnv` — inverse-envelope variant of `DM_mBlip` that wraps the `F_Blip` embedding in a `Cut_SinInv` envelope (controlled by the new `RCUTENV` blueprint key).
- `descriptor::reinit()` and `DC_Selector::reinit()` — in-place parameter updates for HPO.
- Single-ρ support and generalisation of the `F_` (embedding) function family.
- DM_Poly test suite.
- Doxygen / LaTeX documentation for basis functions and `DM_Blip`.
- `TADAH_ENABLE_ASAN` CMake option.

### Changed
- **`Cut_PT` normalised so that f(0) = 1.** Previously the polynomial cutoff was off by a constant, which silently rescaled descriptors that used it. Energies/forces/stresses computed with `Cut_PT` against v1.1.1 will differ. See the umbrella migration guide.
- `Kind` enum and `val_at_rcut()` added to the `Cut_Base` interface (used by load-time validation and per-species bias precomputation).
- `Cut_Base` factory now uses `Registry::instance()` instead of a static map (removes static-init order bugs).
- Per-species intercept layout threaded through descriptor / regressor (matches CORE's blueprint update).
- C++ standard upgraded from C++14 to **C++17**.

### Fixed
- `DM_mEAD`: warns on the rotational-invariance defect for L > 0 with F ≠ F_SQ; warning text now references the new `DM_PS_*` descriptors as the recommended replacement.
- `DM_mJoin` corrections and `drfidx` added.
- Incorrect blip argument in the density formula (LaTeX docs and the corresponding code).

### Performance
- `final` qualifier on `F_RLR` / `F_SQRT` / `F_SQ` / `F_Blip` `calc_F` and `calc_dF` — guarantees devirtualisation of these hot virtual calls.

### Internal
- `dm_impl_utils.h` — shared utilities (`StackOrHeapBuffer`, `FACT_TBL`, `pow_i06`, `precompute_powers`) deduplicated out of `dm_blip.cpp` and `dm_mEAD.hpp`.
- Substantive README rewrite (Session B).
- Community / governance pointer files added.
- CI/CD: lint + release stages, VERSION wired from blueprint TOML.

## [1.2.0] and earlier (condensed)

The `v1.0.0` → `v1.1.1` releases were the "internal playground" era. MODELS held the descriptor and regressor zoo in flux: iterations focused on `FIX...` feature support, basis-function additions, and ad-hoc patches between releases. No published changelog covers this era; treat anything before `1.3.0-beta.1` as undocumented in the modern sense.
