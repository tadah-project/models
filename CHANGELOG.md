# Changelog

All notable changes to **Tadah.MODELS** (descriptors, cutoffs, basis functions, regressors) are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Fixed
- **`D2_mJoin` / `DM_mJoin` no longer round their grids.** `D_mJoin::parse_context()` re-rendered every split CGRID/SGRID/CEMBFUNC/SEMBFUNC/AMPGRIDMB/KDEG\* value with NUMERIC (default 12) fixed decimals on every construction and every `reinit()`: values closer than 1e-12 collapsed, values below 5e-13 became 0, and an HPO step that small on an mJoin grid was invisible to the descriptor. The split is lossless now, so standalone mJoin numerics shift by up to 5e-13 absolute. Neither `parse_context()` nor the evidence-approximation solver (`ea.h`) inserts NUMERIC into the caller's context any more (one source of the `NUMERIC 12` line in HPO potentials). Test: `test_mjoin_precision.cpp`.
- `D_Base::get_grid` accepts subnormal grid values (e.g. `1e-310` written by an HPO step) instead of failing with "Invalid grid specification".
- **`normalize_knot_rcut` (train, predict, `refit --edit`, hpo startup) sets a knot channel's RCUT2B/RCUTMB to exactly its largest knot.** A stated value within 1e-9 relative was kept, while hpo re-derives the cutoff exactly, so the same potential could run with slightly different cutoffs in train/predict and in hpo. A gap above 1e-9 relative is still a WARNING; a smaller one (e.g. the ~3e-13 left in knot potentials written by <= 1.5.0) is fixed at INFO, except under a smooth cutoff: `pair_style tadah` reads the stated RCUT, so LAMMPS then differs from train/predict/hpo, and the WARNING says to re-write the potential (train, `refit --edit` or hpo `best_pot.tadah`). Test: `test_pot_curve_edit.cpp`.

## [1.4.0] — 2026-07-14

### Added
- **Bump-support invariant plumbing** (`c + w/eta <= RCUT` for the compact bump families under a hard Cut_Dummy cutoff — a truncated overhang is an energy/force step at the cutoff): `BasisTraits::bump_halfwidth` (2 Blip / 3 Blip5) + `bump_halfwidth(BasisKind)` + the shared gate tolerance `bump_rcut_tol` (1e-6 relative, so hpotarget-quantized bounds sitting exactly on the cap still pass) in `basis_kind.h`; `potedit::bump_channel_halfwidth(ctx, suffix)` (Context-level twin of `derived_knot_rcut` for bumps) and `potedit::channel_cutoff_is_discontinuous(ctx, suffix)` (RCTYPE\* absent = Cut_Dummy, matching the descriptor ctor default) in `pot_curve_edit`.
- **Read-path warn parity for bumps**: `D_Base::warn_if_bump_support_beyond_rcut` (warn-once, discontinuous cutoffs only), called from the ctor + `reinit()` of `D2_Blip`, `D2_Blip5` and `DM_REAM` (density grids) — the bump analogue of the knot families' `warn_if_knots_beyond_rcut`. Reading a legacy spilled pot stays warn-only: its truncation is part of the fitted model.
- EAM-collapse capability API on `DM_Base`: `eam_form()`, `eam_density()` (scalar pair density and r-derivative, cutoff- and WATOM-weighted), `eam_embed()` (embedding value/derivative collapsed against the full model weight vector, block offsets handled internally). Implemented by `DM_REAM`, whose many-body energy is exactly `F(rho_bar)`. Lets MD drivers (pair_style tadah's `EAMFAST` fast path) collapse linear EAM-form models to classic scalar pair/density/embedding functions. Default: not EAM-form (no behavioural change for other descriptors).

### Changed
- **`potedit::add_basis`: a bump support CROSSING rcut is now a hard error under a hard (Cut_Dummy) cutoff** (was a warning; the message gives the minimal `--eta` that fits). Smooth cutoffs keep the warning (the overhang is damped to zero at rcut). New defaulted parameter `hard_cutoff` (callers with a pot Context pass `channel_cutoff_is_discontinuous`).
- **Dense-mode `EAMRefit::make_grid` / `make_grid_windowed`: radial bump grids end their OUTERMOST support exactly at rcut** (spacing `(hi-lo)/(n-1+w)`, eta `1/dx`; previously the outermost bump sat AT rcut with half its support truncated — every dense bump pot violated the bump-support gate its own hpotarget is now checked against). The embedding grid (no r cutoff) is unchanged; numerics of dense bump refits shift accordingly.

## [1.3.0-beta.1] — 2026-05-20

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
- **Per-descriptor `meta()` schema for the HPO bound-validation pipeline**: `param_meta.h` declares `ParamMeta` / `DescriptorMeta` / `SeedView` / `HealthVerdict`, `closed_form_meta.h` provides the shared 2-body factory, and `D_Base::meta()` is the new virtual hook (default `nullopt`). `meta()` implementations land for D2_LJ/MIE/ZBL/Blip/BP/EAM and DM_Blip/EAD/EAM/mEAD/Poly. The `DM_mEAD<F_mEnv<Cut_SinInv>>` specialisation (registered as `DM_mBlipEnv`) carries the ρ-collapse cross-param check that would have caught the BASINS silent zero-block failure.

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
- Unused `suggested_value()` descriptor hook dropped from the descriptor interface (no external callers).
- Test build routed through `cmake-tools/cmake/TadahTesting.cmake` helpers (`tadah_declare_testing` / `tadah_add_tests`).

## [1.2.0] and earlier (condensed)

The `v1.0.0` → `v1.1.1` releases were the "internal playground" era. MODELS held the descriptor and regressor zoo in flux: iterations focused on `FIX...` feature support, basis-function additions, and ad-hoc patches between releases. No published changelog covers this era; treat anything before `1.3.0-beta.1` as undocumented in the modern sense.
