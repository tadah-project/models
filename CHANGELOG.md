# Changelog

All notable changes to **Tadah.MODELS** (descriptors, cutoffs, basis functions, regressors) are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- `embf0_absent_consequence(bool linear)` (`embf0.h`): the text that explains what a trained DM_REAM potential without an `EMBF0` line means under the default `ZERO` and how to migrate (`EMBF0 FREE`). DM_REAM's once-per-process notice now builds its WARNING from it (output byte-identical), and the LAMMPS pair style prints the same text in its `pair_coeff` WARNING (BUG-08).

## [1.6.0] — 2026-09-25

### Added
- `require_dimer_element_pairs(context, d2, dm)` (`descriptors/dimer_element_pairs.h`) and `D_Base::applies_to_every_element_pair()`. The `DIMER` evaluators (MLIP's `DescriptorsCalc::calc_dimer`, `pair_style tadah`'s `compute_dimers`) evaluate every atom as Z=1, so a `TYPE2B`/`TYPEMB` descriptor declared for specific element pairs (`N N`) applied to no pair they see, and the model predicted from zero descriptors without a word (`tadah predict` gave 0 eV for the published Tadah!N2 potential with `N N` pairs). With `DIMER` on, each switched-on, non-empty descriptor must apply to every element pair, i.e. be declared with `* *`; otherwise the function throws `DIMER evaluates every atom as Z=1; declare TYPE2B/TYPEMB with * * element pairs` and names the line. A `D2_mJoin`/`DM_mJoin` applies to every pair when each descriptor it joins does. MLIP and Tadah.LAMMPS call it when they load a potential. Test: `test_dimer_element_pairs.cpp`.

### Changed
- Framework version `1.6.0`, read from CORE's `blueprint/blueprint_tadah.toml` after the Tadah dependencies are fetched, so the module reports the right version in every layout (sibling checkout or FetchContent). A missing blueprint is now a configure error instead of a silent `0.0.0`.
- Tests take Catch2 and its test discovery from cmake-tools (`tadah_setup_tests()`, `include(Catch)`); the vendored `tests/catch2/catch.hpp`, `tests/Catch.cmake` and `tests/CatchAddTests.cmake` are removed. Every discovered CTest entry now selects exactly its own test case and passes `--warn NoTests`, so an entry that matches no test case fails instead of passing.
- cmake-tools is fetched at `v1-dev` instead of `main`.
- Comments and a test name no longer refer to the removed BASINS HPO strategy: the `DM_mBlipEnv` meta test is `DM_mBlipEnv::meta: CGRIDMB is critical and the cross-param check fails a collapsed ρ`, and the comments describe the ρ collapse itself (a radial GTO centre outside the data's neighbour distances).
- The README badges and links and the pointers to the umbrella's community files (`CITATION.md`, `CONTRIBUTING.md`, …) point to the `v1-dev` line: GitLab `v1-dev`, `https://tadah.readthedocs.io/v1-dev/`. They pointed to `develop`.

### Performance
- Blip5 value and derivative are evaluated in one fused pass.

### Removed
- The static `D_Base::get_init_atoms(const Context &, std::string)`. Nothing called it, and it left one parameter among the element tokens and did not check the line's length. Descriptors use `D2_Base::get_init_atoms` and `DM_Base::get_init_atoms`.

### Fixed
- A `TYPE2B` or `TYPEMB` line without the parameters its descriptor takes is now an error that names the line. The key's format does not know how many parameters a descriptor takes, so two kinds of line reached the descriptor:
  - a line with fewer tokens than parameters, such as `TYPE2B D2_Blip` (the form of potential files written before descriptors named their parameters), crashed the descriptor's constructor with a vector erase past its end. It now fails with `TYPE2B D2_Blip takes 2 parameters before its element pairs; the line has 0: TYPE2B D2_Blip`;
  - a line with element pairs where the grid sizes belong, such as `TYPE2B D2_Blip * *` or `TYPEMB DM_Blip 1 N N`, was read as parameters followed by no element pair, so the descriptor applied to no pair of elements and a potential predicted wrong energies without an error. It now fails with `TYPE2B D2_Blip takes 2 parameters before its element pairs; parameter 1, the number of CGRID2B values, must be a non-negative integer, not "*": TYPE2B D2_Blip * *`.

  The check is `D_Base::element_tokens`, used by `D2_Base` and `DM_Base`: the line must have every parameter, and each size parameter (one per list key of the descriptor, the tokens `D_mJoin` reads as sizes) must be a non-negative integer. The descriptors read and convert their other parameters themselves. Test: `test_descriptor_type_line.cpp`.

## [1.5.1] — 2026-09-25

### Fixed
- **`D2_mJoin` / `DM_mJoin` no longer round their grids.** `D_mJoin::parse_context()` re-rendered every split CGRID/SGRID/CEMBFUNC/SEMBFUNC/AMPGRIDMB/KDEG\* value with NUMERIC (default 12) fixed decimals on every construction and every `reinit()`: values closer than 1e-12 collapsed, values below 5e-13 became 0, and an HPO step that small on an mJoin grid was invisible to the descriptor. The split is lossless now, so standalone mJoin numerics shift by up to 5e-13 absolute. Neither `parse_context()` nor the evidence-approximation solver (`ea.h`) inserts NUMERIC into the caller's context any more (one source of the `NUMERIC 12` line in HPO potentials). Test: `test_mjoin_precision.cpp`.
- `D_Base::get_grid` accepts subnormal grid values (e.g. `1e-310` written by an HPO step) instead of failing with "Invalid grid specification".
- **`normalize_knot_rcut` (train, predict, `refit --edit`, hpo startup) sets a knot channel's RCUT2B/RCUTMB to exactly its largest knot.** A stated value within 1e-9 relative was kept, while hpo re-derives the cutoff exactly, so the same potential could run with slightly different cutoffs in train/predict and in hpo. A gap above 1e-9 relative is still a WARNING; a smaller one (e.g. the ~3e-13 left in knot potentials written by <= 1.5.0) is fixed at INFO, except under a smooth cutoff: `pair_style tadah` reads the stated RCUT, so LAMMPS then differs from train/predict/hpo, and the WARNING says to re-write the potential (train, `refit --edit` or hpo `best_pot.tadah`). Test: `test_pot_curve_edit.cpp`.

### Changed
- **BREAKING — DM_REAM evaluates F(0) = 0 by construction (`EMBF0 ZERO`, the new default).** Every embedding basis function is zero-referenced, g~_m(rho) = g_m(rho) - g_m(0), in `calc_aed` AND `eam_embed` (the LAMMPS EAMFAST path), with g_m(0) recomputed in the constructor and in every `reinit()` from the CEMBFUNC/SEMBFUNC/KDEGEMB values actually evaluated, and an explicit rho == 0 short-circuit so an isolated atom gives exactly 0.0 regardless of FP contraction. For linear models forces and stresses are unchanged (g'_m is untouched) and every atomic energy moves by the same constant; for non-linear models (BF_Polynomial2, non-linear KRR kernels) energies and forces change. `EMBF0 FREE` keeps the raw basis of 1.5.0 (textually unchanged code path, bit-identical). A DM_REAM model built from a context without the key prints a one-time notice (WARNING for a trained potential). New `embf0.h` helpers (`stored_embf0`, `context_has_ream`, `normalise_embf0_key`, `check_embf0_precedence`, ...); `DM_Base::embf0_mode()` / `raw_embedding_f0()` (DM_mJoin forwards).
- **`tadah refit` fits what the potential evaluates.** Under `EMBF0 ZERO` the dense `EAMRefit` fits the zero-referenced embedding against the zero-referenced table F_tab(rho) - F_tab(0) (`F_ref`; the table spline is never modified), `F_fit`/`F_setfl` are zero-referenced (refit.eam starts at F = 0.0 exactly), the soft F(0) pin is not used (exact by construction), and `build_pot_context` writes the `EMBF0` line. `EMBF0 FREE` keeps the 1.5.0 fit bit-identical.
- **`refit --edit` (potedit) edits F in the potential's reference.** `CurveView` evaluates the zero-referenced F under ZERO; `remove_basis` compensates with zero-referenced functions and widens the region to `[0, largest support]` when the removed function has g(0) != 0 (reported); `rebuild_pot_context` writes the `EMBF0` line (DM_REAM only).

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
