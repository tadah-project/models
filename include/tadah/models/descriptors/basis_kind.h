#ifndef TADAH_MODELS_DESCRIPTORS_BASIS_KIND_H
#define TADAH_MODELS_DESCRIPTORS_BASIS_KIND_H

// Single source of truth for the descriptor RADIAL/EMBEDDING basis family.
//
// A descriptor's basis is one of a small, fixed set of shapes. Historically the
// choice was encoded three incompatible ways: (1) a per-class name (D2_Blip vs
// D2_BP vs D2_Blip5 vs D2_Knot5), (2) a compile-time template policy
// (radial_policies.h BlipRadial/GaussianRadial), and (3) a runtime integer code
// baked into DM_REAM's TYPEMB token with an ad-hoc private enum + scattered
// if-chains. `BasisKind` unifies the {enum <-> wire-code <-> name <-> traits}
// relationship in ONE place so every consumer (descriptors, refit, the basis
// analysis engine, HPO capability queries) shares it.
//
// The actual per-basis math lives in d_basis_functions.h; the runtime dispatch
// that maps a BasisKind to those primitives lives in basis_dispatch.h.

#include <cstdint>
#include <stdexcept>
#include <string>

namespace tadah {
namespace models {

/** \brief The radial/embedding basis family of a descriptor.
 *
 * The underlying integer values are the FROZEN wire codes written into the
 * TYPEMB/TYPE2B `<basis>` token (they must match DM_REAM potential files and
 * `tadah refit` output). Do not renumber.
 */
enum class BasisKind : std::uint8_t {
  Blip     = 0,  ///< cubic B-spline (compact support, C^2)
  Gaussian = 1,  ///< GTO exp(-eta (r-rs)^2) (NON-compact)
  Knot5    = 2,  ///< one-sided quintic/cubic truncated power (compact; eta gauge)
  Blip5    = 3,  ///< quintic M6 B-spline bump (compact, C^4)
};

/** \brief Static properties of a basis family used across the codebase. */
struct BasisTraits {
  /** Basis self-terminates (finite support) so a hard Cut_Dummy cutoff is
   *  sound: the density/pair function vanishes smoothly at the cutoff BY
   *  CONSTRUCTION. False for Gaussian (only a real, smooth cutoff avoids a
   *  force discontinuity — the "no-cutoff smell"). */
  bool compact_support;
  /** The width channel (SGRIDMB / SEMBFUNC) is a pure gauge and is ABSENT
   *  from the token stream (knot family): (c, eta) -> (c s^k, eta/s) leaves
   *  the function invariant, so fitters pin eta = 1. */
  bool eta_is_gauge;
  /** Human-readable name (matches the canonical `tadah refit` REFIT_BASIS). */
  const char *name;
};

/** \brief Compile-time trait lookup (no runtime cost). */
constexpr BasisTraits basis_traits(BasisKind k) {
  switch (k) {
    case BasisKind::Blip:     return {true,  false, "Blip"};
    case BasisKind::Gaussian: return {false, false, "Gaussian"};
    case BasisKind::Knot5:    return {true,  true,  "Knot5"};
    case BasisKind::Blip5:    return {true,  false, "Blip5"};
  }
  return {true, false, "Blip"};  // unreachable (enum is exhaustive above)
}

inline int          basis_code(BasisKind k) { return static_cast<int>(k); }
inline const char  *basis_name(BasisKind k) { return basis_traits(k).name; }

/** \brief Map the FROZEN TYPExB integer `<basis>` code to a BasisKind.
 *
 * Throws on an unknown code rather than silently degrading to Blip — a pot
 * written by a newer Tadah must not be evaluated with the wrong basis without
 * any warning. Replaces DM_REAM's private `basis_from_bcode`.
 */
inline BasisKind basis_kind_from_code(int code) {
  switch (code) {
    case 0: return BasisKind::Blip;
    case 1: return BasisKind::Gaussian;
    case 2: return BasisKind::Knot5;
    case 3: return BasisKind::Blip5;
    default:
      throw std::runtime_error(
          "BasisKind: unknown basis code " + std::to_string(code) +
          " (known: 0=Blip, 1=Gaussian, 2=Knot5, 3=Blip5).");
  }
}

/** \brief Map a user-facing basis string (REFIT_BASIS and its aliases) to a
 *  BasisKind. Lenient by design (unknown -> Blip, matching the historical
 *  `tadah refit` default); the exact wire contract is `basis_kind_from_code`.
 *  Hoisted verbatim from the former EAMRefit ctor string ladder. Note the
 *  Knot3 / KnotMix / Mixed aliases all map to the Knot5 family — the per-atom
 *  KDEG* degree arrays (3|5) then select cubic vs quintic knot atoms. */
inline BasisKind basis_kind_from_string(const std::string &s) {
  if (s == "Gaussian" || s == "gaussian" || s == "BP" || s == "D2_BP")
    return BasisKind::Gaussian;
  if (s == "Knot5" || s == "KNOT5" || s == "knot5" ||
      s == "Knot3" || s == "KNOT3" || s == "knot3" ||
      s == "KnotMix" || s == "KNOTMIX" || s == "knotmix" ||
      s == "Mixed" || s == "MIXED")
    return BasisKind::Knot5;
  if (s == "Blip5" || s == "BLIP5" || s == "blip5" || s == "M6")
    return BasisKind::Blip5;
  return BasisKind::Blip;
}

}  // namespace models
}  // namespace tadah

#endif  // TADAH_MODELS_DESCRIPTORS_BASIS_KIND_H
