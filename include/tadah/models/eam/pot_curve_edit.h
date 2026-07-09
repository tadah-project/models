#ifndef TADAH_MODELS_EAM_POT_CURVE_EDIT_H
#define TADAH_MODELS_EAM_POT_CURVE_EDIT_H

// Editing primitives for EXISTING `tadah refit` potentials (pot.tadah with a
// D2_Knot5 / D2_Blip / D2_Blip5 two-body descriptor and a DM_REAM many-body
// descriptor): parse the three EAM curves into per-channel views, add or
// remove ONE basis function, and rebuild a consistent potential Context
// (grids, TYPE2B/TYPEMB count tokens, WEIGHTS/FIXWEIGHT/FIXINDEX, RCUT).
//
// Deliberately self-contained over the shared basis primitives
// (d_basis_functions.h via basis_dispatch.h) so it works in plain MLIP
// builds — it does NOT reuse HPO's Ceres-gated BlipFit (which also carries an
// incompatible BasisKind enum).
//
// RCUT contract (family-aware; the load-bearing asymmetry):
// - Knot channels: the basis vanishes C^4-smoothly at its own largest knot,
//   so the TRUE effective cutoff IS max(knot centre) and RCUT2B/RCUTMB are
//   DERIVED quantities — `derived_knot_rcut` computes them and
//   `normalize_knot_rcut` warns + rewrites a stale stated value. Moving RCUT
//   anywhere >= max knot is prediction-exact.
// - Blip channels: the hard Cut_Dummy truncation at RCUT was part of the fit
//   (the basis overhangs the cutoff by design), so RCUT is AUTHORITATIVE and
//   is never auto-changed — only warned about.

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <tadah/core/context.h>
#include <tadah/models/descriptors/basis_kind.h>

namespace tadah {
namespace models {
namespace potedit {

enum class Channel { Phi, Rho, F };

/** \brief Inner sampling bound of the removal-compensation window as a
 *  fraction of rcut. The pot file does not carry the fit's own r_lo, so this
 *  mirrors the inspection domain used by the refit curve dumps. Heuristic —
 *  documented in refit.rst. */
constexpr double kCompensationRLoFrac = 0.15;

/** \brief One parsed EAM curve: basis family + parallel per-BF arrays.
 *
 * Invariants: mu/coeff always sized n(); eta sized n() for bump families and
 * EMPTY for knots (eta = 1 gauge); deg EMPTY (= all quintic) or sized n().
 * coeff maps to the WEIGHTS 2b block (Phi), AMPGRIDMB (Rho) or the WEIGHTS
 * embedding block (F). rcut > 0 only for the radial channels.
 */
struct CurveView {
  BasisKind family = BasisKind::Blip;
  bool right_active = false;  ///< true only for F with the knot family
  std::vector<double> mu;
  std::vector<double> eta;    ///< empty for the knot family
  std::vector<double> coeff;
  std::vector<int> deg;       ///< empty = all quintic
  double rcut = 0.0;          ///< RCUT2B/RCUTMB; 0 for F (no r cutoff)

  std::size_t n() const { return mu.size(); }
  int degree(std::size_t k) const {
    return deg.empty() ? 5 : deg[k];
  }
  double width(std::size_t k) const {  ///< eta of BF k (1 for knots)
    return eta.empty() ? 1.0 : eta[k];
  }

  /// Curve value/derivatives at x: sum_k coeff_k g_k(x). Radial channels
  /// apply the hard Cut_Dummy semantics (0 beyond rcut), matching runtime.
  double value(double x) const;
  double deriv(double x) const;
  double deriv2(double x) const;

  /// Active x-interval of BF k. One-sided knots return a half-line encoded
  /// with +-support_infinity() on the open side.
  std::pair<double, double> support(std::size_t k) const;
  static double support_infinity();
};

/** \brief A refit potential parsed into its three curves + retained Context. */
struct PotEAMView {
  core::Context ctx;  ///< the parsed original; untouched keys preserved
  CurveView phi, rho, F;
  std::string element;
};

/** \brief Parse a refit-shaped potential Context into per-channel views.
 *  Throws std::runtime_error with an actionable message on any other shape
 *  (D2_mJoin, Gaussian/D2_BP, non-DM_REAM many-body, inconsistent arrays). */
PotEAMView parse_pot_view(const core::Context &ctx);

/** \brief Rebuild a full potential Context from an (edited) view: grids,
 *  KDEG / SGRID / AMPGRIDMB arrays, TYPE2B/TYPEMB count tokens, WEIGHTS =
 *  [phi|F], FIXINDEX "1-K" + FIXWEIGHT, and RCUT2B/RCUTMB (derived for knot
 *  channels, preserved for bump channels). Self-checks by constructing the
 *  D2/DM descriptors from the result before returning. */
core::Context rebuild_pot_context(const PotEAMView &v);

/** \brief Effective cutoff of a knot channel, from the pot Context alone.
 *  suffix is "2B" (TYPE2B family + CGRID2B) or "MB" (DM_REAM bcode +
 *  CGRIDMB). Returns max knot centre for knot channels; std::nullopt for
 *  bump families or when the channel keys are absent/non-refit. */
std::optional<double> derived_knot_rcut(const core::Context &ctx,
                                        const std::string &suffix);

/** \brief F2 normalisation shared by every consumer (refit --edit, train,
 *  predict, hpo startup): for each KNOT channel compare the stated
 *  RCUT2B/RCUTMB against `derived_knot_rcut` and, on mismatch, log a WARNING
 *  (prefixed with `who`) and rewrite the key in ctx (RCUT*MAX follows via
 *  Context internal-key derivation). Bump channels are never modified.
 *  Returns true if anything changed. */
bool normalize_knot_rcut(core::Context &ctx, const std::string &who);

struct AddSpec {
  double centre = 0.0;
  std::optional<double> eta;  ///< bump families only
  std::optional<int> deg;     ///< knot family only: 3|5 (default 5)
};

struct AddReport {
  std::size_t index = 0;   ///< position of the new BF in the arrays (appended)
  double eta_used = 1.0;
  int deg_used = 5;
  std::vector<std::string> warnings;
};

/** \brief Insert one basis function with coefficient 0 (curve, derivatives
 *  and extrapolation bit-identical until the coefficient is varied).
 *  `radial` selects the validation rules (centre > 0, rcut interplay).
 *  Throws on invalid spec (eta on knots, deg on bumps, bad centre, bump
 *  support entirely beyond rcut). */
AddReport add_basis(CurveView &c, const AddSpec &s, bool radial);

struct RemoveReport {
  std::size_t removed_idx = 0;
  double removed_mu = 0.0;
  double removed_coeff = 0.0;
  std::vector<std::size_t> adjusted;  ///< compensator indices (post-removal order NOT applied; indices refer to the pre-removal arrays)
  std::vector<double> delta_c;
  double resid_plain = 0.0;  ///< stacked value/d1/d2 residual, no compensation
  double resid_comp = 0.0;   ///< same after compensation (<= resid_plain)
  double vresid_plain = 0.0; ///< value-only RMS over the affected region
  double vresid_comp = 0.0;
  bool compensated = false;
  std::vector<std::string> warnings;
};

/** \brief Remove the basis function whose centre is nearest centre_request
 *  (ties -> lower index). Unless its coefficient is exactly 0, the (up to)
 *  two remaining nearest-in-centre coefficients are least-squares adjusted
 *  against the removed contribution's value + 1st + 2nd derivative on a
 *  dense sample of the affected region (per-order rows normalised by the
 *  removed term's RMS so the orders are balanced). Throws when the channel
 *  has a single basis function. */
RemoveReport remove_basis(CurveView &c, double centre_request, bool radial);

}  // namespace potedit
}  // namespace models
}  // namespace tadah

#endif  // TADAH_MODELS_EAM_POT_CURVE_EDIT_H
