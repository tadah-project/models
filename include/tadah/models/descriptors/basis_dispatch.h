#ifndef TADAH_MODELS_DESCRIPTORS_BASIS_DISPATCH_H
#define TADAH_MODELS_DESCRIPTORS_BASIS_DISPATCH_H

// The ONE runtime dispatch that maps a BasisKind to the hand-tuned basis
// primitives in d_basis_functions.h. Every runtime "which basis?" decision in
// the codebase funnels through here instead of a per-site if/switch on a
// descriptor name or a private enum (the anti-pattern the design removes).
//
// Cost note (goal d, "never regress speed"): each function is a single
// `static inline` switch over a 4-value enum -> a jump table. It replaces the
// former per-call if-ladders (up to 3 branch compares) and forwards to the
// exact same inlined primitive, so codegen is neutral-to-better versus the
// hand-written dispatch it subsumes. The basis is validated once at descriptor
// construction (basis_kind_from_code throws on bad codes), so the trailing
// unreachable/throw is a cold, never-taken path.
//
// Two orientations, matching the descriptor contract:
//   * radial_*   : LEFT-active knots, cutoff-weighted (density / two-body).
//   * embedding_*: RIGHT-active knots, NO cutoff (embedding F(rho)).
//   * *_bare     : no cutoff weight (basis analysis / plotting).

#include <string>

#include <tadah/models/descriptors/basis_kind.h>
#include <tadah/models/descriptors/d_basis_functions.h>

namespace tadah {
namespace models {
namespace basis_detail {

[[noreturn]] inline void unknown_basis(BasisKind k) {
  throw std::runtime_error(
      "basis_dispatch: unhandled BasisKind code " +
      std::to_string(static_cast<int>(k)));
}

}  // namespace basis_detail

// --------------------------------------------------------------------------
//  Radial basis  psi(r) * f_c(r)   (LEFT-active knots)
//  s = width (eta / SGRID; ignored by the knot family), c = centre (mu / CGRID),
//  deg = knot degree 3|5 (ignored by non-knot families).
// --------------------------------------------------------------------------

/** \brief psi_c(r) * f_c(r) for the given basis. */
inline double radial_value(BasisKind k, double r, double fc,
                           double s, double c, int deg = 5) {
  switch (k) {
    case BasisKind::Blip:     return B (s * (r - c), fc);
    case BasisKind::Blip5:    return B5(s * (r - c), fc);
    case BasisKind::Knot5:    return deg == 3 ? K3L(r - c, fc) : K5L(r - c, fc);
    case BasisKind::Gaussian: return G (r, s, c, fc);
  }
  basis_detail::unknown_basis(k);
}

/** \brief d/dr [ psi_c(r) * f_c(r) ] for the given basis (fcp = f_c'(r)). */
inline double radial_deriv(BasisKind k, double r, double fc, double fcp,
                           double s, double c, int deg = 5) {
  switch (k) {
    case BasisKind::Blip:     return dB (s * (r - c), s, fc, fcp);
    case BasisKind::Blip5:    return dB5(s * (r - c), s, fc, fcp);
    case BasisKind::Knot5:    return deg == 3 ? dK3L(r - c, 1.0, fc, fcp)
                                              : dK5L(r - c, 1.0, fc, fcp);
    case BasisKind::Gaussian: return dG (r, s, c, fc, fcp);
  }
  basis_detail::unknown_basis(k);
}

// --------------------------------------------------------------------------
//  Embedding basis  g_m(rho)   (RIGHT-active knots, NO cutoff)
// --------------------------------------------------------------------------

/** \brief g_m(rho) for the given basis (embedding, no cutoff). */
inline double embedding_value(BasisKind k, double rho,
                              double s, double c, int deg = 5) {
  switch (k) {
    case BasisKind::Blip:     return B (s * (rho - c));
    case BasisKind::Blip5:    return B5(s * (rho - c));
    case BasisKind::Knot5:    return deg == 3 ? K3R(rho - c) : K5R(rho - c);
    case BasisKind::Gaussian: return G (rho, s, c);
  }
  basis_detail::unknown_basis(k);
}

/** \brief g'_m(rho) w.r.t. rho for the given basis (embedding, no cutoff). */
inline double embedding_deriv(BasisKind k, double rho,
                              double s, double c, int deg = 5) {
  switch (k) {
    case BasisKind::Blip:     return dB (s * (rho - c), s);
    case BasisKind::Blip5:    return dB5(s * (rho - c), s);
    case BasisKind::Knot5:    return deg == 3 ? dK3R(rho - c, 1.0)
                                              : dK5R(rho - c, 1.0);
    case BasisKind::Gaussian: return dG (rho, s, c);
  }
  basis_detail::unknown_basis(k);
}

// --------------------------------------------------------------------------
//  Bare radial basis  psi_c(r)  (no cutoff weight) — basis analysis / plots.
//  Matches the legacy `double(*)(double,double,double)` shape of models::B /
//  models::G etc. exactly, so it is a drop-in for the analysis engine.
// --------------------------------------------------------------------------

/** \brief psi_c(r) with no cutoff weight (LEFT-active knots). */
inline double radial_value_bare(BasisKind k, double r,
                                double s, double c, int deg = 5) {
  switch (k) {
    case BasisKind::Blip:     return B (s * (r - c));   // == B(r, s, c)
    case BasisKind::Blip5:    return B5(s * (r - c));
    case BasisKind::Knot5:    return deg == 3 ? K3L(r - c) : K5L(r - c);
    case BasisKind::Gaussian: return G (r, s, c);
  }
  basis_detail::unknown_basis(k);
}

/** \brief d/dr psi_c(r) with no cutoff weight (LEFT-active knots). */
inline double radial_deriv_bare(BasisKind k, double r,
                                double s, double c, int deg = 5) {
  switch (k) {
    case BasisKind::Blip:     return dB (s * (r - c), s);  // == dB(r, s, c)
    case BasisKind::Blip5:    return dB5(s * (r - c), s);
    case BasisKind::Knot5:    return deg == 3 ? dK3L(r - c, 1.0)
                                              : dK5L(r - c, 1.0);
    case BasisKind::Gaussian: return dG (r, s, c);
  }
  basis_detail::unknown_basis(k);
}

}  // namespace models
}  // namespace tadah

#endif  // TADAH_MODELS_DESCRIPTORS_BASIS_DISPATCH_H
