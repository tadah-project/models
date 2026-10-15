#ifndef TADAH_MODELS_DESCRIPTORS_PARAM_META_H
#define TADAH_MODELS_DESCRIPTORS_PARAM_META_H

// Per-descriptor parameter metadata schema.
//
// Each registered descriptor (D2_*, DM_*) optionally returns a
// DescriptorMeta from its overridden `D_Base::meta() const`. The
// HPO bound-validation pipeline reads this metadata to:
//
//   * Suggest data-driven OPTIM bounds (BoundSuggester).
//   * Diagnose degenerate parameter combinations (HealthCheck).
//   * Flag missing OPTIM blocks on Critical hyperparameters.
//   * Predict design-matrix column counts before construction.
//
// Default-`nullopt` semantics keep every existing descriptor working
// unchanged.

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <tadah/core/context.h>
#include <tadah/models/descriptors/dm/embed_spec.h>
#include <tadah/core/stats/data_stats.h>

namespace tadah {
namespace models {

/** \brief Physical scale category of a parameter — drives which dataset
 *  statistic the suggester anchors to.
 */
enum class Scale {
  Distance,           // Å (e.g. RCUT, CGRIDxB)
  InverseDistanceSq,  // Å^-2 (e.g. SGRIDxB — Gaussian inverse width)
  InverseDistance,    // Å^-1 (compact bump eta)
  Density,            // ρ (rho-units; e.g. CEMBFUNC, RCUTENV)
  InverseDensity,     // ρ^-1 (e.g. SEMBFUNC)
  InverseDensitySq,   // ρ^-2 (Gaussian embedding eta)
  Dimensionless,      // unitless (e.g. exponents)
  Energy,             // eV
  Force,              // eV/Å
  Count               // dimensionless integer count
};

/** \brief Hard physical limits — used by the audit only to refuse
 *  obviously impossible values (negative widths, etc.). NEVER tuning
 *  thresholds.
 */
struct HardLimit {
  std::optional<double> lo{};
  std::optional<double> hi{};
};

/** \brief Whether a parameter must be optimised.
 *  - Critical: a fixed value that misses the data can make the
 *    descriptor vanish (e.g. a CGRIDMB centre outside the neighbour
 *    distances collapses ρ); used by the pre-flight audit to warn when
 *    no OPTIM block is present.
 *  - Recommended: optimisation usually helps; not strictly necessary.
 *  - Optional: routine to fix.
 */
enum class Criticality { Critical, Recommended, Optional };

/** \brief Whether a parameter is shared across all element pairs or
 *  declared per (element-i, element-j) pair.
 */
enum class PairScope { Global, PerElementPair };

// ---------------------------------------------------------------------
//  HealthVerdict — Ok / Warn / Fail
// ---------------------------------------------------------------------

struct HealthOk {};

struct HealthWarn {
  std::string message;
};

struct HealthFail {
  std::string message;
};

using HealthVerdict = std::variant<HealthOk, HealthWarn, HealthFail>;

inline bool is_ok  (const HealthVerdict& v) { return std::holds_alternative<HealthOk>(v); }
inline bool is_warn(const HealthVerdict& v) { return std::holds_alternative<HealthWarn>(v); }
inline bool is_fail(const HealthVerdict& v) { return std::holds_alternative<HealthFail>(v); }

inline std::string verdict_message(const HealthVerdict& v) {
  if (auto* w = std::get_if<HealthWarn>(&v)) return w->message;
  if (auto* f = std::get_if<HealthFail>(&v)) return f->message;
  return {};
}

// ---------------------------------------------------------------------
//  SeedView — read-only view over the seed Context
// ---------------------------------------------------------------------

/** \brief Thin read-only wrapper over `tadah::core::Context` giving
 *  metadata closures access to the current parameter values without
 *  allowing mutation.
 *
 *  Constructed implicitly from a `const Context&`; the underlying
 *  Context must outlive the SeedView.
 */
class SeedView {
 public:
  explicit SeedView(const tadah::core::Context& c) : ctx_(&c) {}

  template <typename T>
  T get(const std::string& key) const {
    if(key.find('.')!=std::string::npos) return get<T>(key,0);
    return ctx_->get<T>(key);
  }

  template <typename T>
  T get(const std::string& key, std::size_t arg_pos) const {
    if(key.find('.')!=std::string::npos) {
      double v=embed_binding(*ctx_,key).values.at(arg_pos);
      if constexpr(std::is_arithmetic_v<T>) return static_cast<T>(v);
      else if constexpr(std::is_same_v<T,std::string>) return core::format_roundtrip(v);
    }
    return ctx_->get<T>(key, arg_pos);
  }

  std::size_t size(const std::string& key) const {
    if(key.find('.')!=std::string::npos) return embed_binding(*ctx_,key).values.size();
    return ctx_->size(key);
  }

  bool contains(const std::string& key) const {
    if(key.find('.')!=std::string::npos) {
      for(const auto &b:embed_bindings(*ctx_)) if(b.path==core::to_upper(key)) return true;
      return false;
    }
    return ctx_->contains(key);
  }

  /** Read a v_type / vector<double> in one shot. */
  std::vector<double> vec(const std::string& key) const {
    std::vector<double> out;
    if (!contains(key)) return out;
    const std::size_t n = size(key);
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      out.push_back(get<double>(key, i));
    }
    return out;
  }

  const tadah::core::Context& context() const { return *ctx_; }

 private:
  const tadah::core::Context* ctx_;
};

// ---------------------------------------------------------------------
//  Closures
// ---------------------------------------------------------------------

/** \brief Suggest a `(lo, hi)` bound from data + seed.
 *
 *  arg_pos: -1 for scalar params; 0..n-1 for vector entries
 *  (the framework calls the suggester per vector entry).
 *
 *  Closures are written to read DataStats — never to embed magic
 *  numbers. The only constants permitted in a suggester body are
 *  small structural ones (e.g. "double the FWHM") that are documented
 *  at the declaration site.
 */
using BoundSuggester = std::function<std::pair<double, double>(
    const tadah::core::stats::DataStats&,
    const SeedView&,
    int arg_pos)>;

/** \brief Validate (or warn about) a parameter or parameter combination
 *  given the seed and the data.
 */
using HealthCheck = std::function<HealthVerdict(
    const tadah::core::stats::DataStats&,
    const SeedView&)>;

/** \brief Predict the number of design-matrix columns this descriptor
 *  occupies, given a seed. Used by ColumnLayout to populate
 *  the per-block layout without instantiating descriptors.
 */
using ColumnCounter = std::function<std::size_t(const SeedView&)>;

// ---------------------------------------------------------------------
//  ParamMeta / DescriptorMeta
// ---------------------------------------------------------------------

struct ParamMeta {
  std::string  key{};
  /** Position within the key's value list (1-based, matching OPTIM (n)
   *  syntax used in HPO target files). Use -1 for scalar parameters
   *  and 0 for "vector entry whose index is filled by the framework
   *  when iterating". */
  int          arg_pos{-1};

  Scale        scale{Scale::Dimensionless};
  std::string  units{};
  std::string  role{};

  HardLimit    hard{};
  Criticality  crit{Criticality::Optional};
  PairScope    scope{PairScope::Global};

  BoundSuggester suggest{};
  HealthCheck    check{};

  /** Per-parameter aggregate thresholds used by the per-iter health
   *  monitor. All have rationale-bearing defaults; descriptors
   *  may override.
   *
   *  drift_to_bound_threshold: a parameter that sits within
   *      (1 - this) * (bound width) of either edge for >= window evals
   *      is flagged. Default 0.20 — i.e. "stuck within 20% of a bound".
   *      Rationale: 1/5 of the bound is the smallest "near-edge"
   *      window where local-search noise alone is unlikely to keep a
   *      parameter pinned across hundreds of evaluations.
   *
   *  weight_collinearity_ratio: max(|w|)/median(|w|) > this in this
   *      block triggers a warning. Default 1e3 — "1000x stretch
   *      between the largest and the typical weight is the canonical
   *      signature of a near-rank-deficient block."
   *
   *  weight_extrapolation_max: max(|w|) > this triggers a warning.
   *      Default 1e4 — extrapolation outside the training pool's
   *      energy scale is unreliable beyond this magnitude.
   *
   *  All three are *defaults* — the framework never compares against a
   *  literal in the audit code.
   */
  double drift_to_bound_threshold{0.20};
  double weight_collinearity_ratio{1.0e3};
  double weight_extrapolation_max{1.0e4};
};

struct DescriptorMeta {
  std::string                    label{};
  std::vector<ParamMeta>         params{};
  /** Whole-descriptor checks (e.g. "ρ collapse" — uses several
   *  parameters in concert). */
  std::vector<HealthCheck>       cross_param_checks{};
  /** Predict column count in the design matrix from the seed values. */
  ColumnCounter                  expected_columns{};
  /** Block-dead column-norm-ratio threshold, in [0, 1). A column block
   *  whose median column norm is < `block_dead_norm_ratio * max_block_norm`
   *  is flagged "dead" by the pre-flight audit. The default: a column
   *  block 6 orders of magnitude smaller than another block in the same
   *  matrix has no meaningful signal. */
  double block_dead_norm_ratio{1.0e-6};
};

}  // namespace models
}  // namespace tadah

#endif  // TADAH_MODELS_DESCRIPTORS_PARAM_META_H
