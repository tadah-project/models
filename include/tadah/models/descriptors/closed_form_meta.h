#ifndef TADAH_MODELS_DESCRIPTORS_CLOSED_FORM_META_H
#define TADAH_MODELS_DESCRIPTORS_CLOSED_FORM_META_H

// Helper for descriptors with no learnable continuous parameters of
// their own (D2_LJ, D2_MIE, D2_ZBL, D2_EAM, DM_EAM). Each has only
// RCUT2B (or RCUTMB) under the user's control; the audit watches for
// the cutoff being too small to find any neighbours. D2_EAM and DM_EAM
// take the cutoff of their setfl file (they refuse another), so for them
// it is not a search target.

#include <optional>
#include <string>
#include <utility>

#include <tadah/models/descriptors/param_meta.h>
#include <tadah/core/stats/data_stats.h>

namespace tadah {
namespace models {

/** \brief Minimal `DescriptorMeta` for closed-form descriptors that
 *  expose only a real-space cutoff. `cutoff_key` is one of
 *  `"RCUT2B"` / `"RCUTMB"`. The descriptor `label` is taken at
 *  runtime from the caller.
 *
 *  `cutoff_from_setfl`: the descriptor takes the cutoff of its setfl file
 *  (D2_EAM, DM_EAM; D_Base::cutoff_fixed_by()). The cutoff is then
 *  Optional and has no bound suggester, so neither the critical-OPTIM
 *  audit nor `tadah analysis suggest_bounds` proposes searching it; the
 *  health check (a cutoff that finds no neighbours) stays. Otherwise it is
 *  Critical, with a suggester. */
inline DescriptorMeta make_closed_form_meta(
    const std::string& descriptor_label,
    const std::string& cutoff_key,
    bool cutoff_from_setfl = false) {
  using DataStats = tadah::core::stats::DataStats;
  static constexpr std::size_t I_P50 = tadah::core::stats::IDX_P50;
  static constexpr std::size_t I_P99 = tadah::core::stats::IDX_P99;

  DescriptorMeta dm;
  dm.label = descriptor_label;

  ParamMeta p_rcut;
  p_rcut.key      = cutoff_key;
  p_rcut.arg_pos  = -1;
  p_rcut.scale    = Scale::Distance;
  p_rcut.units    = "Å";
  p_rcut.role     = cutoff_from_setfl
      ? "Real-space cutoff: the setfl file's (SETFL), not a search target."
      : "Real-space cutoff (closed-form descriptor: only knob).";
  p_rcut.hard.lo  = 0.0;
  p_rcut.crit     = cutoff_from_setfl ? Criticality::Optional
                                      : Criticality::Critical;
  p_rcut.scope    = PairScope::Global;

  // Capture cutoff_key by value so the closure doesn't dangle on the
  // caller's string.
  const std::string ck = cutoff_key;

  if (!cutoff_from_setfl)
    p_rcut.suggest = [](const DataStats& d, const SeedView&, int)
        -> std::pair<double, double> {
      if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
      const auto& pp = d.per_pair[0];
      double low = pp.r_percentiles[I_P50];
      if (pp.r_peaks.size() >= 2) low = pp.r_peaks[1].position * 1.05;
      return std::make_pair(low, pp.r_percentiles[I_P99]);
    };
  p_rcut.check = [ck](const DataStats& d, const SeedView& s) -> HealthVerdict {
    if (!s.contains(ck)) return HealthOk{};
    const double rc = s.get<double>(ck);
    if (d.per_pair.empty() && d.n_atoms > 0)
      return HealthFail{ck + "=" + std::to_string(rc) +
                        " is too small — no neighbours found."};
    if (d.per_pair.empty()) return HealthOk{};
    if (rc < d.per_pair[0].r_min)
      return HealthFail{ck + "=" + std::to_string(rc) +
                        " is below dataset r_min=" +
                        std::to_string(d.per_pair[0].r_min) + "."};
    return HealthOk{};
  };
  dm.params.push_back(std::move(p_rcut));

  dm.expected_columns = [](const SeedView&) -> std::size_t { return 1; };
  return dm;
}

}  // namespace models
}  // namespace tadah

#endif  // TADAH_MODELS_DESCRIPTORS_CLOSED_FORM_META_H
