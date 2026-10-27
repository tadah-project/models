#include <tadah/models/descriptors/dm/rhoform.h>
#include <tadah/models/descriptors/d_mjoin.h>
#include <tadah/models/descriptors/dm/embf0.h>
#include <tadah/models/descriptors/dm/ream_density.h>
#include <tadah/models/eam/pot_curve_edit.h>

#include <algorithm>
#include <cctype>
#include <iterator>

namespace tadah {
namespace models {

namespace {
std::string upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char ch) { return std::toupper(ch); });
  return s;
}
// The cutoffs measured non-negative and non-increasing on [0, rcut]
// (rhoform_cutoff_whitelisted, rhoform_cutoff_whitelist_text).
const char *const kMonotoneCutoffs[] = {"Cut_Dummy", "Cut_Cos",  "Cut_CosS",
                                        "Cut_Tanh",  "Cut_PT",   "Cut_Poly",
                                        "Cut_PolyS"};
}  // namespace

const char *rhoform_name(RhoForm f) {
  switch (f) {
    case RhoForm::Signed:   return "SIGNED";
    case RhoForm::NonNeg:   return "NONNEG";
    case RhoForm::Monotone: return "MONOTONE";
  }
  return "SIGNED";  // unreachable (enum is exhaustive above)
}

bool rhoform_is_positive(RhoForm f) { return f != RhoForm::Signed; }

RhoForm parse_rhoform(const std::string &value) {
  const std::string v = upper(value);
  if (v == "SIGNED") return RhoForm::Signed;
  if (v == "NONNEG") return RhoForm::NonNeg;
  if (v == "MONOTONE") return RhoForm::Monotone;
  throw DensityFormError("RHOFORM must be SIGNED, NONNEG or MONOTONE (got '" +
                         value + "').");
}

namespace {
// The three forms, shared by the absent-key and the invalid-value messages.
const char *const kFormList =
    "SIGNED (rho = W f_c sum a_n psi_n with signed a_n, the density before "
    "tadah 1.6.0; add 'RHOFORM SIGNED' to use an existing potential "
    "unchanged), "
    "NONNEG (the same with every a_n >= 0: rho >= 0 by construction) or "
    "MONOTONE (rho = W f_c sum a_n S_n(r), S_n the integrated bump, "
    "a_n >= 0: rho >= 0 and non-increasing by construction).";
}  // namespace

std::string rhoform_absent_message() {
  return std::string("DM_REAM: RHOFORM is required: ") + kFormList;
}

std::optional<RhoForm> stored_rhoform(const tadah::core::Context &ctx) {
  const auto v = ctx.stored_value("RHOFORM");
  if (!v) return std::nullopt;
  return parse_rhoform(*v);
}

RhoForm require_rhoform(const tadah::core::Context &ctx,
                        const std::string &who) {
  const std::string prefix = who.empty() ? std::string() : who + ": ";
  const auto v = ctx.stored_value("RHOFORM");
  if (!v) throw DensityFormError(prefix + rhoform_absent_message());
  try {
    return parse_rhoform(*v);
  } catch (const DensityFormError &) {
    throw DensityFormError(prefix + "DM_REAM: RHOFORM '" + *v +
                           "' is not a density form: use " + kFormList);
  }
}

void normalise_rhoform_key(tadah::core::Context &c) {
  if (context_has_ream(c)) {
    const RhoForm f = require_rhoform(c, "");
    c.remove("RHOFORM");
    c.add("RHOFORM", std::string(rhoform_name(f)));
    return;
  }
  if (c.has_stored("RHOFORM")) {
    report_ream_only_key("RHOFORM");
    c.remove("RHOFORM");
  }
}

void check_rhoform_precedence(const tadah::core::Context &pot,
                              const tadah::core::Context &config,
                              const std::string &pot_name) {
  if (!context_has_ream(pot)) {
    // a key on either side has no effect (validated when it was read)
    if (pot.has_stored("RHOFORM") || config.has_stored("RHOFORM"))
      report_ream_only_key("RHOFORM");
    return;
  }
  // The potential must declare its own form: a configuration value never
  // supplies it.
  if (!pot.has_stored("RHOFORM") && config.has_stored("RHOFORM"))
    throw DensityFormError(
        pot_name + ": the RHOFORM of the configuration file does not apply "
        "(the form belongs to the potential, whose density amplitudes were "
        "fitted in it); the RHOFORM line belongs in " + pot_name + ". " +
        rhoform_absent_message());
  const RhoForm pv = require_rhoform(pot, pot_name);
  const auto cv = stored_rhoform(config);
  if (cv && *cv != pv)
    throw DensityFormError(
        "RHOFORM conflict: " + pot_name + " declares RHOFORM " +
        rhoform_name(pv) + " and the configuration file sets RHOFORM " +
        rhoform_name(*cv) +
        ". RHOFORM is a property of the potential (its density amplitudes "
        "were fitted in that form): the potential file sets it, not the "
        "configuration file.");
}

std::string rhoform_subject(RhoForm f) {
  return std::string("DM_REAM: RHOFORM ") + rhoform_name(f);
}

bool rhoform_cutoff_whitelisted(const std::string &label) {
  for (const char *l : kMonotoneCutoffs)
    if (label == l) return true;
  return false;
}

std::string rhoform_cutoff_whitelist_text() {
  const std::size_t n = std::size(kMonotoneCutoffs);
  std::string s;
  for (std::size_t i = 0; i < n; ++i)
    s += std::string(i == 0 ? "" : i + 1 == n ? " or " : ", ") +
         kMonotoneCutoffs[i];
  return s;
}

const char *monotone_knot_refusal_text() {
  return "a non-negative knot density is already non-increasing with a "
         "non-increasing cutoff";
}

OptimLowerBound optim_key_lower_bound(const std::string &key,
                                      std::optional<RhoForm> form) {
  if (form && rhoform_is_positive(*form) && upper(key) == "AMPGRIDMB")
    return OptimLowerBound::NonNegative;
  return OptimLowerBound::Free;
}

std::optional<OptimLowBoundRule> optim_low_bound_rule(
    const std::string &key, std::optional<RhoForm> form,
    std::optional<BasisKind> basis) {
  if (!form) return std::nullopt;
  const std::string by = std::string("RHOFORM ") + rhoform_name(*form);
  if (optim_key_lower_bound(key, form) == OptimLowerBound::NonNegative)
    return OptimLowBoundRule{false, "every AMPGRIDMB amplitude >= 0", by};
  if (!rhoform_is_positive(*form)) return std::nullopt;
  const std::string k = upper(key);
  if (k == "WATOM")
    return OptimLowBoundRule{true, "every density pair weight W > 0", by};
  if (k == "SGRIDMB" && *form == RhoForm::Monotone)
    return OptimLowBoundRule{true,
                             "every SGRIDMB width > 0 (S_n steps down across "
                             "the bump s_n (r - c_n))",
                             by};
  // ReamDensity::validate: NONNEG needs s > 0 for the Gaussian basis only
  if (k == "SGRIDMB" && basis == BasisKind::Gaussian)
    return OptimLowBoundRule{true,
                             "every SGRIDMB width > 0 for the Gaussian basis "
                             "(exp(-s (r - c)^2) grows without bound for "
                             "s < 0)",
                             by};
  return std::nullopt;
}

namespace {
/* optim_low_bound_rule(key, ctx) and, with @p entry, of that entry
 * (0-based) of the key. */
std::optional<OptimLowBoundRule> low_bound_rule_of(
    const std::string &key, const tadah::core::Context &ctx,
    std::optional<std::size_t> entry) {
  std::optional<RhoForm> form;
  std::optional<BasisKind> basis;
  if (context_has_ream(ctx)) {
    form = require_rhoform(ctx, "the OPTIM box rules");
    basis = optim_rule_basis(ctx);
  }
  if (auto r = optim_low_bound_rule(key, form, basis)) return r;
  // A radial bump channel: every width > 0 in every form.
  const std::string k = upper(key);
  if (k == "SGRIDMB" && basis && bump_halfwidth(*basis) > 0.0)
    return OptimLowBoundRule{
        true,
        "every SGRIDMB width > 0 (a radial bump spans [c - w/s, c + w/s]; a "
        "width <= 0 has no support edge)",
        std::string("the DM_REAM ") + basis_name(*basis) + " density basis"};
  if (k == "SGRID2B") {
    if (potedit::bump_channel_halfwidth(ctx, "2B")) {
      const auto label = ctx.stored_value("TYPE2B", 0);
      return OptimLowBoundRule{
          true,
          "every SGRID2B width > 0 (a radial bump spans [c - w/s, c + w/s]; "
          "a width <= 0 has no support edge)",
          "the " + label.value_or("TYPE2B") + " pair basis"};
    }
    // A D2_mJoin: a D2_Blip or D2_Blip5 descriptor holds its own SGRID2B
    // entries (D_mJoin::parse_context gives each descriptor its share, in
    // the order of the line), the entries of another descriptor keep that
    // descriptor's rule (a D2_BP's, none). Without an entry, the rule of
    // the first bump descriptor's entries.
    const auto line = ctx.stored_value("TYPE2B", 0);
    if (line && *line == "D2_mJoin") {
      std::size_t at = 0;  // the first entry of the descriptor
      for (const auto &child : D_mJoin::parse_context(ctx, "TYPE2B")) {
        const std::size_t n =
            child.contains("SGRID2B") ? child.size("SGRID2B") : 0;
        const bool bump =
            potedit::bump_channel_halfwidth(child, "2B").has_value();
        const bool owns = entry && *entry >= at && *entry < at + n;
        at += n;
        if (entry && !owns) continue;
        if (!bump) {
          if (entry) return std::nullopt;
          continue;
        }
        return OptimLowBoundRule{
            true,
            "each of its SGRID2B widths > 0 (a radial bump spans [c - w/s, "
            "c + w/s]; a width <= 0 has no support edge)",
            "the " + child.stored_value("TYPE2B", 0).value_or("TYPE2B") +
                " pair basis of the D2_mJoin line"};
      }
    }
  }
  return std::nullopt;
}
}  // namespace

std::optional<OptimLowBoundRule> optim_low_bound_rule(
    const std::string &key, const tadah::core::Context &ctx) {
  return low_bound_rule_of(key, ctx, std::nullopt);
}

std::optional<OptimLowBoundRule> optim_low_bound_rule(
    const std::string &key, const tadah::core::Context &ctx,
    std::size_t entry) {
  return low_bound_rule_of(key, ctx, entry);
}

std::optional<BasisKind> optim_rule_basis(const tadah::core::Context &ctx) {
  if (!context_has_ream(ctx)) return std::nullopt;
  auto basis_of = [](const tadah::core::Context &c) {
    try {
      return ReamDensity::basis_from_context(c);
    } catch (const std::runtime_error &e) {
      throw std::runtime_error(std::string("DM_REAM: ") + e.what());
    }
  };
  const auto label = ctx.stored_value("TYPEMB", 0);
  if (label && *label == "DM_REAM") return basis_of(ctx);
  // a DM_mJoin: its DM_REAM children
  std::optional<BasisKind> first;
  for (const auto &child : D_mJoin::parse_context(ctx, "TYPEMB")) {
    const auto l = child.stored_value("TYPEMB", 0);
    if (!l || *l != "DM_REAM") continue;
    const BasisKind b = basis_of(child);
    if (b == BasisKind::Gaussian) return b;
    if (!first) first = b;
  }
  return first;
}

bool OptimLowBoundRule::admits(double lo) const {
  return strict ? lo > 0.0 : lo >= 0.0;
}

}  // namespace models
}  // namespace tadah
