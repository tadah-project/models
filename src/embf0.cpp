#include <tadah/models/descriptors/dm/embf0.h>
#include <tadah/core/logger.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <stdexcept>

namespace tadah {
namespace models {

namespace {
std::string upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char ch) { return std::toupper(ch); });
  return s;
}
std::atomic<bool> warned_stray{false};
}  // namespace

const char *embf0_name(EmbF0Mode m) {
  return m == EmbF0Mode::Free ? "FREE" : "ZERO";
}

EmbF0Mode parse_embf0(const std::string &value) {
  const std::string v = upper(value);
  if (v == "ZERO") return EmbF0Mode::Zero;
  if (v == "FREE") return EmbF0Mode::Free;
  throw std::runtime_error(
      "EMBF0 must be ZERO or FREE (got '" + value +
      "'). ZERO (default) evaluates the DM_REAM embedding with F(0) = 0 by "
      "construction; FREE is the raw embedding of tadah 1.5.0 and earlier.");
}

EmbF0Mode embf0_from_context(const tadah::core::Context &ctx) {
  return parse_embf0(ctx.get<std::string>("EMBF0"));
}

std::optional<EmbF0Mode> stored_embf0(const tadah::core::Context &ctx) {
  const auto v = ctx.stored_value("EMBF0");
  if (!v) return std::nullopt;
  return parse_embf0(*v);
}

bool context_has_ream(const tadah::core::Context &ctx) {
  const auto initmb = ctx.stored_value("INITMB");
  if (!initmb || upper(*initmb) != "TRUE") return false;
  for (std::size_t n = 0;; ++n) {
    const auto tok = ctx.stored_value("TYPEMB", n);
    if (!tok) return false;
    if (*tok == "DM_REAM") return true;
  }
}

bool model_is_linear(const tadah::core::Context &ctx) {
  const auto bf = ctx.stored_value("MODEL", 1);
  return bf && (*bf == "BF_Linear" || *bf == "Kern_Linear");
}

std::string embf0_absent_consequence(bool linear) {
  if (linear)
    return "this potential is evaluated with the default EMBF0 ZERO (F(0) = 0 "
           "by construction). A potential written by tadah 1.5.0 or earlier "
           "now gives every atomic energy shifted by -F(0) (forces and "
           "stresses are unchanged for this linear model); add 'EMBF0 FREE' "
           "to the potential to reproduce its energies.";
  return "this potential is evaluated with the default EMBF0 ZERO (F(0) = 0 "
         "by construction). The model is non-linear in the descriptor, so for "
         "a potential written by tadah 1.5.0 or earlier its energies AND "
         "forces change; add 'EMBF0 FREE' to the potential to reproduce them.";
}

void normalise_embf0_key(tadah::core::Context &c) {
  if (context_has_ream(c)) {
    const EmbF0Mode m = stored_embf0(c).value_or(EmbF0Mode::Zero);
    c.remove("EMBF0");
    c.add("EMBF0", std::string(embf0_name(m)));
    return;
  }
  if (c.has_stored("EMBF0")) {
    c.remove("EMBF0");
    if (!warned_stray.exchange(true))
      tadah::core::Logger::instance().log(
          tadah::core::LogLevel::WARNING,
          "EMBF0 is ignored and not written: the potential has no DM_REAM "
          "many-body block.");
  }
}

void check_embf0_precedence(const tadah::core::Context &pot,
                            const tadah::core::Context &task,
                            const std::string &pot_name) {
  const auto tv = stored_embf0(task);
  if (!tv) return;
  if (!context_has_ream(pot)) {
    tadah::core::Logger::instance().log(
        tadah::core::LogLevel::WARNING,
        "EMBF0 " + std::string(embf0_name(*tv)) + " is ignored: " + pot_name +
            " has no DM_REAM many-body block.");
    return;
  }
  const auto pv = stored_embf0(pot);
  if (pv && *pv != *tv)
    throw std::runtime_error(
        "EMBF0 conflict: " + pot_name + " declares EMBF0 " +
        embf0_name(*pv) + " but the task configuration sets EMBF0 " +
        embf0_name(*tv) +
        ". EMBF0 is a property of the potential (its weights were fitted in "
        "that embedding reference); edit the potential file instead.");
  if (!pv)
    tadah::core::Logger::instance().log(
        tadah::core::LogLevel::INFO,
        "EMBF0 " + std::string(embf0_name(*tv)) +
            " from the task configuration applies to " + pot_name +
            " (the potential does not declare EMBF0).");
}

}  // namespace models
}  // namespace tadah
