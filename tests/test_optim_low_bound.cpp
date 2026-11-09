// Catch2 v2 tests for the low-bound rule of an OPTIM box that comes from
// the lower limits the descriptors declare (ParamMeta::hard.lo):
//   * declared_low_bound_rule(key): the limit every registered descriptor
//     declares for the values of the key, one per key;
//   * optim_box_low_bound_rule(key, ctx[, entry]): the density-form and
//     bump rules (optim_low_bound_rule) first, then the limit the
//     descriptor of the configuration that reads the key declares;
//   * OptimLowBoundRule::clamp_low: where an AUTO box that reaches below
//     the rule is clamped.

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/registry.h>
#include <tadah/models/descriptors/d_base.h>
#include <tadah/models/descriptors/dm/rhoform.h>
#include <tadah/models/descriptors/param_meta.h>

#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>

using namespace tadah::models;
using Catch::Matchers::Contains;
using tadah::core::Context;

namespace {

/* A configuration of one two-body descriptor @p type2b ("D2_BP 2 2 H H")
 * and, when @p typemb is not empty, one many-body descriptor. */
Context lines(const std::string &type2b, const std::string &typemb = "") {
  Context c;
  c.add("ATOM", "H");
  c.add("WATOM", 1.0);
  c.add("INIT2B", "true");
  c.add("TYPE2B", type2b);
  c.add("CGRID2B", "1.0 2.0");
  c.add("SGRID2B", "0.5 1.0");
  c.add("RCTYPE2B", "Cut_Cos");
  c.add("RCUT2B", 5.0);
  if (typemb.empty()) {
    c.add("INITMB", "false");
    return c;
  }
  c.add("INITMB", "true");
  c.add("TYPEMB", typemb);
  c.add("CGRIDMB", 0.0);
  c.add("SGRIDMB", 0.2);
  c.add("CEMBFUNC", "1.0 2.0");
  c.add("SEMBFUNC", "1.0 1.0");
  c.add("RCUTENV", 1.5);
  c.add("RCTYPEMB", "Cut_Cos");
  c.add("RCUTMB", 5.0);
  return c;
}

/* The bits of @p x, so that a comparison holds under -ffast-math. */
std::uint64_t bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

}  // namespace

/* Pins: the lower limits the registered descriptors declare for the values
 * of a key (ParamMeta::hard.lo of a parameter with arg_pos <= 0, the
 * values of the key itself and not of a TYPE line) are one value per key,
 * none below 0: 0 for the distance cutoffs RCUT2B and RCUTMB and the
 * centres CGRID2B and CGRIDMB, kHardLoPositive (1e-6, the limit of a value
 * that must be > 0) for the widths SGRID2B, SGRIDMB and SEMBFUNC and for
 * RCUTENV, and 3 for the knot degrees. declared_low_bound_rule(key) gives
 * the rule of the key (case-insensitive): > 0 (strict, bound 0, "every
 * RCUTENV value > 0") for a key whose limit is kHardLoPositive, >= its
 * limit for another ("every CGRID2B value >= 0"), and none for a key no
 * descriptor limits (CEMBFUNC, AMPGRIDMB, WATOM, FIXWEIGHT, ESIGMA); its
 * subject names the first registered descriptor that declares it.
 * Before the fix: declared_low_bound_rule does not exist, and DM_REAM
 * declared -1e30 for CEMBFUNC and AMPGRIDMB (no limit, as a number). */
TEST_CASE("The descriptors declare one lower limit per key, which gives the "
          "OPTIM box rule of the key",
          "[metadata][optim_rule][declared]") {
  std::map<std::string, std::set<double>> declared;
  for (const auto &label : tadah::core::Registry<D_Base>::names()) {
    std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>(label));
    const auto meta = d->meta();
    if (!meta.has_value()) continue;
    for (const auto &pm : meta->params)
      if (pm.arg_pos <= 0 && pm.hard.lo) declared[pm.key].insert(*pm.hard.lo);
  }
  const std::map<std::string, double> expected = {
      {"RCUT2B", 0.0},          {"RCUTMB", 0.0},
      {"CGRID2B", 0.0},         {"CGRIDMB", 0.0},
      {"SGRID2B", kHardLoPositive}, {"SGRIDMB", kHardLoPositive},
      {"SEMBFUNC", kHardLoPositive}, {"RCUTENV", kHardLoPositive},
      {"KDEG2B", 3.0},          {"KDEGMB", 3.0},
      {"KDEGEMB", 3.0}};
  for (const auto &[key, values] : declared) {
    CAPTURE(key);
    CHECK(values.size() == 1);
    REQUIRE(expected.count(key) == 1);
    CHECK(bits(*values.begin()) == bits(expected.at(key)));
  }
  CHECK(declared.size() == expected.size());
  CHECK(bits(kHardLoPositive) == bits(1.0e-6));

  for (const char *key : {"RCUTENV", "rcutenv", "SGRID2B", "SGRIDMB",
                          "SEMBFUNC"}) {
    CAPTURE(key);
    const auto r = declared_low_bound_rule(key);
    REQUIRE(r);
    CHECK(r->strict);
    CHECK(bits(r->bound) == bits(0.0));
    CHECK_FALSE(r->admits(0.0));
    CHECK_FALSE(r->admits(-1.0));
    CHECK(r->admits(1e-9));
  }
  CHECK(declared_low_bound_rule("RCUTENV")->requirement ==
        "every RCUTENV value > 0");
  CHECK(declared_low_bound_rule("RCUTENV")->by ==
        "the DM_mBlipEnv descriptor");
  for (const char *key : {"RCUT2B", "RCUTMB", "CGRID2B", "CGRIDMB"}) {
    CAPTURE(key);
    const auto r = declared_low_bound_rule(key);
    REQUIRE(r);
    CHECK_FALSE(r->strict);
    CHECK(bits(r->bound) == bits(0.0));
    CHECK(r->admits(0.0));
    CHECK_FALSE(r->admits(-1e-9));
    CHECK(r->requirement == std::string("every ") + key + " value >= 0");
  }
  const auto kdeg = declared_low_bound_rule("KDEG2B");
  REQUIRE(kdeg);
  CHECK_FALSE(kdeg->strict);
  CHECK(kdeg->admits(3.0));
  CHECK_FALSE(kdeg->admits(2.5));
  CHECK(kdeg->requirement == "every KDEG2B value >= 3");
  for (const char *key :
       {"CEMBFUNC", "AMPGRIDMB", "WATOM", "FIXWEIGHT", "ESIGMA", "LAMBDA"}) {
    CAPTURE(key);
    CHECK_FALSE(declared_low_bound_rule(key));
  }
}

/* Pins: optim_box_low_bound_rule(key, ctx) is the density-form or bump
 * rule of the key (optim_low_bound_rule(key, ctx)) when there is one, else
 * the limit that a descriptor of the configuration that reads the key
 * declares, its subject naming that descriptor:
 *   - DM_mBlipEnv: RCUTENV, SGRIDMB and SEMBFUNC > 0 ("the DM_mBlipEnv
 *     descriptor requires every RCUTENV value > 0"), CGRIDMB and RCUTMB
 *     >= 0, CEMBFUNC none;
 *   - D2_BP: SGRID2B > 0 by the D2_BP descriptor, CGRID2B >= 0; D2_Blip:
 *     SGRID2B by the bump rule of its pair basis;
 *   - a key that no descriptor of the configuration reads (SGRIDMB without
 *     a many-body term) has none, whatever other descriptors declare;
 *   - on a D2_mJoin line of a D2_BP and a D2_Blip, the entry form gives
 *     each entry the rule of its descriptor: the D2_BP's SGRID2B entries
 *     the D2_BP limit, the D2_Blip's the bump rule, and the rule of the
 *     key without an entry is the bump rule.
 * Before the fix: optim_box_low_bound_rule does not exist; tadah hpo took
 * the boxes of RCUTENV, SGRIDMB of DM_mBlipEnv, SGRID2B of D2_BP and the
 * centres below 0 without a word. */
TEST_CASE("optim_box_low_bound_rule: the form and bump rules first, then "
          "the limit the descriptor of the configuration declares",
          "[metadata][optim_rule][declared]") {
  const Context env = lines("D2_BP 2 2 H H", "DM_mBlipEnv 0 1 1 2 2 H H");
  for (const char *key : {"RCUTENV", "SGRIDMB", "SEMBFUNC", "SGRID2B"}) {
    CAPTURE(key);
    const auto r = optim_box_low_bound_rule(key, env);
    REQUIRE(r);
    CHECK(r->strict);
    CHECK(r->requirement == std::string("every ") + key + " value > 0");
    CHECK(r->by == std::string(std::strcmp(key, "SGRID2B") == 0
                                   ? "the D2_BP descriptor"
                                   : "the DM_mBlipEnv descriptor"));
  }
  for (const char *key : {"CGRIDMB", "RCUTMB", "CGRID2B", "RCUT2B"}) {
    CAPTURE(key);
    const auto r = optim_box_low_bound_rule(key, env);
    REQUIRE(r);
    CHECK_FALSE(r->strict);
    CHECK(r->requirement == std::string("every ") + key + " value >= 0");
  }
  CHECK_FALSE(optim_box_low_bound_rule("CEMBFUNC", env));
  CHECK(optim_box_low_bound_rule("rcutenv", env));

  const Context bp = lines("D2_BP 2 2 H H");
  CHECK_FALSE(optim_box_low_bound_rule("SGRIDMB", bp));
  CHECK_FALSE(optim_box_low_bound_rule("RCUTENV", bp));

  const Context blip = lines("D2_Blip 2 2 H H");
  const auto b = optim_box_low_bound_rule("SGRID2B", blip);
  REQUIRE(b);
  CHECK(b->strict);
  CHECK(b->by == "the D2_Blip pair basis");

  Context join = lines("D2_mJoin D2_BP 2 2 H H D2_Blip 2 2 H H");
  join.remove("CGRID2B");
  join.remove("SGRID2B");
  join.remove("RCUT2B");
  join.remove("RCTYPE2B");
  join.add("CGRID2B", "0 0 4.0 3.0");
  join.add("SGRID2B", "0.5 1.0 2.0 1.0");
  join.add("RCUT2B", "5.0 5.0");
  join.add("RCTYPE2B", "Cut_Cos Cut_Cos");
  for (std::size_t entry : {0u, 1u}) {
    const auto r = optim_box_low_bound_rule("SGRID2B", join, entry);
    REQUIRE(r);
    CHECK(r->strict);
    CHECK(r->by == "the D2_BP descriptor");
  }
  for (std::size_t entry : {2u, 3u}) {
    const auto r = optim_box_low_bound_rule("SGRID2B", join, entry);
    REQUIRE(r);
    CHECK(r->by == "the D2_Blip pair basis of the D2_mJoin line");
  }
  const auto key_rule = optim_box_low_bound_rule("SGRID2B", join);
  REQUIRE(key_rule);
  CHECK(key_rule->by == "the D2_Blip pair basis of the D2_mJoin line");
  const auto centre = optim_box_low_bound_rule("CGRID2B", join, 2);
  REQUIRE(centre);
  CHECK_FALSE(centre->strict);
  CHECK(centre->by == "the D2_Blip descriptor");
}

/* Pins: OptimLowBoundRule::clamp_low(seed), the low end of an AUTO box of
 * the seed that reaches below the rule: the bound of a >= rule (0 for
 * CGRID2B), and for a > 0 rule kHardLoPositive (1e-6), or half the seed
 * for a seed at or below 2e-6, so the clamped box holds the seed; admits
 * takes NaN for no rule.
 * Before the fix: clamp_low and the bound of a rule do not exist (an AUTO
 * box below a rule was refused). */
TEST_CASE("OptimLowBoundRule::clamp_low: at the bound, at 1e-6 for a value "
          "that must be > 0",
          "[optim_rule][declared]") {
  const auto pos = declared_low_bound_rule("RCUTENV");
  REQUIRE(pos);
  CHECK(bits(pos->clamp_low(5.0)) == bits(1e-6));
  CHECK(bits(pos->clamp_low(0.2)) == bits(1e-6));
  CHECK(bits(pos->clamp_low(1e-6)) == bits(5e-7));
  CHECK(pos->admits(pos->clamp_low(1e-6)));
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(pos->admits(nan));
  const auto nonneg = declared_low_bound_rule("CGRID2B");
  REQUIRE(nonneg);
  CHECK(bits(nonneg->clamp_low(0.3)) == bits(0.0));
  CHECK(nonneg->admits(nonneg->clamp_low(0.3)));
  CHECK_FALSE(nonneg->admits(nan));
  const auto kdeg = declared_low_bound_rule("KDEGMB");
  REQUIRE(kdeg);
  CHECK(bits(kdeg->clamp_low(5.0)) == bits(3.0));
}
