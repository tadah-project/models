#include "catch2/catch.hpp"
#include <tadah/models/descriptors/d_base.h>
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/periodic_table.h>
#include <tadah/core/core_user_config_parser.h>
#include <tadah/core/number_format.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

namespace {
struct WatomTestFixture {
  WatomTestFixture() { PeriodicTable::initialize(); }
};

Context make_minimal_blip_context_watom_tests() {
  Context c;
  c.add("INIT2B", "true");
  c.add("VERBOSE", 0);
  c.add("SGRID2B", "2.0 1.0");
  c.add("CGRID2B", "0.5 1.4");
  c.add("RCUT2B", 3.5);
  c.add("TYPE2B", "D2_Blip");
  c.add("TYPE2B", 2);
  c.add("TYPE2B", 2);
  c.add("TYPE2B", "H"); c.add("TYPE2B", "H");
  return c;
}
}  // namespace
TEST_CASE( "Testing D_Base", "[D_Base]" ) {

  v_type sgrid(5);
  v_type cgrid(6);
  Context c;
  c.add("SGRID2B", "LIN 5 0.1 12");
  c.add("SGRIDMB", "LIN 5 0.1 12");
  c.add("CGRID2B", "LIN 5 0.1 12");
  c.add("CGRIDMB", "LIN 5 0.1 12");
  D_Base::get_grid(c, "SGRID2B", sgrid );
  std::vector<double> sg1 = {0.1,3.075,6.05,9.025,12};
  REQUIRE_THAT(sgrid, Catch::Matchers::Equals(sg1));
  D_Base::get_grid(c, "SGRIDMB", sgrid );


  D_Base::get_grid(c, "CGRID2B", cgrid );
  D_Base::get_grid(c, "CGRIDMB", cgrid );
  std::vector<double> cg1 = {0.1,3.075,6.05,9.025,12};
  REQUIRE_THAT(cgrid, Catch::Matchers::Equals(cg1));

  double margin = 1e-10;
  Context c2;
  c2.add("SGRID2B", "GEOM 5 0.1 12");
  c2.add("SGRIDMB", "GEOM 5 0.1 12");
  c2.add("CGRID2B", "GEOM 5 0.1 12");
  c2.add("CGRIDMB", "GEOM 5 0.1 12");
  D_Base::get_grid(c2, "SGRID2B", sgrid );
  D_Base::get_grid(c2, "SGRIDMB", sgrid );
  std::vector<double> sg2 = {0.1, 0.330975092, 1.095445115, 3.6256504768, 12.0 };
  REQUIRE_THAT(sgrid, Catch::Matchers::Approx(sg2).margin(margin));

  D_Base::get_grid(c2, "CGRID2B", cgrid );
  D_Base::get_grid(c2, "CGRIDMB", cgrid );
  std::vector<double> cg2 = {0.1, 0.330975092, 1.095445115, 3.6256504768, 12.0 };
  REQUIRE_THAT(cgrid, Catch::Matchers::Approx(cg2).margin(margin));

  Context c3;
  c3.add("SGRID2B", "-2.1 5 0.1 12");
  D_Base::get_grid(c3, "SGRID2B", sgrid );

  //throws
  Context c_err;
  c_err.add("SGRID2B", "LOG 5 0.1 12");
  REQUIRE_THROWS(D_Base::get_grid(c_err,"SGRID2B",sgrid));
}

// =====================================================================
// Dual-format WATOM tests. Exercise D_Base::update_weights /
// get_atom_Z_for_arg_pos / has_watom_arg_pos / dweight_dw_n through a
// concrete D2_Blip descriptor.
// =====================================================================

TEST_CASE_METHOD(WatomTestFixture, "WATOM per-species single element: diagonal write",
          "[D_Base][watom][dual_format]") {
  Context c = make_minimal_blip_context_watom_tests();
  c.add("ATOM", "H");
  c.add("WATOM", 0.5);
  D2_Blip d(c);

  CHECK(d.weights[1][1] == Approx(0.5));
  CHECK(d.weights[8][8] == Approx(1.0));  // unrelated element stays at default
  CHECK(d.has_watom_arg_pos(1));
  CHECK_FALSE(d.has_watom_arg_pos(2));
  CHECK(d.get_atom_Z_for_arg_pos(1) == 1);
  CHECK(d.dweight_dw_n(1, 1, 1) == Approx(1.0));
}

TEST_CASE_METHOD(WatomTestFixture, "WATOM per-species two elements: geometric-mean off-diagonal",
          "[D_Base][watom][dual_format]") {
  Context c = make_minimal_blip_context_watom_tests();
  c.add("ATOM", "H");
  c.add("ATOM", "O");
  c.add("WATOM", 0.5);
  c.add("WATOM", 2.0);
  D2_Blip d(c);

  CHECK(d.weights[1][1] == Approx(0.5));
  CHECK(d.weights[8][8] == Approx(2.0));
  const double cross = std::sqrt(0.5 * 2.0);
  CHECK(d.weights[1][8] == Approx(cross));
  CHECK(d.weights[8][1] == Approx(cross));
  CHECK(d.has_watom_arg_pos(1));
  CHECK(d.has_watom_arg_pos(2));
  CHECK_FALSE(d.has_watom_arg_pos(3));
  CHECK(d.get_atom_Z_for_arg_pos(1) == 1);
  CHECK(d.get_atom_Z_for_arg_pos(2) == 8);

  // Chain factors at arg_pos=1 (Z_target=H, w_H=0.5, w_O=2.0):
  //   d w[H][H] / d w_H = 1
  //   d w[O][O] / d w_H = 0
  //   d w[H][O] / d w_H = 0.5 * sqrt(w_O / w_H) = 0.5 * sqrt(4) = 1.0
  CHECK(d.dweight_dw_n(1, 1, 1) == Approx(1.0));
  CHECK(d.dweight_dw_n(8, 8, 1) == Approx(0.0));
  CHECK(d.dweight_dw_n(1, 8, 1) == Approx(0.5 * std::sqrt(2.0 / 0.5)));
  CHECK(d.dweight_dw_n(8, 1, 1) == Approx(0.5 * std::sqrt(2.0 / 0.5)));
  // arg_pos=2 (Z_target=O):
  CHECK(d.dweight_dw_n(8, 8, 2) == Approx(1.0));
  CHECK(d.dweight_dw_n(1, 1, 2) == Approx(0.0));
  CHECK(d.dweight_dw_n(1, 8, 2) == Approx(0.5 * std::sqrt(0.5 / 2.0)));
}

TEST_CASE_METHOD(WatomTestFixture, "WATOM pair-triple: direct symmetric writes, no geometric mean",
          "[D_Base][watom][dual_format]") {
  Context c = make_minimal_blip_context_watom_tests();
  c.add("ATOM", "H");
  c.add("ATOM", "O");
  // Pair triples in arbitrary order.
  c.add("WATOM", "H"); c.add("WATOM", "H"); c.add("WATOM", 0.7);
  c.add("WATOM", "O"); c.add("WATOM", "O"); c.add("WATOM", 1.3);
  c.add("WATOM", "H"); c.add("WATOM", "O"); c.add("WATOM", 0.9);
  D2_Blip d(c);

  CHECK(d.weights[1][1] == Approx(0.7));
  CHECK(d.weights[8][8] == Approx(1.3));
  CHECK(d.weights[1][8] == Approx(0.9));  // direct, NOT sqrt(0.7*1.3)
  CHECK(d.weights[8][1] == Approx(0.9));
  CHECK(d.weights[6][6] == Approx(1.0));  // unrelated element default

  CHECK(d.has_watom_arg_pos(1));
  CHECK(d.has_watom_arg_pos(2));
  CHECK(d.has_watom_arg_pos(3));
  CHECK_FALSE(d.has_watom_arg_pos(4));

  // PairTriple mode: arg_pos→Z mapping is empty.
  CHECK(d.get_atom_Z_for_arg_pos(1) == 0);

  // Chain factors in pair-triple mode:
  //   arg_pos=1 → (H,H); 1.0 only on the (H,H) pair.
  CHECK(d.dweight_dw_n(1, 1, 1) == Approx(1.0));
  CHECK(d.dweight_dw_n(8, 8, 1) == Approx(0.0));
  CHECK(d.dweight_dw_n(1, 8, 1) == Approx(0.0));
  //   arg_pos=3 → (H,O); 1.0 on (H,O) and (O,H) only.
  CHECK(d.dweight_dw_n(1, 8, 3) == Approx(1.0));
  CHECK(d.dweight_dw_n(8, 1, 3) == Approx(1.0));
  CHECK(d.dweight_dw_n(1, 1, 3) == Approx(0.0));
}

TEST_CASE_METHOD(WatomTestFixture, "WATOM pair-triple: bad token count throws",
          "[D_Base][watom][dual_format]") {
  Context c = make_minimal_blip_context_watom_tests();
  c.add("ATOM", "H");
  c.add("ATOM", "O");
  // Two tokens — not a multiple of three.
  c.add("WATOM", "H"); c.add("WATOM", "H");
  REQUIRE_THROWS_AS(D2_Blip(c), std::runtime_error);
}

TEST_CASE_METHOD(WatomTestFixture, "WATOM default (no WATOM key) leaves pair weights at 1.0",
          "[D_Base][watom][dual_format]") {
  Context c = make_minimal_blip_context_watom_tests();
  // No ATOM, no WATOM: update_weights returns early, default-filled.
  D2_Blip d(c);
  for (int i = 1; i < 119; ++i) {
    CHECK(d.weights[i][i] == Approx(1.0));
  }
  CHECK_FALSE(d.has_watom_arg_pos(1));
  CHECK(d.dweight_dw_n(1, 1, 1) == Approx(0.0));  // Empty mode
}

TEST_CASE_METHOD(WatomTestFixture, "WATOM single-species H reduces to old per-Zj semantics",
          "[D_Base][watom][dual_format][regression]") {
  // Phase 2C-WATOM regression: single H species → diagonal chain == 1.0,
  // recovers the old δ(Zj, Z_target) gate.
  Context c = make_minimal_blip_context_watom_tests();
  c.add("ATOM", "H");
  c.add("WATOM", 1.0);
  D2_Blip d(c);
  CHECK(d.dweight_dw_n(1, 1, 1) == Approx(1.0));
  CHECK(d.weights[1][1] == Approx(1.0));
}

TEST_CASE("D_Base::get_grid reads the subnormal values the lossless writers "
          "produce", "[D_Base][lossless]") {
  // Context::add(double) and write_potential() now keep a subnormal (the
  // old fixed-decimal rendering flushed it to 0); std::stod threw
  // out_of_range on it, so get_grid rejected the grid.
  const double sub = 1e-310;
  const double dmin = std::numeric_limits<double>::denorm_min();
  SECTION("plain tokens") {
    const v_type g = D_Base::get_grid(std::vector<std::string>{
        "1.0", format_roundtrip(sub), "2.0", format_roundtrip(-dmin)});
    REQUIRE(g.size() == 4);
    CHECK(g[1] == sub);
    CHECK(g[3] == -dmin);
  }
  SECTION("LIN/GEOM endpoints") {
    const v_type g = D_Base::get_grid(
        std::vector<std::string>{"LIN", "3", "0", "1e-310"});
    REQUIRE(g.size() == 3);
    CHECK(g.back() == sub);
  }
  SECTION("Context add -> write_potential -> reparse -> get_grid") {
    Context c;
    c.add("CGRID2B", std::vector<double>{1.0, sub, 2.0});
    std::ostringstream os;
    c.write_potential(os, 18);
    std::istringstream in(os.str());
    const Context back = CoreUserConfigParser().parse(in);
    v_type g;
    D_Base::get_grid(back, "CGRID2B", g);
    REQUIRE(g.size() == 3);
    CHECK(g[0] == 1.0);
    CHECK(g[1] == sub);
    CHECK(g[2] == 2.0);
  }
  SECTION("non-numbers and overflow are still rejected") {
    CHECK_THROWS_AS(D_Base::get_grid(std::vector<std::string>{"abc"}),
                    std::invalid_argument);
    CHECK_THROWS_AS(D_Base::get_grid(std::vector<std::string>{"1e400"}),
                    std::invalid_argument);
    CHECK_THROWS_AS(D_Base::get_grid(std::vector<std::string>{"1.5x"}),
                    std::invalid_argument);
  }
}
