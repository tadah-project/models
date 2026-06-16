#include "catch2/catch.hpp"
#include <tadah/models/eam/eam_lattice.h>

#include <algorithm>
#include <cmath>

using namespace tadah::models;

TEST_CASE("eam_lattice nearest-neighbour distances", "[eam_lattice]") {
  const double a = 3.30;
  REQUIRE(eam_lattice::nn_distance("bcc", a) ==
          Approx(a * std::sqrt(3.0) / 2.0));
  REQUIRE(eam_lattice::nn_distance("fcc", a) == Approx(a / std::sqrt(2.0)));
  REQUIRE(eam_lattice::nn_distance("hcp", a) == Approx(a));
  // unknown falls back to bcc
  REQUIRE(eam_lattice::nn_distance("dia", a) ==
          Approx(a * std::sqrt(3.0) / 2.0));
}

TEST_CASE("eam_lattice neighbour multiplicities and counts", "[eam_lattice]") {
  const double a = 3.0;
  auto count_at = [](const std::vector<double> &d, double r) {
    return std::count_if(d.begin(), d.end(), [&](double x) {
      return std::fabs(x - r) < 1e-6;
    });
  };

  SECTION("bcc: 8 first + 6 second neighbours") {
    auto d = eam_lattice::neighbour_distances("bcc", a, 1.05 * a);
    REQUIRE(count_at(d, a * std::sqrt(3.0) / 2.0) == 8);
    REQUIRE(count_at(d, a) == 6);
    REQUIRE(d.size() == 14);
  }
  SECTION("fcc: 12 first neighbours") {
    auto d = eam_lattice::neighbour_distances("fcc", a, 0.99 * a);
    REQUIRE(count_at(d, a / std::sqrt(2.0)) == 12);
    REQUIRE(d.size() == 12);
  }
  SECTION("hcp (ideal c/a): 12 first neighbours") {
    auto d = eam_lattice::neighbour_distances("hcp", a, 1.01 * a);
    REQUIRE(count_at(d, a) == 12);
    REQUIRE(d.size() == 12);
  }
  SECTION("all distances within cutoff and positive") {
    for (const char *lat : {"bcc", "fcc", "hcp"}) {
      auto d = eam_lattice::neighbour_distances(lat, a, 2.5 * a);
      REQUIRE(!d.empty());
      for (double x : d) {
        REQUIRE(x > 0.0);
        REQUIRE(x <= 2.5 * a);
      }
    }
  }
}

TEST_CASE("eam_lattice supported()", "[eam_lattice]") {
  REQUIRE(eam_lattice::supported("bcc"));
  REQUIRE(eam_lattice::supported("fcc"));
  REQUIRE(eam_lattice::supported("hcp"));
  REQUIRE(!eam_lattice::supported("dia"));
  REQUIRE(!eam_lattice::supported(""));
  // supported() expects the canonical lower-case form (the setfl reader
  // normalises the header label); raw upper case is intentionally not matched.
  REQUIRE(!eam_lattice::supported("FCC"));
}

TEST_CASE("eam_lattice crystal_epa", "[eam_lattice]") {
  const double a = 3.0;
  auto one = [](double) { return 1.0; };   // phi = rho = 1
  auto idF = [](double x) { return x; };    // F(x) = x
  // With phi=rho=1 and F=id, E/atom = 0.5*sum(phi) + F(sum(rho)) = 1.5 * Nnbr,
  // so crystal_epa reduces to a clean neighbour-count check.
  SECTION("fcc: 12 first neighbours -> 1.5*12") {
    REQUIRE(eam_lattice::crystal_epa("fcc", a, 0.99 * a, one, one, idF) ==
            Approx(1.5 * 12));
  }
  SECTION("bcc: 8+6 neighbours -> 1.5*14") {
    REQUIRE(eam_lattice::crystal_epa("bcc", a, 1.05 * a, one, one, idF) ==
            Approx(1.5 * 14));
  }
  SECTION("hcp (ideal c/a): 12 first neighbours -> 1.5*12") {
    REQUIRE(eam_lattice::crystal_epa("hcp", a, 1.01 * a, one, one, idF) ==
            Approx(1.5 * 12));
  }
  SECTION("unsupported lattice -> NaN") {
    REQUIRE(std::isnan(eam_lattice::crystal_epa("dia", a, 2.0 * a, one, one,
                                                idF)));
  }
}

TEST_CASE("eam_lattice require_crystal", "[eam_lattice]") {
  // accepts fcc/bcc/hcp with a0>0, case-insensitively
  REQUIRE_NOTHROW(eam_lattice::require_crystal("fcc", 3.6));
  REQUIRE_NOTHROW(eam_lattice::require_crystal("FCC", 3.6));
  REQUIRE_NOTHROW(eam_lattice::require_crystal("BCC", 2.8));
  REQUIRE_NOTHROW(eam_lattice::require_crystal("HCP", 2.9));
  // rejects unspecified / unsupported lattice, and non-positive a0
  REQUIRE_THROWS(eam_lattice::require_crystal("NONE", 0.0));
  REQUIRE_THROWS(eam_lattice::require_crystal("", 3.6));
  REQUIRE_THROWS(eam_lattice::require_crystal("dia", 5.0));
  REQUIRE_THROWS(eam_lattice::require_crystal("fcc", 0.0));
}
