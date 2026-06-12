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
}
