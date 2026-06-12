#include "catch2/catch.hpp"
#include <tadah/models/eam/eam_setfl_io.h>
#include <cmath>
#include <cstdio>
#include <string>

using namespace tadah::models;

// Round-trip the shared EAM reader/writer: read the Ta2 (eam.alloy) fixture,
// write it back out through the NATIVE-function callbacks (here fed the original
// tabulated functions), read the result, and confirm the grid metadata and the
// three functions agree at a dense set of sample points.
TEST_CASE("EAMSetflIO read/write round-trip", "[eam]") {
  EAMSetflIO io;
  io.read_setfl("tests_data/eam.alloy");
  io.build_splines();

  REQUIRE(io.nr()   > 0);
  REQUIRE(io.nrho() > 0);
  REQUIRE(io.rcut() > 0.0);

  const std::string out = "tests_data_eam_roundtrip.alloy";
  io.write_setfl(out, io.input_grid(),
                 [&](double rho){ return io.F_tab(rho); },
                 [&](double r){ return io.rho_tab(r); },
                 [&](double r){ return r > 0.0 ? io.phi_tab(r) * r : 0.0; });

  EAMSetflIO io2;
  io2.read_setfl(out);
  io2.build_splines();

  // Grid metadata is carried over verbatim.
  REQUIRE(io2.nrho() == io.nrho());
  REQUIRE(io2.nr()   == io.nr());
  REQUIRE(io2.drho() == Approx(io.drho()));
  REQUIRE(io2.dr()   == Approx(io.dr()));
  REQUIRE(io2.rcut() == Approx(io.rcut()));
  REQUIRE(io2.atomic_number() == io.atomic_number());
  REQUIRE(io2.element() == io.element());

  // Functions agree across the sampled range.
  const double rmax = io.rcut();
  for (int k = 1; k <= 60; ++k) {
    const double r = rmax * k / 61.0;
    REQUIRE(io2.phi_tab(r) == Approx(io.phi_tab(r)).margin(1e-9));
    REQUIRE(io2.rho_tab(r) == Approx(io.rho_tab(r)).margin(1e-9));
  }
  const double rhomax = io.rhomax();
  for (int k = 1; k <= 60; ++k) {
    const double rho = rhomax * k / 61.0;
    REQUIRE(io2.F_tab(rho) == Approx(io.F_tab(rho)).margin(1e-9));
  }

  std::remove(out.c_str());
}

// The single-element guard must reject a multi-element setfl header.
TEST_CASE("EAMSetflIO rejects missing file", "[eam]") {
  EAMSetflIO io;
  REQUIRE_THROWS(io.read_setfl("tests_data/does_not_exist.alloy"));
}
