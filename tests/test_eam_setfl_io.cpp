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

// The analytic spline derivatives (eval_deriv / eval_deriv2, and the phi=z2r/r
// quotient rule) must match central finite differences of the value functions
// — sampled at cell midpoints with a small step so [r-h, r+h] stays inside one
// cubic cell (where the 2nd-difference is essentially exact and the
// 1st-difference truncation is tiny). Guards against sign / factor / missing
// 1/r mistakes.
TEST_CASE("EAMSetflIO analytic derivatives match finite differences", "[eam]") {
  EAMSetflIO io;
  io.read_setfl("tests_data/eam.alloy");
  io.build_splines();

  const double dr = io.dr();
  const double h = 0.02 * dr;
  int checked = 0;
  for (int i = 3; i < io.nr() - 3; i += 5) {
    const double r = (i + 0.5) * dr;          // mid-cell
    if (r - h <= 0.0 || r + h >= io.rcut()) continue;
    const double drho_fd = (io.rho_tab(r + h) - io.rho_tab(r - h)) / (2 * h);
    const double d2rho_fd =
        (io.rho_tab(r + h) - 2 * io.rho_tab(r) + io.rho_tab(r - h)) / (h * h);
    REQUIRE(io.drho_tab(r) == Approx(drho_fd).epsilon(0.02).margin(1e-6));
    REQUIRE(io.d2rho_tab(r) == Approx(d2rho_fd).epsilon(0.02).margin(1e-5));
    const double dphi_fd = (io.phi_tab(r + h) - io.phi_tab(r - h)) / (2 * h);
    const double d2phi_fd =
        (io.phi_tab(r + h) - 2 * io.phi_tab(r) + io.phi_tab(r - h)) / (h * h);
    REQUIRE(io.dphi_tab(r) == Approx(dphi_fd).epsilon(0.02).margin(1e-6));
    REQUIRE(io.d2phi_tab(r) == Approx(d2phi_fd).epsilon(0.02).margin(1e-5));
    ++checked;
  }
  REQUIRE(checked > 5);

  const double dq = io.drho();
  const double hq = 0.02 * dq;
  for (int i = 3; i < io.nrho() - 3; i += 5) {
    const double rho = (i + 0.5) * dq;
    if (rho - hq <= 0.0 || rho + hq >= io.rhomax()) continue;
    const double dF_fd = (io.F_tab(rho + hq) - io.F_tab(rho - hq)) / (2 * hq);
    const double d2F_fd =
        (io.F_tab(rho + hq) - 2 * io.F_tab(rho) + io.F_tab(rho - hq)) / (hq * hq);
    REQUIRE(io.dF_tab(rho) == Approx(dF_fd).epsilon(0.02).margin(1e-6));
    REQUIRE(io.d2F_tab(rho) == Approx(d2F_fd).epsilon(0.02).margin(1e-5));
  }
}

// The single-element guard must reject a multi-element setfl header.
TEST_CASE("EAMSetflIO rejects missing file", "[eam]") {
  EAMSetflIO io;
  REQUIRE_THROWS(io.read_setfl("tests_data/does_not_exist.alloy"));
}
