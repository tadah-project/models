#include "catch2/catch.hpp"
#include <tadah/models/eam/eam_setfl_io.h>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

// scaled_grid(scale): the output grid of tadah refit's OUTPUT_GRID_SCALE.
TEST_CASE("EAMSetflIO::scaled_grid refuses a grid a setfl table cannot have",
          "[eam]") {
  /* Pins: scaled_grid rounds nr * scale and nrho * scale, computed in
   * double, to the point counts and keeps the physical spans (n - 1) * d;
   * the smallest table it gives has 3 points (the spline that LAMMPS and
   * tadah fit through a setfl table reads its third point). A scale whose
   * counts fall below 3 or above 2147483647 (an int) throws
   * std::invalid_argument naming the scale and both counts, and a scale
   * that is not a positive finite number throws naming it.
   * Before the fix: each count was std::max(2, int(lround(n * scale))), so
   * a count below 2 became 2 silently and one above INT_MAX wrapped
   * negative and became 2 as well (lround is out of range above LONG_MAX):
   * scaled_grid(1e-4), (1e7) and (1e30) returned a 2-point grid, (0) and
   * (-1) too, and no exception is thrown. */
  EAMSetflIO io;
  io.read_setfl("tests_data/eam.alloy");
  REQUIRE(io.nr() == 5000);
  REQUIRE(io.nrho() == 5000);

  const EAMGrid g = io.scaled_grid(2.0);
  CHECK(g.nr == 10000);
  CHECK(g.nrho == 10000);
  CHECK(g.dr * (g.nr - 1) == Approx(io.dr() * (io.nr() - 1)));
  CHECK(g.drho * (g.nrho - 1) == Approx(io.drho() * (io.nrho() - 1)));
  CHECK(g.rcut == io.rcut());
  // 5000 x 0.0005 = 2.5 rounds to 3: the smallest table
  const EAMGrid g3 = io.scaled_grid(0.0005);
  CHECK(g3.nr == 3);
  CHECK(g3.nrho == 3);
  CHECK(g3.dr * 2 == Approx(io.dr() * (io.nr() - 1)));

  const auto refusal = [&](double scale) -> std::string {
    try {
      (void)io.scaled_grid(scale);
    } catch (const std::invalid_argument &e) {
      return e.what();
    }
    return "";
  };
  const std::string needs =
      " (the input table has 5000 and 5000); a setfl table needs 3 to "
      "2147483647 points (the spline that LAMMPS and tadah fit through it "
      "reads its third point)";
  CHECK(refusal(1e-4) == "EAMSetflIO: a grid scale of 0.0001 gives nrho = 1 "
                         "and nr = 1 points" + needs);
  // 5000 x 0.0003 = 1.4999999999999998 rounds to 1
  CHECK(refusal(0.0003) == "EAMSetflIO: a grid scale of 0.0003 gives nrho = "
                           "1 and nr = 1 points" + needs);
  CHECK(refusal(1e7) == "EAMSetflIO: a grid scale of 10000000 gives nrho = "
                        "50000000000 and nr = 50000000000 points" + needs);
  CHECK(refusal(1e30).rfind("EAMSetflIO: a grid scale of 1e+30 gives nrho = "
                            "5", 0) == 0);
  for (const double bad : {0.0, -1.0}) {
    INFO(bad);
    CHECK(refusal(bad).rfind("EAMSetflIO: a grid scale must be a positive "
                             "finite number, got ", 0) == 0);
  }
}

// ---- read_setfl refuses a malformed setfl file -----------------------------

namespace {

// A valid single-element setfl file with 3-point tables (the smallest the
// spline reads), one line per header record and per table: line 4 the
// element line, 5 the grid, 6 the element record, 7-9 F(rho), rho(r) and
// z2r(r).
std::vector<std::string> setfl3_lines() {
  return {"comment 1",
          "comment 2",
          "comment 3",
          "1 Ta",
          "3 0.5 3 1.0 2.0",
          "73 180.95 3.304 bcc",
          "0 -1 -1.5",
          "1 0.5 0",
          "2 1 0"};
}

std::string joined(const std::vector<std::string> &lines) {
  std::string text;
  for (const std::string &l : lines) text += l + "\n";
  return text;
}

// The message read_setfl_text throws for `text` (named "bad.alloy"), or
// "no exception".
std::string setfl_refusal(const std::string &text) {
  EAMSetflIO io;
  try {
    io.read_setfl_text(text, "bad.alloy");
  } catch (const std::exception &e) {
    return e.what();
  }
  return "no exception";
}

}  // namespace

TEST_CASE("EAMSetflIO reads the smallest setfl table, 3 points",
          "[eam][setfl_read]") {
  /* Pins: a 3-point setfl file reads (nrho, drho, nr, dr, rcut, Z, mass,
   * a, lattice, the tables) and splines, and the splines pass through the
   * table values (F at rho = drho is the second F value, -1; rho(r) at
   * r = 2 dr the third, 0).
   * Before the fix: passes as well (a pin of the floor the fix keeps). */
  EAMSetflIO io;
  io.read_setfl_text(joined(setfl3_lines()), "ok.alloy");
  CHECK(io.nrho() == 3);
  CHECK(io.drho() == 0.5);
  CHECK(io.nr() == 3);
  CHECK(io.dr() == 1.0);
  CHECK(io.rcut() == 2.0);
  CHECK(io.atomic_number() == 73);
  CHECK(io.atomic_mass() == 180.95);
  CHECK(io.lattice_param() == 3.304);
  CHECK(io.lattice() == "bcc");
  CHECK(io.element() == "Ta");
  CHECK(io.frho() == std::vector<double>{0, -1, -1.5});
  CHECK(io.rhor() == std::vector<double>{1, 0.5, 0});
  CHECK(io.z2r() == std::vector<double>{2, 1, 0});
  io.build_splines();
  CHECK(io.F_tab(0.5) == -1.0);
  CHECK(io.rho_tab(2.0) == 0.0);
}

TEST_CASE("EAMSetflIO skips '#' comments and blank lines after line 3, as "
          "LAMMPS does",
          "[eam][setfl_read]") {
  /* Pins: from line 4 on, '#' starts a comment and a line with no words is
   * skipped (LAMMPS' TextFileReader::next_line, which pair eam/alloy reads
   * the file with); lines 1-3 are comments whatever they hold. A file with
   * comment lines, trailing comments, blank lines and CRLF line ends reads
   * as the plain file does, value for value.
   * Before the fix: the commented file was refused, its line 4 being a
   * comment: "EAMSetflIO: malformed setfl header in 'commented.alloy' (line
   * 4 must start with the number of element types)." (the CRLF file read). */
  EAMSetflIO plain;
  plain.read_setfl_text(joined(setfl3_lines()), "plain.alloy");

  const std::string commented =
      "comment 1 # lines 1-3 are comments whatever they hold\n"
      "\n"
      "\n"
      "# the element line follows\n"
      "1 Ta   # Ta only\n"
      "\n"
      "3 0.5 3 1.0 2.0\n"
      "73 180.95 3.304 bcc # Z mass a lattice\n"
      "# F(rho)\n"
      "0 -1 -1.5\n"
      "\n"
      "1 0.5 0 # rho(r)\n"
      "2 1 0\n"
      "\n";
  std::string crlf;
  for (const std::string &l : setfl3_lines()) crlf += l + "\r\n";

  for (const auto &[name, text] :
       {std::pair<std::string, std::string>{"commented.alloy", commented},
        {"crlf.alloy", crlf}}) {
    INFO(name);
    EAMSetflIO io;
    REQUIRE_NOTHROW(io.read_setfl_text(text, name));
    CHECK(io.nrho() == plain.nrho());
    CHECK(io.drho() == plain.drho());
    CHECK(io.nr() == plain.nr());
    CHECK(io.dr() == plain.dr());
    CHECK(io.rcut() == plain.rcut());
    CHECK(io.atomic_number() == plain.atomic_number());
    CHECK(io.atomic_mass() == plain.atomic_mass());
    CHECK(io.lattice_param() == plain.lattice_param());
    CHECK(io.lattice() == plain.lattice());
    CHECK(io.element() == plain.element());
    CHECK(io.frho() == plain.frho());
    CHECK(io.rhor() == plain.rhor());
    CHECK(io.z2r() == plain.z2r());
  }
}

TEST_CASE("EAMSetflIO reads a header record that spans lines, as LAMMPS "
          "joins it",
          "[eam][setfl_read]") {
  /* Pins: a header record may span lines, as LAMMPS' next_values joins
   * lines until they hold the record's words: the grid record written as
   * '3 0.5 3' and '1.0 2.0', and the element record as '73' and '180.95
   * 3.304 bcc', read the values of the plain file, value for value.
   * Before the fix: passes as well (a pin: operator>> read across lines). */
  EAMSetflIO plain;
  plain.read_setfl_text(joined(setfl3_lines()), "plain.alloy");
  std::vector<std::string> l = setfl3_lines();
  l[4] = "3 0.5 3";
  l.insert(l.begin() + 5, "1.0 2.0");
  l[6] = "73";
  l.insert(l.begin() + 7, "180.95 3.304 bcc");
  REQUIRE(l.size() == 11);
  EAMSetflIO io;
  REQUIRE_NOTHROW(io.read_setfl_text(joined(l), "split.alloy"));
  CHECK(io.nrho() == plain.nrho());
  CHECK(io.drho() == plain.drho());
  CHECK(io.nr() == plain.nr());
  CHECK(io.dr() == plain.dr());
  CHECK(io.rcut() == plain.rcut());
  CHECK(io.atomic_number() == plain.atomic_number());
  CHECK(io.atomic_mass() == plain.atomic_mass());
  CHECK(io.lattice_param() == plain.lattice_param());
  CHECK(io.lattice() == plain.lattice());
  CHECK(io.element() == plain.element());
  CHECK(io.frho() == plain.frho());
  CHECK(io.rhor() == plain.rhor());
  CHECK(io.z2r() == plain.z2r());
}

TEST_CASE("EAMSetflIO: a refused file leaves an object that holds a file as "
          "it was",
          "[eam][setfl_read]") {
  /* Pins: read_setfl_text of a refused file into an EAMSetflIO that holds a
   * file throws and leaves it as it was: after the 3-point Ta file
   * (ok.alloy), a file of another element and grid (Nb, Z 41, nrho 3 drho
   * 0.25, nr 3 dr 0.5 rcut 1.0) that ends before its z2r table is refused,
   * and the object still has ok.alloy's path, element, header and tables,
   * and splines them.
   * Before the fix: the refused read had written the path, the element and
   * the header into the object before it threw ("EAMSetflIO: truncated
   * tabulated data in: bad.alloy"): path 'bad.alloy', element 'Nb', Z 41,
   * drho 0.25 and so on, with the tables of ok.alloy. */
  EAMSetflIO io;
  io.read_setfl_text(joined(setfl3_lines()), "ok.alloy");
  std::vector<std::string> bad = setfl3_lines();
  bad[3] = "1 Nb";
  bad[4] = "3 0.25 3 0.5 1.0";
  bad[5] = "41 92.906 3.30 bcc";
  bad.pop_back();  // no z2r table
  CHECK_THROWS_AS(io.read_setfl_text(joined(bad), "bad.alloy"),
                  std::runtime_error);
  CHECK(io.path() == "ok.alloy");
  CHECK(io.element() == "Ta");
  CHECK(io.atomic_number() == 73);
  CHECK(io.atomic_mass() == 180.95);
  CHECK(io.lattice_param() == 3.304);
  CHECK(io.lattice() == "bcc");
  CHECK(io.nrho() == 3);
  CHECK(io.drho() == 0.5);
  CHECK(io.nr() == 3);
  CHECK(io.dr() == 1.0);
  CHECK(io.rcut() == 2.0);
  CHECK(io.frho() == std::vector<double>{0, -1, -1.5});
  CHECK(io.rhor() == std::vector<double>{1, 0.5, 0});
  CHECK(io.z2r() == std::vector<double>{2, 1, 0});
  io.build_splines();
  CHECK(io.F_tab(0.5) == -1.0);
}

TEST_CASE("EAMSetflIO refuses a malformed setfl file, naming the file, the "
          "line and the field",
          "[eam][setfl_read]") {
  /* Pins: read_setfl_text throws std::runtime_error with exactly these
   * messages for each defect of the 3-point file: a count that is not a
   * whole number from 3 to 2147483647 (the spline that LAMMPS and tadah fit
   * through a table reads its third point), a spacing or cutoff that is not
   * a positive finite number, a header record with a field too many or too
   * few, a table value that is not a finite decimal number, a table shorter
   * than its count (also a count absurd for the file, without allocating
   * it), and a table whose last line holds values of the next one (LAMMPS
   * drops them). Lines after the z2r table are not read, as LAMMPS does not
   * read them.
   * Before the fix: read without an exception: nrho 2 and nr 1 (a 1- or
   * 2-point table makes build_splines read past its end), drho 0, dr -1.0,
   * rcut 0, an element line '1 Ta W', '1' (no element) or '1.0 Ta', an
   * element record with a fifth word (read as the first F value: every
   * table shifted by one value) and the values past a table's end; '3.0',
   * '2147483648', nan, inf, 1e400, '1.0D-03', '0x1p0', a 6-word grid line,
   * a 3-word element record and a file that ends failed as "EAMSetflIO:
   * malformed setfl header in: bad.alloy" or "EAMSetflIO: truncated
   * tabulated data in: bad.alloy" (a file of 3 lines: "... (line 4 must
   * start with the number of element types)."), naming neither the line
   * nor the field; nr -3 threw std::length_error ("vector::_M_default_append").
   * A line after the z2r table was not read either (that check passes). */
  const std::string at = "EAMSetflIO: setfl file 'bad.alloy', line ";
  const std::string spline =
      " (the cubic spline that LAMMPS and tadah fit through a table reads "
      "its third point)";
  struct Case {
    int line;               // 1-based line of setfl3_lines() replaced
    std::string text;       // its new text
    std::string message;
  };
  const std::vector<Case> cases = {
      {4, "1 Ta W",
       at + "4: the element line of a file of 1 element type is '1 <symbol>', "
            "it has 3 words"},
      {4, "1",
       at + "4: the element line of a file of 1 element type is '1 <symbol>', "
            "it has 1 word"},
      {4, "1.0 Ta",
       at + "4: the element line must start with the number of element "
            "types, got '1.0'"},
      {5, "2 0.5 3 1.0 2.0",
       at + "5: nrho must be a whole number from 3 to 2147483647" + spline +
           ", got '2'"},
      {5, "3 0.5 1 1.0 2.0",
       at + "5: nr must be a whole number from 3 to 2147483647" + spline +
           ", got '1'"},
      {5, "3 0.5 -3 1.0 2.0",
       at + "5: nr must be a whole number from 3 to 2147483647" + spline +
           ", got '-3'"},
      {5, "3.0 0.5 3 1.0 2.0",
       at + "5: nrho must be a whole number from 3 to 2147483647" + spline +
           ", got '3.0'"},
      {5, "2147483648 0.5 3 1.0 2.0",
       at + "5: nrho must be a whole number from 3 to 2147483647" + spline +
           ", got '2147483648'"},
      {5, "3 0 3 1.0 2.0",
       at + "5: drho must be a positive finite number, got '0'"},
      {5, "3 0.5 3 -1.0 2.0",
       at + "5: dr must be a positive finite number, got '-1.0'"},
      {5, "3 0.5 3 nan 2.0",
       at + "5: dr must be a positive finite number, got 'nan'"},
      {5, "3 0.5 3 1.0 0",
       at + "5: rcut must be a positive finite number, got '0'"},
      {5, "3 0.5 3 1.0 1e400",
       at + "5: rcut must be a positive finite number, got '1e400'"},
      {5, "3 0.5 3 1.0 2.0 7",
       at + "5: the grid line is 'nrho drho nr dr rcut', it has 6 words"},
      {6, "73 180.95 3.304",
       at + "6: the element record is 'Z mass a lattice', it has 6 words "
            "on lines 6-7"},
      {6, "73 180.95 3.304 bcc 1",
       at + "6: the element record is 'Z mass a lattice', it has 5 words"},
      {6, "73.0 180.95 3.304 bcc",
       at + "6: the atomic number must be a whole number, got '73.0'"},
      {6, "73 heavy 3.304 bcc",
       at + "6: the mass must be a finite number, got 'heavy'"},
      {6, "73 180.95 inf bcc",
       at + "6: the lattice constant must be a finite number, got 'inf'"},
      {7, "0 nan -1.5",
       at + "7: value 2 of the F(rho) table, 'nan', is not a finite number"},
      {8, "1 0.5 1.0D-03",
       at + "8: value 3 of the rho(r) table, '1.0D-03', is not a finite "
            "number"},
      {9, "2 0x1p0 0",
       at + "9: value 2 of the z2r(r) table, '0x1p0', is not a finite "
            "number"},
      {9, "2 1 1e400",
       at + "9: value 3 of the z2r(r) table, '1e400', is not a finite "
            "number"},
      {7, "0 -1 -1.5 1",
       at + "7: the F(rho) table of nrho = 3 values ends with word 3 of this "
            "line, which has 4: each table starts on a new line (LAMMPS "
            "drops the rest of a table's last line)"},
      {8, "1 0.5 0 2 1 0",
       at + "8: the rho(r) table of nr = 3 values ends with word 3 of this "
            "line, which has 6: each table starts on a new line (LAMMPS "
            "drops the rest of a table's last line)"},
      {9, "2 1 0 0",
       at + "9: the z2r(r) table of nr = 3 values ends with word 3 of this "
            "line, which has 4: each table starts on a new line (LAMMPS "
            "drops the rest of a table's last line)"},
  };
  for (const Case &c : cases) {
    std::vector<std::string> lines = setfl3_lines();
    lines[static_cast<std::size_t>(c.line - 1)] = c.text;
    INFO("line " << c.line << ": '" << c.text << "'");
    CHECK(setfl_refusal(joined(lines)) == c.message);
  }

  // The file ends: inside a table, before a header record, and with a count
  // the file cannot hold (no table of that size is allocated).
  const std::string ends = "EAMSetflIO: setfl file 'bad.alloy': the file ends ";
  std::vector<std::string> lines = setfl3_lines();
  lines.pop_back();
  CHECK(setfl_refusal(joined(lines)) ==
        ends + "inside the z2r(r) table, after 0 of its nr = 3 values");
  lines = setfl3_lines();
  lines[8] = "2 1";
  CHECK(setfl_refusal(joined(lines)) ==
        ends + "inside the z2r(r) table, after 2 of its nr = 3 values");
  lines = setfl3_lines();
  lines[4] = "1000000 0.5 3 1.0 2.0";
  CHECK(setfl_refusal(joined(lines)) ==
        ends + "inside the F(rho) table, after 9 of its nrho = 1000000 "
               "values");
  lines = setfl3_lines();
  lines.resize(4);
  CHECK(setfl_refusal(joined(lines)) ==
        ends + "before the grid line ('nrho drho nr dr rcut')");
  lines.resize(3);
  CHECK(setfl_refusal(joined(lines)) ==
        ends + "before the element line (the number of element types and "
               "their names)");
  lines = setfl3_lines();
  lines.resize(5);
  lines.push_back("73 180.95");
  CHECK(setfl_refusal(joined(lines)) ==
        ends + "inside the element record ('Z mass a lattice'), which has "
               "2 of its 4 words");

  // Lines after the z2r table are not read (LAMMPS does not read them).
  lines = setfl3_lines();
  lines.push_back("this line is not read");
  CHECK(setfl_refusal(joined(lines)) == "no exception");
}
