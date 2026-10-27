// test_potential_lines.cpp
//
// A potential is read as it was written: the RCTYPE*, RCUT* and list values
// of its own lines that its descriptors do not take are reported once per
// key and left out (D_mJoin::parse_context, test_potential_list_values.cpp).
// A configuration that includes the potential with CONFIG is the
// concatenation of its lines: a line outside the potential's own lines (a
// line of the including file, above or below the CONFIG line) adds its
// values to the key, and the key is read as in a training configuration:
// values past those the descriptors take are refused, naming the lines that
// set the key and their files. The refusal of a count names those lines for
// a training configuration read from a file as well.

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/core_user_config_parser.h>
#include <tadah/core/key_notice.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/dc_selector.h>
#include <tadah/models/descriptors/d_mjoin.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <vector>
#include <unistd.h>

using tadah::core::Context;
using tadah::core::CoreUserConfigParser;
using tadah::core::IgnoredKeyNoticeCapture;
using tadah::models::D_mJoin;
using tadah::models::DC_Selector;

namespace {

/* A scratch directory unique to this process, removed with its files when
 * the object goes out of scope. */
class Scratch {
public:
  Scratch()
      : dir_(std::filesystem::temp_directory_path() /
             ("tadah_models_potential_lines_" + std::to_string(::getpid()))) {
    std::error_code ec;
    std::filesystem::remove_all(dir_, ec);
    std::filesystem::create_directories(dir_);
  }
  ~Scratch() {
    std::error_code ec;
    std::filesystem::remove_all(dir_, ec);
  }
  Scratch(const Scratch &) = delete;
  Scratch &operator=(const Scratch &) = delete;

  /* Writes @p body to <dir>/<name>; returns its path. */
  std::string write(const std::string &name, const std::string &body) const {
    const auto p = dir_ / name;
    std::ofstream os(p);
    os << body;
    os.close();
    REQUIRE_FALSE(os.fail());
    return p.string();
  }

private:
  std::filesystem::path dir_;
};

/* A single D2_Blip potential: RCUT2B 5 on line 6. */
const char *kSingle2B = "ATOM Ta\n"
                        "INIT2B true\n"
                        "INITMB false\n"
                        "TYPE2B D2_Blip 2 2 Ta Ta\n"
                        "RCTYPE2B Cut_Cos\n"
                        "RCUT2B 5\n"
                        "CGRID2B 2 3\n"
                        "SGRID2B 1 1\n"
                        "WEIGHTS 0.1 0.2 0.3\n";

/* A single DM_EAD potential: RCUTMB 5 on line 6. */
const char *kSingleMB = "ATOM Ta\n"
                        "INIT2B false\n"
                        "INITMB true\n"
                        "TYPEMB DM_EAD 0 2 2 Ta Ta\n"
                        "RCTYPEMB Cut_Cos\n"
                        "RCUTMB 5\n"
                        "CGRIDMB 2 3\n"
                        "SGRIDMB 1 1\n"
                        "WEIGHTS 0.1 0.2 0.3\n";

/* A D2_mJoin potential, D2_Blip with RCUT2B 5.3 and D2_BP with 4.5: RCUT2B
 * on line 8, CGRID2B on line 9. With @p old one more RCUT2B value, 6, than
 * the line takes, as tadah 1.5.0 wrote them. */
std::string joined_2b(bool old = false) {
  return std::string("ATOM Ta\n"
                     "INIT2B true\n"
                     "INITMB false\n"
                     "TYPE2B D2_mJoin\n"
                     "TYPE2B D2_Blip 2 2 Ta Ta\n"
                     "TYPE2B D2_BP 2 2 Ta Ta\n"
                     "RCTYPE2B Cut_Cos Cut_Cos\n") +
         (old ? "RCUT2B 5.3 4.5 6\n" : "RCUT2B 5.3 4.5\n") +
         "CGRID2B 2 3 2.5 3.5\n"
         "SGRID2B 1 1 0.5 0.5\n"
         "WEIGHTS 0.1 0.2 0.3 0.4 0.5\n";
}

/* A DM_mJoin potential of two DM_EAD descriptors: RCTYPEMB on line 7,
 * RCUTMB 4.5 5 on line 8. */
const char *kJoinedMB = "ATOM Ta\n"
                        "INIT2B false\n"
                        "INITMB true\n"
                        "TYPEMB DM_mJoin\n"
                        "TYPEMB DM_EAD 0 2 2 Ta Ta\n"
                        "TYPEMB DM_EAD 0 1 1 Ta Ta\n"
                        "RCTYPEMB Cut_Cos Cut_Cos\n"
                        "RCUTMB 4.5 5\n"
                        "CGRIDMB 2 3 2.5\n"
                        "SGRIDMB 1 1 1\n"
                        "WEIGHTS 0.1 0.2 0.3 0.4\n";

Context parsed(const std::string &file) {
  tadah::core::PeriodicTable::initialize();
  return CoreUserConfigParser().parse(file);
}

/* The message of D_mJoin::expand_grids (tadah train, tadah hpo and the
 * analysis commands) on @p c; empty when it throws none. */
std::string refusal(Context c) {
  try {
    D_mJoin::expand_grids(c);
  } catch (const std::exception &e) {
    return e.what();
  }
  return "";
}

/* The message of a model built from @p c (DC_Selector: tadah predict,
 * pair_style tadah); empty when it throws none. */
std::string load_refusal(const Context &c) {
  try {
    DC_Selector s(c);
  } catch (const std::exception &e) {
    return e.what();
  }
  return "";
}

std::string joined(const std::vector<std::string> &v) {
  std::string s;
  for (const auto &x : v) s += x + "\n";
  return s;
}

const char *kOneRcut2B =
    ": each descriptor takes one RCUT2B value, in the order of the TYPE2B "
    "lines";

} // namespace

/* A line of RCUT2B in a file that includes a single-descriptor potential
 * with CONFIG, above or below the CONFIG line, gives RCUT2B two values, and
 * the descriptor takes one: refused as in a training configuration, naming
 * both lines and their files, for training (expand_grids) and for a model
 * built from the configuration (DC_Selector). A smaller cutoff than the
 * potential's alike.
 * Before the fix: accepted with 'Part of RCUT2B (5) is ignored: TYPE2B
 * D2_Blip takes one RCUT2B value, and its cutoff is the largest, 6; the
 * potential works unchanged with RCUT2B 6 alone': the largest cutoff was
 * the model's (4 below the CONFIG line was ignored). */
TEST_CASE("A line of RCUT2B in a file that includes a single-descriptor "
          "potential is refused naming both lines",
          "[D_mJoin][potential][config_include][potential_lines]") {
  Scratch dir;
  const std::string pot = dir.write("single.tadah", kSingle2B);
  SECTION("above the CONFIG line") {
    const std::string above =
        dir.write("above.train", "RCUT2B 6\nCONFIG " + pot + "\n");
    const std::string expected =
        "TYPE2B names 1 descriptor (D2_Blip) but RCUT2B has 2 values, set "
        "on line 1 of " + above + " and on line 6 of " + pot +
        " (CONFIG on line 2 of " + above + ")" + kOneRcut2B;
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(above)) == expected);
    CHECK(load_refusal(parsed(above)) == expected);
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
  SECTION("below the CONFIG line, a smaller cutoff") {
    const std::string below =
        dir.write("below.train", "CONFIG " + pot + "\nRCUT2B 4\n");
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(below)) ==
          "TYPE2B names 1 descriptor (D2_Blip) but RCUT2B has 2 values, set "
          "on line 6 of " + pot + " (CONFIG on line 1 of " + below +
              ") and on line 2 of " + below + kOneRcut2B);
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
  SECTION("a list key the descriptor does not read") {
    const std::string kdeg =
        dir.write("kdeg.train", "KDEG2B 5 5\nCONFIG " + pot + "\n");
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(kdeg)) ==
          "TYPE2B D2_Blip takes 0 KDEG2B values, but KDEG2B has 2 values, set "
          "on line 1 of " + kdeg + ": TYPE2B D2_Blip 2 2 Ta Ta");
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
}

/* The same for a single many-body descriptor: RCUTMB below the CONFIG line.
 * Before the fix: 'Part of RCUTMB (4.5) is ignored: TYPEMB DM_EAD takes one
 * RCUTMB value, and its cutoff is the largest, 5; ...', accepted. */
TEST_CASE("A line of RCUTMB in a file that includes a single many-body "
          "potential is refused naming both lines",
          "[D_mJoin][potential][config_include][potential_lines]") {
  Scratch dir;
  const std::string pot = dir.write("single_mb.tadah", kSingleMB);
  const std::string below =
      dir.write("below_mb.train", "CONFIG " + pot + "\nRCUTMB 4.5\n");
  IgnoredKeyNoticeCapture cap;
  CHECK(refusal(parsed(below)) ==
        "TYPEMB names 1 descriptor (DM_EAD) but RCUTMB has 2 values, set on "
        "line 6 of " + pot + " (CONFIG on line 1 of " + below +
            ") and on line 2 of " + below +
            ": each descriptor takes one RCUTMB value, in the order of the "
            "TYPEMB lines");
  INFO(joined(cap.notices()));
  CHECK(cap.notices().empty());
}

/* A D2_mJoin potential (D2_Blip 5.3, D2_BP 4.5) included with CONFIG:
 * 'RCUT2B 6' above the CONFIG line would give the descriptors 6 and 5.3,
 * below it 6 would be read by none; a CGRID2B line below it gives the line
 * two values past those its descriptors take. Each is refused, naming the
 * lines that set the key.
 * Before the fix: accepted, 'Part of RCUT2B (4.5) is ignored: TYPE2B names
 * 2 descriptors (D2_Blip, D2_BP), which take its first 2 values; the
 * potential works unchanged without it' (D2_BP's cutoff 5.3), 'Part of
 * RCUT2B (6) ...' and 'Part of CGRID2B (4 5) ...'. */
TEST_CASE("A line of a descriptor key in a file that includes a D2_mJoin "
          "potential is refused naming the lines",
          "[D_mJoin][potential][config_include][potential_lines]") {
  Scratch dir;
  const std::string pot = dir.write("joined.tadah", joined_2b());
  const std::string named =
      "TYPE2B names 2 descriptors (D2_Blip, D2_BP)";
  SECTION("RCUT2B above the CONFIG line") {
    const std::string above =
        dir.write("above.train", "RCUT2B 6.0\nCONFIG " + pot + "\n");
    const std::string expected =
        named + " but RCUT2B has 3 values, set on line 1 of " + above +
        " and on line 8 of " + pot + " (CONFIG on line 2 of " + above + ")" +
        kOneRcut2B;
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(above)) == expected);
    CHECK(load_refusal(parsed(above)) == expected);
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
  SECTION("RCUT2B below the CONFIG line") {
    const std::string below =
        dir.write("below.train", "CONFIG " + pot + "\nRCUT2B 6.0\n");
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(below)) ==
          named + " but RCUT2B has 3 values, set on line 8 of " + pot +
              " (CONFIG on line 1 of " + below + ") and on line 2 of " +
              below + kOneRcut2B);
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
  SECTION("CGRID2B below the CONFIG line") {
    const std::string grid =
        dir.write("grid.train", "CONFIG " + pot + "\nCGRID2B 4 5\n");
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(grid)) ==
          named + ", which take 4 CGRID2B values, but CGRID2B has 6 values, "
                  "set on line 9 of " + pot + " (CONFIG on line 1 of " +
              grid + ") and on line 2 of " + grid +
              ": each descriptor takes its CGRID2B values, in the order of "
              "the TYPE2B lines");
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
}

/* A DM_mJoin potential (two DM_EAD) included with CONFIG: an RCUTMB line
 * below the CONFIG line and an RCTYPEMB line above it are refused naming
 * the lines.
 * Before the fix: accepted, 'Part of RCUTMB (6) is ignored: ...' and 'Part
 * of RCTYPEMB (Cut_Cos) is ignored: ...' (the descriptors took Cut_Tanh and
 * the potential's first Cut_Cos). */
TEST_CASE("A line of a descriptor key in a file that includes a DM_mJoin "
          "potential is refused naming the lines",
          "[D_mJoin][potential][config_include][potential_lines]") {
  Scratch dir;
  const std::string pot = dir.write("joined_mb.tadah", kJoinedMB);
  const std::string named = "TYPEMB names 2 descriptors (DM_EAD, DM_EAD)";
  SECTION("RCUTMB below the CONFIG line") {
    const std::string below =
        dir.write("below_mb.train", "CONFIG " + pot + "\nRCUTMB 6\n");
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(below)) ==
          named + " but RCUTMB has 3 values, set on line 8 of " + pot +
              " (CONFIG on line 1 of " + below + ") and on line 2 of " +
              below +
              ": each descriptor takes one RCUTMB value, in the order of the "
              "TYPEMB lines");
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
  SECTION("RCTYPEMB above the CONFIG line") {
    const std::string above =
        dir.write("above_mb.train", "RCTYPEMB Cut_Tanh\nCONFIG " + pot + "\n");
    IgnoredKeyNoticeCapture cap;
    CHECK(refusal(parsed(above)) ==
          named + " but RCTYPEMB has 3 values, set on line 1 of " + above +
              " and on line 7 of " + pot + " (CONFIG on line 2 of " + above +
              "): each descriptor takes one RCTYPEMB value, in the order of "
              "the TYPEMB lines");
    INFO(joined(cap.notices()));
    CHECK(cap.notices().empty());
  }
}

/* A TYPE2B line below the CONFIG line of a D2_mJoin potential of tadah
 * 1.5.0 (RCUT2B 5.3 4.5 6 7, two values past those its two descriptors
 * take) adds a third descriptor, with its RCTYPE2B, CGRID2B and SGRID2B
 * lines: the potential's RCUT2B values are then read by descriptors that
 * are not the potential's own, as in a configuration, and four values for
 * three descriptors are refused naming the potential's line.
 * Before the fix: accepted with 'Part of RCUT2B (7) is ignored: TYPE2B
 * names 3 descriptors (D2_Blip, D2_BP, D2_BP), which take its first 3
 * values; the potential works unchanged without it', the third descriptor
 * taking the potential's unread 6 as its cutoff. */
TEST_CASE("A descriptor that a file including a potential adds to its line "
          "reads the potential's values as a configuration does",
          "[D_mJoin][potential][config_include][potential_lines]") {
  Scratch dir;
  std::string body = joined_2b(true);
  body.replace(body.find("RCUT2B 5.3 4.5 6\n"), 17, "RCUT2B 5.3 4.5 6 7\n");
  const std::string pot = dir.write("old2.tadah", body);
  const std::string ext = dir.write(
      "ext.train", "CONFIG " + pot +
                       "\nTYPE2B D2_BP 2 2 Ta Ta\nRCTYPE2B Cut_Cos\n"
                       "CGRID2B 2 3\nSGRID2B 1 1\n");
  IgnoredKeyNoticeCapture cap;
  CHECK(refusal(parsed(ext)) ==
        "TYPE2B names 3 descriptors (D2_Blip, D2_BP, D2_BP) but RCUT2B has 4 "
        "values, set on line 8 of " + pot + " (CONFIG on line 1 of " + ext +
            ")" + kOneRcut2B);
  INFO(joined(cap.notices()));
  CHECK(cap.notices().empty());
}

/* Pins: the values of a potential's own lines that its descriptors do not
 * take are reported and ignored, read alone (tadah predict -p, pair_style
 * tadah) or included with CONFIG by a file whose lines set other keys
 * (training from the potential): a D2_mJoin potential of tadah 1.5.0 with
 * RCUT2B 5.3 4.5 6 loads with 'Part of RCUT2B (6) is ignored: TYPE2B names
 * 2 descriptors (D2_Blip, D2_BP), which take its first 2 values; the
 * potential works unchanged without it', and training from it keeps
 * RCUT2B 5.3 4.5; a single D2_Blip potential with RCUT2B 5 6 loads with its
 * largest cutoff, 6.
 * Before the fix: passes as well. */
TEST_CASE("A potential's own values that its descriptors do not take are "
          "reported, read alone or included with CONFIG",
          "[D_mJoin][potential][config_include][potential_lines]") {
  Scratch dir;
  const std::string pot = dir.write("old.tadah", joined_2b(true));
  const std::string notice =
      "Part of RCUT2B (6) is ignored: TYPE2B names 2 descriptors (D2_Blip, "
      "D2_BP), which take its first 2 values; the potential works unchanged "
      "without it";
  SECTION("read alone") {
    IgnoredKeyNoticeCapture cap;
    CHECK(load_refusal(parsed(pot)).empty());
    CHECK(cap.notices() == std::vector<std::string>{notice});
  }
  SECTION("included with CONFIG") {
    const std::string run = dir.write(
        "run.train", "DBFILE db.tadah\nCONFIG " + pot + "\nLAMBDA 0\n");
    IgnoredKeyNoticeCapture cap;
    Context c = parsed(run);
    REQUIRE_NOTHROW(D_mJoin::expand_grids(c));
    CHECK(c("RCUT2B") == std::vector<std::string>{"5.3", "4.5"});
    CHECK(cap.notices() == std::vector<std::string>{notice});
  }
  SECTION("a single descriptor read alone") {
    std::string body = kSingle2B;
    body.replace(body.find("RCUT2B 5\n"), 9, "RCUT2B 5 6\n");
    const std::string single = dir.write("old_single.tadah", body);
    IgnoredKeyNoticeCapture cap;
    tadah::core::PeriodicTable::initialize();
    const Context c = CoreUserConfigParser().parse(single);
    DC_Selector s(c);
    CHECK(s.d2b->get_rcut() == Approx(6.0));
    CHECK(cap.notices() ==
          std::vector<std::string>{
              "Part of RCUT2B (5) is ignored: TYPE2B D2_Blip takes one RCUT2B "
              "value, and its cutoff is the largest, 6; the potential works "
              "unchanged with RCUT2B 6 alone"});
  }
}

/* A training configuration read from a file, without a potential: the
 * refusal of a count names the line that sets the key and its file (a line
 * of no file, its number).
 * Before the fix: 'TYPE2B names 1 descriptor (D2_Blip) but RCUT2B has 2
 * values: each descriptor ...' and 'TYPE2B D2_Blip takes 2 CGRID2B values,
 * but CGRID2B has 3 values: TYPE2B D2_Blip 2 2 Ta Ta', without the line. */
TEST_CASE("The refusal of a count of a training configuration names the "
          "lines of the key",
          "[D_mJoin][list_keys][potential_lines]") {
  Scratch dir;
  std::string body = kSingle2B;
  body.erase(body.find("WEIGHTS"));
  SECTION("RCUT2B") {
    std::string two = body;
    two.replace(two.find("RCUT2B 5\n"), 9, "RCUT2B 5 6\n");
    const std::string cfg = dir.write("two.train", two);
    CHECK(refusal(parsed(cfg)) ==
          "TYPE2B names 1 descriptor (D2_Blip) but RCUT2B has 2 values, set "
          "on line 6 of " + cfg + kOneRcut2B);
  }
  SECTION("CGRID2B, by lines of no file") {
    std::vector<std::string> lines;
    for (std::size_t at = 0, end; (end = body.find('\n', at)) !=
                                  std::string::npos;
         at = end + 1)
      lines.push_back(body.substr(at, end - at));
    lines.push_back("CGRID2B 4");
    tadah::core::PeriodicTable::initialize();
    CHECK(refusal(CoreUserConfigParser().parse(lines)) ==
          "TYPE2B D2_Blip takes 2 CGRID2B values, but CGRID2B has 3 values, "
          "set on line 7 and on line 9: TYPE2B D2_Blip 2 2 Ta Ta");
  }
}
