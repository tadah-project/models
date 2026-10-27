/* Keys that change nothing a fit gives (inert_keys.h). A fit -- tadah
 * train, tadah hpo, tadah refit and a program's own training -- refuses a
 * configuration that has a line of one, naming the key, the line and the
 * file:
 *  - a deprecated key (MPIWPCKG, MBLOCK, NBLOCK, SIGMA), whose lines the
 *    parse records (core::Context::deprecated_keys());
 *  - EMBF0 and RHOFORM when the configuration has no DM_REAM block;
 *  - WATOM with a value other than 1 when no descriptor of the
 *    configuration reads it (D2_EAM and DM_EAM take their functions from
 *    the SETFL file, unweighted).
 * A potential loaded to predict with gives a notice for each instead
 * (embf0.h, the notice of a deprecated key, and DC_Selector for WATOM).
 */

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/core_user_config_parser.h>
#include <tadah/core/logger.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/inert_keys.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>
#include <unistd.h>

using tadah::core::Context;
using tadah::core::CoreUserConfigParser;
using namespace tadah::models;

namespace {

/* A scratch directory unique to this process, removed with its files. */
class ScratchDir {
public:
  ScratchDir()
      : dir_(std::filesystem::temp_directory_path() /
             ("tadah_models_inert_keys_" + std::to_string(::getpid()))) {
    std::filesystem::create_directories(dir_);
  }
  ~ScratchDir() {
    std::error_code ec;
    std::filesystem::remove_all(dir_, ec);
  }
  ScratchDir(const ScratchDir &) = delete;
  ScratchDir &operator=(const ScratchDir &) = delete;
  std::string write_file(const std::string &name,
                         const std::string &body) const {
    const auto p = dir_ / name;
    std::ofstream(p) << body;
    return p.string();
  }

private:
  std::filesystem::path dir_;
};

/* Standard error and output silenced while a parse writes its notices. */
class Quiet {
public:
  Quiet()
      : out_(std::cout.rdbuf(sink_.rdbuf())),
        err_(std::cerr.rdbuf(sink_.rdbuf())) {}
  ~Quiet() {
    std::cout.rdbuf(out_);
    std::cerr.rdbuf(err_);
  }
  Quiet(const Quiet &) = delete;
  Quiet &operator=(const Quiet &) = delete;

private:
  std::ostringstream sink_;
  std::streambuf *out_, *err_;
};

/* A context built as a program builds one (no file, no lines). */
Context ctx_of(const std::vector<std::pair<std::string, std::string>> &kv) {
  Context c;
  for (const auto &p : kv) {
    c.remove(p.first);
    c.add(p.first, p.second);
  }
  return c;
}

const char *const kLj =
    "ATOM Ta\nINIT2B true\nINITMB false\nTYPE2B D2_LJ Ta Ta\n"
    "RCTYPE2B Cut_Cos\nRCUT2B 5.0\nMODEL M_BLR BF_Linear\n";
// D2_EAM and DM_EAM, without ATOM and WATOM (the cases give them).
const char *const kEam =
    "INIT2B true\nTYPE2B D2_EAM 0 Ta Ta\nINITMB true\n"
    "TYPEMB DM_EAM 0 Ta Ta\nSETFL eam.alloy\nRCUT2B 5.58190742024\n"
    "RCUTMB 5.58190742024\nMODEL M_BLR BF_Linear\n";

std::string refusal(const Context &c, const InertKeyRules &rules = {}) {
  try {
    require_no_inert_keys(c, "tadah train", rules);
  } catch (const std::runtime_error &e) {
    return e.what();
  }
  return "";
}

}  // namespace

/* Pins: context_reads_watom is true when an initialised TYPE2B or TYPEMB
 * line names a descriptor that weights its pairs by WATOM, or a component
 * (an @id of a DESCRIPTOR block), and false when every descriptor it names
 * is D2_EAM, DM_EAM or none: D2_EAM + DM_EAM, DM_EAM alone with INIT2B
 * false, INIT2B and INITMB false: false; D2_LJ, D2_EAM with a DM_mJoin of
 * DM_EAM and DM_Blip, a component TYPEMB @ta: true. A D2_LJ line with
 * INIT2B false is not read.
 * Before the fix: does not compile (no inert_keys.h). */
TEST_CASE("context_reads_watom: whether a descriptor of a configuration "
          "weights its pairs by WATOM",
          "[inert_keys][watom]") {
  CHECK_FALSE(context_reads_watom(ctx_of({{"INIT2B", "true"},
                                          {"TYPE2B", "D2_EAM 0 Ta Ta"},
                                          {"INITMB", "true"},
                                          {"TYPEMB", "DM_EAM 0 Ta Ta"}})));
  CHECK_FALSE(context_reads_watom(ctx_of({{"INIT2B", "false"},
                                          {"TYPE2B", "D2_LJ Ta Ta"},
                                          {"INITMB", "true"},
                                          {"TYPEMB", "DM_EAM 0 Ta Ta"}})));
  CHECK_FALSE(context_reads_watom(
      ctx_of({{"INIT2B", "false"}, {"INITMB", "false"}})));
  CHECK(context_reads_watom(
      ctx_of({{"INIT2B", "true"}, {"TYPE2B", "D2_LJ Ta Ta"}})));
  CHECK(context_reads_watom(ctx_of(
      {{"INIT2B", "true"},
       {"TYPE2B", "D2_EAM 0 Ta Ta"},
       {"INITMB", "true"},
       {"TYPEMB", "DM_mJoin DM_EAM 0 Ta Ta DM_Blip 0 1 1 Ta Ta"}})));
  CHECK(context_reads_watom(
      ctx_of({{"INIT2B", "false"}, {"INITMB", "true"}, {"TYPEMB", "@ta"}})));
}

/* Pins: require_no_inert_keys refuses a configuration with a line that
 * changes nothing a fit gives, naming who refuses, the key with its values
 * and the line and file of the line: MBLOCK 64, NBLOCK 8 and MPIWPCKG 50
 * (deprecated); EMBF0 FREE and RHOFORM SIGNED with no DM_REAM block;
 * WATOM 2.5, and WATOM Ta Ta 2.5 (a pair triple), with D2_EAM and DM_EAM ("WATOM
 * 2.5 (line 2 of eam.train)"). Several such lines give one refusal that
 * lists each. Accepted: the same configurations without those lines, WATOM
 * 1 with D2_EAM and DM_EAM, WATOM 2.5 with D2_LJ, and, with the rules of a
 * fit that reads EMBF0 and RHOFORM itself (InertKeyRules{false}, tadah
 * refit), EMBF0 FREE and RHOFORM SIGNED without a DM_REAM block.
 * Before the fix: does not compile (no inert_keys.h): a fit accepted each
 * line, the deprecated keys and EMBF0 with a WARNING, WATOM without one. */
TEST_CASE("require_no_inert_keys: a fit refuses a line that changes nothing "
          "it gives, naming the key, the line and the file",
          "[inert_keys]") {
  tadah::core::PeriodicTable::initialize();
  const ScratchDir dir;
  Quiet quiet;
  auto parse = [&](const std::string &name, const std::string &body) {
    return CoreUserConfigParser().parse(dir.write_file(name, body));
  };
  struct Case {
    std::string name, extra, named;
  };
  for (const Case &k : std::vector<Case>{
           {"mblock.train", "MBLOCK 64\n", "MBLOCK 64 (line 8 of "},
           {"nblock.train", "NBLOCK 8\n", "NBLOCK 8 (line 8 of "},
           {"mpiwpckg.train", "MPIWPCKG 50\n", "MPIWPCKG 50 (line 8 of "},
           {"embf0.train", "EMBF0 FREE\n", "EMBF0 FREE (line 8 of "},
           {"rhoform.train", "RHOFORM SIGNED\n",
            "RHOFORM SIGNED (line 8 of "}}) {
    INFO(k.name);
    const std::string path = dir.write_file(k.name, std::string(kLj) + k.extra);
    const std::string what = refusal(CoreUserConfigParser().parse(path));
    INFO("refusal: " << what);
    CHECK_THAT(what, Catch::StartsWith("tadah train: "));
    CHECK_THAT(what, Catch::Contains(k.named + path + ")"));
  }
  for (const std::string watom : {"2.5", "Ta Ta 2.5"}) {
    INFO("WATOM " << watom);
    const std::string path = dir.write_file(
        "eam.train", "ATOM Ta\nWATOM " + watom + "\n" + kEam);
    const Context c = CoreUserConfigParser().parse(path);
    CHECK_FALSE(context_reads_watom(c));
    const std::string what = refusal(c);
    INFO("refusal: " << what);
    CHECK_THAT(what, Catch::Contains("WATOM " + watom + " (line 2 of " +
                                     path + ")"));
    CHECK_THAT(what, Catch::Contains("D2_EAM"));
  }
  {
    const std::string path =
        dir.write_file("several.train", std::string(kLj) +
                                            "MBLOCK 64\nEMBF0 FREE\n");
    const std::string what = refusal(CoreUserConfigParser().parse(path));
    INFO("refusal: " << what);
    CHECK_THAT(what, Catch::Contains("MBLOCK 64 (line 8 of " + path + ")"));
    CHECK_THAT(what, Catch::Contains("EMBF0 FREE (line 9 of " + path + ")"));
  }
  CHECK(refusal(parse("lj.train", kLj)).empty());
  CHECK(refusal(parse("eam1.train", std::string("ATOM Ta\nWATOM 1\n") + kEam))
            .empty());
  CHECK(refusal(parse("lj_watom.train", std::string(kLj) + "WATOM 2.5\n"))
            .empty());
  {
    const Context c =
        parse("refit.train", std::string(kLj) + "EMBF0 FREE\nRHOFORM SIGNED\n");
    CHECK(refusal(c, InertKeyRules{false}).empty());
    CHECK_FALSE(refusal(c).empty());
  }
}
