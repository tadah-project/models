// The weights of a model and the descriptors it predicts from.
//
// Function_Base::weights_count(dsize) is the number of model weights (the
// design-matrix width) of each basis function and kernel for a descriptor
// of dsize entries, and basis_vector_size() the length of the basis
// vectors of a kernel with a basis. require_weights_fit() refuses a
// potential whose WEIGHTS, NMEAN/NSTDEV (NORM true) or basis vectors do not
// fit the descriptor of its configuration, before anything is predicted:
// the dot products of the predictors run over the weights, so a potential
// with weights for another descriptor layout predicted silently wrong
// energies and forces (a potential trained by tadah before 1.6.0 with
// BIAS true and two or more species has one intercept column shared by all
// species, while pair_style tadah builds one per species).

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/dc_selector.h>
#include <tadah/models/functions/f_all.h>
#include <tadah/models/functions/weights_fit.h>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

namespace {

/* BIAS and ATOM, NORM, and n WEIGHTS / NMEAN / NSTDEV values. */
Context pot(bool bias, const std::vector<std::string> &atoms, bool norm,
            std::size_t n_weights, std::size_t n_mean = 0,
            std::size_t n_stdev = 0) {
  Context c;
  c.add("BIAS", bias ? "true" : "false");
  for (const auto &a : atoms) c.add("ATOM", a);
  c.add("NORM", norm ? "true" : "false");
  for (std::size_t i = 0; i < n_weights; ++i) c.add("WEIGHTS", 0.1 * i);
  for (std::size_t i = 0; i < n_mean; ++i) c.add("NMEAN", 0.2 * i);
  for (std::size_t i = 0; i < n_stdev; ++i) c.add("NSTDEV", 1.0 + i);
  return c;
}

Matrix basis_of(std::size_t length, std::size_t count) {
  Matrix b(length, count);
  b.set_zero();
  for (std::size_t j = 0; j < count; ++j) b(j % length, j) = 1.0;
  return b;
}

std::string refusal(const Context &c, std::size_t dsize,
                    const Function_Base &f) {
  try {
    require_weights_fit(c, dsize, f);
  } catch (const std::runtime_error &e) {
    return e.what();
  }
  return "<accepted>";
}

const char *kLegacy = "one intercept column shared by all species";

}  // namespace

// Pins: weights_count(dsize) is dsize for BF_Linear and Kern_Linear,
// dsize (dsize + 1) / 2 for BF_Polynomial2, and the number of basis vectors
// for the kernels with a basis (Kern_RBF, Kern_Quadratic, Kern_Polynomial,
// Kern_LQ, Kern_Sigmoid), whatever dsize is; basis_vector_size() is the
// length of those basis vectors, and 0 for BF_Linear, BF_Polynomial2 and
// Kern_Linear (no basis).
// Before the fix: weights_count and basis_vector_size did not exist (run
// with stubs that throw std::logic_error: every section fails on that
// exception).
TEST_CASE("Function_Base::weights_count: the model weights of each basis "
          "function and kernel for a descriptor of dsize entries",
          "[weights_fit]") {
  SECTION("basis functions and the linear kernel") {
    CHECK(BF_Linear().weights_count(4) == 4);
    CHECK(BF_Linear().weights_count(1) == 1);
    CHECK(BF_Polynomial2().weights_count(4) == 10);
    CHECK(BF_Polynomial2().weights_count(3) == 6);
    CHECK(Kern_Linear().weights_count(4) == 4);
    CHECK(BF_Linear().basis_vector_size() == 0);
    CHECK(BF_Polynomial2().basis_vector_size() == 0);
    Kern_Linear kl;
    kl.set_basis(basis_of(4, 5));   // Kern_Linear has no basis
    CHECK(kl.basis_vector_size() == 0);
    CHECK(kl.weights_count(4) == 4);
  }
  SECTION("kernels with a basis") {
    std::vector<std::unique_ptr<Kern_Base>> kernels;
    kernels.emplace_back(new Kern_RBF());
    kernels.emplace_back(new Kern_Quadratic());
    kernels.emplace_back(new Kern_Polynomial());
    kernels.emplace_back(new Kern_LQ());
    kernels.emplace_back(new Kern_Sigmoid());
    for (auto &k : kernels) {
      CAPTURE(k->get_label());
      CHECK(k->weights_count(4) == 0);   // no basis yet
      CHECK(k->basis_vector_size() == 0);
      k->set_basis(basis_of(4, 5));
      CHECK(k->weights_count(4) == 5);
      CHECK(k->weights_count(9) == 5);
      CHECK(k->basis_vector_size() == 4);
    }
  }
}

// Pins: require_weights_fit accepts potentials whose weights fit the
// descriptor: BF_Linear with ATOM Nb Ta and BIAS true on DSIZE 4 with 4
// WEIGHTS (and with NORM true, 4 NMEAN and 4 NSTDEV); BIAS false on DSIZE 2
// with 2 WEIGHTS; BF_Polynomial2 on DSIZE 4 with 10 WEIGHTS; Kern_RBF with
// 5 basis vectors of length 4 and 5 WEIGHTS.
// Before the fix: require_weights_fit did not exist (run with a stub that
// throws std::logic_error: every check fails on that exception).
TEST_CASE("require_weights_fit accepts weights that fit the descriptor",
          "[weights_fit]") {
  PeriodicTable::initialize();
  CHECK(refusal(pot(true, {"Nb", "Ta"}, false, 4), 4, BF_Linear()) ==
        "<accepted>");
  CHECK(refusal(pot(true, {"Nb", "Ta"}, true, 4, 4, 4), 4, BF_Linear()) ==
        "<accepted>");
  CHECK(refusal(pot(false, {"Nb", "Ta"}, false, 2), 2, BF_Linear()) ==
        "<accepted>");
  CHECK(refusal(pot(true, {"Nb", "Ta"}, false, 10), 4, BF_Polynomial2()) ==
        "<accepted>");
  Kern_RBF rbf;
  rbf.set_basis(basis_of(4, 5));
  CHECK(refusal(pot(true, {"Nb", "Ta"}, false, 5), 4, rbf) == "<accepted>");
  // The overload with explicit sizes, as a trained model calls it.
  CHECK_NOTHROW(require_weights_fit(pot(true, {"Nb", "Ta"}, true, 0), 4,
                                    BF_Linear(), TrainedSizes{4, 4, 4}));
}

// Pins: require_weights_fit refuses, with std::runtime_error naming
// WEIGHTS, DSIZE, BIAS and ATOM:
//   - BF_Linear, ATOM Nb Ta, BIAS true, DSIZE 4 with 3 WEIGHTS (the layout
//     of a potential trained before 1.6.0: one intercept column shared by
//     the two species), with the hint that names that layout and says to
//     retrain; with 5 WEIGHTS without the hint;
//   - BIAS false, DSIZE 2 with 3 WEIGHTS ("BIAS false"), no hint;
//   - ATOM Ta (one species), DSIZE 3 with 2 WEIGHTS, no hint;
//   - BF_Polynomial2, ATOM Nb Ta, DSIZE 4: 6 WEIGHTS (the count for the
//     one-intercept layout, DSIZE 3) with the hint, 9 without it;
//   - the explicit-size overload with sizes that do not fit.
// Before the fix: require_weights_fit did not exist (run with a stub that
// throws std::logic_error: every check fails on the stub's message).
TEST_CASE("require_weights_fit refuses WEIGHTS that do not fit the "
          "descriptor, naming WEIGHTS, DSIZE, BIAS and ATOM",
          "[weights_fit]") {
  PeriodicTable::initialize();
  SECTION("BF_Linear, two species, the one-intercept layout") {
    const std::string m = refusal(pot(true, {"Nb", "Ta"}, false, 3), 4,
                                  BF_Linear());
    INFO(m);
    CHECK(m.find("WEIGHTS has 3 values") != std::string::npos);
    CHECK(m.find("DSIZE 4") != std::string::npos);
    CHECK(m.find("BIAS true") != std::string::npos);
    CHECK(m.find("ATOM Nb Ta") != std::string::npos);
    CHECK(m.find(kLegacy) != std::string::npos);
    CHECK(m.find("retrain") != std::string::npos);
  }
  SECTION("BF_Linear, two species, another count") {
    const std::string m = refusal(pot(true, {"Nb", "Ta"}, false, 5), 4,
                                  BF_Linear());
    INFO(m);
    CHECK(m.find("WEIGHTS has 5 values") != std::string::npos);
    CHECK(m.find("DSIZE 4") != std::string::npos);
    CHECK(m.find(kLegacy) == std::string::npos);
  }
  SECTION("BIAS false") {
    const std::string m = refusal(pot(false, {"Nb", "Ta"}, false, 3), 2,
                                  BF_Linear());
    INFO(m);
    CHECK(m.find("WEIGHTS has 3 values") != std::string::npos);
    CHECK(m.find("BIAS false") != std::string::npos);
    CHECK(m.find(kLegacy) == std::string::npos);
  }
  SECTION("one species") {
    const std::string m = refusal(pot(true, {"Ta"}, false, 2), 3,
                                  BF_Linear());
    INFO(m);
    CHECK(m.find("WEIGHTS has 2 values") != std::string::npos);
    CHECK(m.find("ATOM Ta") != std::string::npos);
    CHECK(m.find(kLegacy) == std::string::npos);
  }
  SECTION("BF_Polynomial2") {
    const std::string legacy = refusal(pot(true, {"Nb", "Ta"}, false, 6), 4,
                                       BF_Polynomial2());
    INFO(legacy);
    CHECK(legacy.find("WEIGHTS has 6 values") != std::string::npos);
    CHECK(legacy.find("needs 10") != std::string::npos);
    CHECK(legacy.find(kLegacy) != std::string::npos);
    const std::string other = refusal(pot(true, {"Nb", "Ta"}, false, 9), 4,
                                      BF_Polynomial2());
    INFO(other);
    CHECK(other.find("WEIGHTS has 9 values") != std::string::npos);
    CHECK(other.find(kLegacy) == std::string::npos);
  }
  SECTION("explicit sizes") {
    CHECK_THROWS_WITH(require_weights_fit(pot(true, {"Nb", "Ta"}, false, 4),
                                          4, BF_Linear(),
                                          TrainedSizes{3, 0, 0}),
                      Catch::Contains("WEIGHTS has 3 values") &&
                          Catch::Contains(kLegacy));
  }
}

// Pins: with NORM true, NMEAN and NSTDEV must have one value per
// descriptor entry: NMEAN 3 (the one-intercept layout of DSIZE 4 with two
// species) is refused naming NMEAN, with the hint; NSTDEV 5 naming NSTDEV,
// without it; a kernel whose basis vectors have 3 entries on DSIZE 4 (5
// basis vectors, 5 WEIGHTS) is refused naming BASIS, with the hint.
// Before the fix: require_weights_fit did not exist (run with a stub that
// throws std::logic_error: every check fails on the stub's message).
TEST_CASE("require_weights_fit refuses NMEAN, NSTDEV and basis vectors "
          "that do not fit the descriptor",
          "[weights_fit]") {
  PeriodicTable::initialize();
  SECTION("NMEAN") {
    const std::string m = refusal(pot(true, {"Nb", "Ta"}, true, 4, 3, 4), 4,
                                  BF_Linear());
    INFO(m);
    CHECK(m.find("NMEAN has 3 values") != std::string::npos);
    CHECK(m.find("DSIZE 4") != std::string::npos);
    CHECK(m.find(kLegacy) != std::string::npos);
  }
  SECTION("NSTDEV") {
    const std::string m = refusal(pot(true, {"Nb", "Ta"}, true, 4, 4, 5), 4,
                                  BF_Linear());
    INFO(m);
    CHECK(m.find("NSTDEV has 5 values") != std::string::npos);
    CHECK(m.find(kLegacy) == std::string::npos);
  }
  SECTION("NORM false: NMEAN and NSTDEV are not read") {
    CHECK(refusal(pot(true, {"Nb", "Ta"}, false, 4, 3, 3), 4, BF_Linear()) ==
          "<accepted>");
  }
  SECTION("basis vectors") {
    Kern_RBF rbf;
    rbf.set_basis(basis_of(3, 5));
    const std::string m = refusal(pot(true, {"Nb", "Ta"}, false, 5), 4, rbf);
    INFO(m);
    CHECK(m.find("BASIS") != std::string::npos);
    CHECK(m.find("3 entries") != std::string::npos);
    CHECK(m.find("DSIZE 4") != std::string::npos);
    CHECK(m.find(kLegacy) != std::string::npos);
  }
  SECTION("WEIGHTS of a kernel: one per basis vector") {
    Kern_RBF rbf;
    rbf.set_basis(basis_of(4, 5));
    const std::string m = refusal(pot(true, {"Nb", "Ta"}, false, 4), 4, rbf);
    INFO(m);
    CHECK(m.find("WEIGHTS has 4 values") != std::string::npos);
    CHECK(m.find("needs 5") != std::string::npos);
  }
}

// Pins: DC_Selector, the descriptor builder of pair_style tadah, gives
// BIAS true one intercept column per ATOM species (ATOM Nb Ta: bias 2, the
// two-body block from column 2, DSIZE 2 + its size) and refuses BIAS true
// without ATOM, naming ATOM.
// Before the fix: passes as well (a pin: the refusal is CORE's
// n_intercept_columns, CORE ad5b2ec; before that DC_Selector built no
// intercept column for BIAS true without ATOM, while MLIP's DescriptorsCalc
// built one).
TEST_CASE("DC_Selector: one intercept column per ATOM species; BIAS true "
          "without ATOM is refused",
          "[weights_fit][dc_selector]") {
  PeriodicTable::initialize();
  Context c;
  c.add("INIT2B", "true");
  c.add("INITMB", "false");
  c.add("TYPE2B", "D2_LJ");
  c.add("TYPE2B", "*");
  c.add("TYPE2B", "*");
  c.add("RCTYPE2B", "Cut_Cos");
  c.add("RCUT2B", 5.0);
  c.add("BIAS", "true");
  Context no_atom = c;
  c.add("ATOM", "Nb");
  c.add("ATOM", "Ta");
  DC_Selector dcs(c);
  CHECK(dcs.bias == 2);
  CHECK(dcs.d2b->get_fidx() == 2);
  CHECK(dcs.descriptor_size() == 2 + dcs.d2b->size());
  CHECK_THROWS_WITH(DC_Selector(no_atom), Catch::Contains("ATOM"));
}
