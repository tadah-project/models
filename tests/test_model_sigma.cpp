// test_model_sigma.cpp
//
// The weight covariance belongs to the fit that computed it. An OALGO 3 fit
// (LinearRegressor::train, the SVD solve) returns it, and M_BLR_Train /
// M_KRR_Train hold it in their Sigma; a fit that computes none (another
// OALGO, M_KRR implementation 2) leaves the model without one, whatever an
// earlier fit of the model computed. LinearRegressor::train writes nothing
// into the Context: the potential gets the covariance from the model
// (MLIP's get_param_file, keys WCOV and WNOISE). Until WCOV replaced it,
// the fit wrote it into the Context as SIGMA (config_add_sigma), and
// read_sigma read it back; the count checks of read_sigma are those of
// the WCOV reader now (MLIP, test_predict_error.cpp and
// test_weight_covariance_io.cpp).
//
// Also: an M_KRR model built from a Context with WEIGHTS (a potential used
// as the configuration) can be trained, as an M_BLR model can.

#include "catch2/catch.hpp"
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/models/functions/basis_functions/bf_linear.h>
#include <tadah/models/functions/kernels/kern_linear.h>
#include <tadah/models/functions/kernels/kern_rbf.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/m_blr_train.h>
#include <tadah/models/m_krr_train.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

namespace {

std::uint64_t bits(double v) {
  std::uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return u;
}

// Deterministic, well-conditioned 8x4 design matrix (as test_fixed_weights).
phi_type make_phi() {
  phi_type Phi(8, 4);
  double v = 0.17;
  for (size_t j = 0; j < Phi.cols(); ++j)
    for (size_t i = 0; i < Phi.rows(); ++i) {
      v = std::fmod(v * 1.618 + 0.31 * static_cast<double>(j + 1), 1.0) - 0.5;
      Phi(i, j) = v + (i == j ? 2.0 : 0.0);
    }
  return Phi;
}

t_type make_target(const phi_type &Phi) {
  t_type w_true(Phi.cols());
  for (size_t j = 0; j < Phi.cols(); ++j)
    w_true[j] = 0.5 * static_cast<double>(j + 1) * (j % 2 ? -1.0 : 1.0);
  t_type T = Phi * w_true;
  for (size_t i = 0; i < T.size(); ++i) T[i] += 0.01 * static_cast<double>(i % 3);
  return T;
}

/* make_phi() as the DesignMatrix gives it to a fit: its 8 rows, then one
 * ridge-augmentation row per column (zero for OALGO 3), so that the 8 rows
 * are the observations. */
phi_type make_phi_aug() {
  const phi_type Phi8 = make_phi();
  phi_type Phi(Phi8.rows() + Phi8.cols(), Phi8.cols());
  Phi.set_zero();
  for (size_t j = 0; j < Phi8.cols(); ++j)
    for (size_t i = 0; i < Phi8.rows(); ++i) Phi(i, j) = Phi8(i, j);
  return Phi;
}

t_type make_target_aug() {
  const t_type T8 = make_target(make_phi());
  t_type T(T8.size() + 4);
  for (size_t i = 0; i < T.size(); ++i) T[i] = i < T8.size() ? T8[i] : 0.0;
  return T;
}

Context make_ctx(const std::string &oalgo) {
  Context c;
  c.add("VERBOSE", "0");
  c.add("LAMBDA", "0");
  c.add("OALGO", oalgo);
  return c;
}

/* A SIGMA as an earlier OALGO 3 fit of a 4-weight model leaves it: N = 4
 * and 16 values. */
void add_earlier_sigma(Context &c) {
  c.add("SIGMA", "4");
  for (int k = 0; k < 16; ++k) c.add("SIGMA", 0.25 * (k + 1));
}

}  // namespace

/* Pins: two OALGO 3 fits (LinearRegressor::train) on one Context (that of
 * a model trained twice, or a potential used as the configuration) write
 * nothing into it: no SIGMA, no WCOV, no WNOISE; each returns its
 * covariance in its Sigma out-parameter, the second fit's equal bit for
 * bit to the first's (4 x 4).
 * Before the fix (WCOV and WNOISE): each fit wrote SIGMA into the Context,
 * N = 4 and the 16 values (config_add_sigma replaced the earlier SIGMA;
 * before fix/hpo-trainer-reuse it appended to it, 34 values after the
 * second fit). */
TEST_CASE("OALGO 3: LinearRegressor::train writes no covariance into the "
          "Context",
          "[sigma]") {
  Context c = make_ctx("3");
  const phi_type Phi0 = make_phi_aug();
  const t_type T0 = make_target_aug();
  Matrix Sigma[2];
  for (int fit = 0; fit < 2; ++fit) {
    phi_type Phi = Phi0;
    t_type T = T0;
    t_type w;
    LinearRegressor::train(c, Phi, T, w, Sigma[fit]);
    INFO("fit " << fit);
    CHECK_FALSE(c.has_stored("SIGMA"));
    CHECK_FALSE(c.has_stored("WCOV"));
    CHECK_FALSE(c.has_stored("WNOISE"));
  }
  REQUIRE(Sigma[0].rows() == 4);
  REQUIRE(Sigma[1].rows() == 4);
  REQUIRE(Sigma[1].cols() == 4);
  for (size_t col = 0; col < 4; ++col)
    for (size_t row = 0; row < 4; ++row)
      CHECK(bits(Sigma[1](row, col)) == bits(Sigma[0](row, col)));
}

/* Pins: M_BLR_Train<BF_Linear> and M_KRR_Train<Kern_Linear> hold the
 * covariance of the fit that computed it and none otherwise, whatever
 * SIGMA the Context carried from an earlier fit (the 17 values of a
 * 4-weight OALGO 3 fit): with OALGO 1 (no FIXINDEX, FIXINDEX 2, FIXINDEX
 * 1-4) Sigma is empty and get_weights_uncertainty() throws naming OALGO
 * 3; with OALGO 3 Sigma is the 4 x 4 covariance of every weight and
 * get_weights_uncertainty() gives 4 standard deviations, also with
 * FIXINDEX 2 (M_Train::sigma_of_all_weights; test_fixed_weights.cpp and
 * test_weight_covariance.cpp check the values) and with FIXINDEX 1-4,
 * every weight pinned and nothing solved (a zero covariance); with every
 * weight pinned the weights are the FIXWEIGHT values. The Context's SIGMA
 * is not read.
 * Before the fix (WCOV and WNOISE): get_weights_uncertainty() tested
 * LAMBDA (0 here), not the fit: it threw "Sigma matrix is not computed."
 * after every fit (M_KRR: "Sigma matrix is only computed for LAMBDA !=
 * 0"), and with every weight pinned an OALGO 3 fit left Sigma empty.
 * (The fit's covariance was the Context's SIGMA, which a fit replaced or
 * removed: the expectation of this case, then named "A fit replaces or
 * removes the SIGMA of an earlier fit", before WCOV.) */
TEST_CASE("A model holds the covariance of its last fit only",
          "[sigma]") {
  struct Pins {
    const char *name;
    const char *fixindex;   // nullptr: no FIXINDEX
    std::vector<double> fixweight;
  };
  const Pins pins[] = {{"no FIXINDEX", nullptr, {}},
                       {"FIXINDEX 2", "2", {0.7}},
                       {"FIXINDEX 1-4", "1-4", {0.1, 0.2, 0.3, 0.4}}};
  for (const char *oalgo : {"1", "3"}) {
    for (const Pins &p : pins) {
      for (const char *model : {"M_BLR BF_Linear", "M_KRR Kern_Linear"}) {
        const bool all_pinned = p.fixindex && std::string(p.fixindex) == "1-4";
        DYNAMIC_SECTION(model << ", OALGO " << oalgo << ", " << p.name) {
          Context c = make_ctx(oalgo);
          if (p.fixindex) {
            c.add("FIXINDEX", p.fixindex);
            for (double w : p.fixweight) c.add("FIXWEIGHT", w);
          }
          add_earlier_sigma(c);
          phi_type Phi = make_phi_aug();
          t_type T = make_target_aug();
          t_type w;
          Matrix S;
          t_type unc;
          std::string msg;
          auto fit = [&](auto &m) {
            REQUIRE_NOTHROW(m.train(Phi, T));
            w = m.get_weights();
            S = m.Sigma;
            try {
              unc = m.get_weights_uncertainty();
            } catch (const std::exception &e) {
              msg = e.what();
            }
          };
          if (std::string(model) == "M_BLR BF_Linear") {
            BF_Linear bf;
            M_BLR_Train<BF_Linear> m(bf, c);
            fit(m);
          } else {
            Kern_Linear k;
            M_KRR_Train<Kern_Linear> m(k, c);
            fit(m);
          }
          INFO("get_weights_uncertainty: " << msg);
          REQUIRE(w.size() == 4);
          if (all_pinned)
            for (size_t i = 0; i < 4; ++i) CHECK(w[i] == p.fixweight[i]);
          if (std::string(oalgo) == "3") {
            // The covariance of this fit, of every weight (a pinned one
            // with no variance).
            CHECK(S.rows() == 4);
            CHECK(S.cols() == 4);
            CHECK(unc.size() == 4);
            CHECK(msg.empty());
          } else {
            CHECK(S.rows() == 0);
            CHECK(unc.size() == 0);
            CHECK(msg.find("OALGO 3") != std::string::npos);
          }
        }
      }
    }
  }
}

/* Pins: M_KRR_Train::train2, the solve of MODEL M_KRR <kernel> 2 (the
 * kernel matrix of the basis vectors plus LAMBDA on its diagonal, no
 * covariance), leaves the model without a covariance after an OALGO 3 fit
 * of the same model computed one (Sigma 4 x 4 after train(), empty after
 * train2(), and get_weights_uncertainty() throws naming OALGO 3), for
 * Kern_Linear and for Kern_RBF (MPARAMS 1.0); it fits one weight per basis
 * vector (4).
 * Before the fix (WCOV and WNOISE): Sigma kept the 4 x 4 covariance of the
 * earlier fit beside train2's weights, and get_weights_uncertainty()
 * returned its diagonal (LAMBDA 1e-3). (train2 removed the Context's
 * SIGMA: the expectation of this case, then named "M_KRR implementation 2
 * (train2) removes the SIGMA of an earlier fit", before WCOV.) */
TEST_CASE("M_KRR implementation 2 (train2) leaves no covariance of an "
          "earlier fit",
          "[sigma][krr]") {
  for (const char *kern : {"Kern_Linear", "Kern_RBF"}) {
    DYNAMIC_SECTION(kern) {
      Context c = make_ctx("3");
      c.remove("LAMBDA");
      c.add("LAMBDA", "1e-3");
      c.add("MPARAMS", "1.0");
      // Four basis vectors (the columns of the 8x4 matrix), one target
      // value each.
      Matrix B = make_phi();
      t_type T(B.cols());
      for (size_t j = 0; j < T.size(); ++j)
        T[j] = 0.3 * static_cast<double>(j + 1);
      t_type w;
      Matrix S_fit, S_train2;
      t_type unc;
      std::string msg;
      auto fit = [&](auto &m) {
        phi_type Phi = make_phi_aug();
        t_type Tf = make_target_aug();
        REQUIRE_NOTHROW(m.train(Phi, Tf));
        S_fit = m.Sigma;
        REQUIRE_NOTHROW(m.train2(B, T));
        w = m.get_weights();
        S_train2 = m.Sigma;
        try {
          unc = m.get_weights_uncertainty();
        } catch (const std::exception &e) {
          msg = e.what();
        }
      };
      if (std::string(kern) == "Kern_Linear") {
        Kern_Linear k;
        M_KRR_Train<Kern_Linear> m(k, c);
        fit(m);
      } else {
        Kern_RBF k(c);
        M_KRR_Train<Kern_RBF> m(k, c);
        fit(m);
      }
      INFO("get_weights_uncertainty: " << msg);
      CHECK(w.size() == 4);
      CHECK(S_fit.rows() == 4);
      CHECK(S_train2.rows() == 0);
      CHECK(unc.size() == 0);
      CHECK(msg.find("OALGO 3") != std::string::npos);
    }
  }
}

/* Pins: an M_KRR model built from a Context that has WEIGHTS (a potential
 * used as the configuration, `tadah train -c pot.tadah`, or an HPO -c
 * potential) holds those weights but is not trained, and train() fits it:
 * the weights are those of the same fit on a Context without WEIGHTS, bit
 * for bit. The same for M_BLR (which already behaved so).
 * Before the fix: M_KRR_Core::init set trained when WEIGHTS were present,
 * so train() threw "This object is already trained!". */
TEST_CASE("An M_KRR model built with WEIGHTS can be trained",
          "[krr][weights]") {
  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);
  t_type w_ref;
  {
    Context c = make_ctx("1");
    Kern_Linear k;
    M_KRR_Train<Kern_Linear> m(k, c);
    phi_type Phi = Phi0;
    t_type T = T0;
    m.train(Phi, T);
    w_ref = m.get_weights();
  }
  REQUIRE(w_ref.size() == 4);
  for (const char *model : {"M_KRR Kern_Linear", "M_BLR BF_Linear"}) {
    DYNAMIC_SECTION(model) {
      Context c = make_ctx("1");
      c.add("WEIGHTS", "9 8 7 6");
      phi_type Phi = Phi0;
      t_type T = T0;
      t_type w;
      bool trained_before = true;
      std::string msg;
      auto fit = [&](auto &m) {
        trained_before = m.is_trained();
        CHECK(m.get_weights().size() == 4);
        try {
          m.train(Phi, T);
        } catch (const std::exception &e) {
          msg = e.what();
        }
        w = m.get_weights();
      };
      if (std::string(model) == "M_KRR Kern_Linear") {
        Kern_Linear k;
        M_KRR_Train<Kern_Linear> m(k, c);
        fit(m);
      } else {
        BF_Linear bf;
        M_BLR_Train<BF_Linear> m(bf, c);
        fit(m);
      }
      INFO("train: " << msg);
      CHECK_FALSE(trained_before);
      CHECK(msg.empty());
      REQUIRE(w.size() == 4);
      for (size_t i = 0; i < 4; ++i) CHECK(bits(w[i]) == bits(w_ref[i]));
    }
  }
}
