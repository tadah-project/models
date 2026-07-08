#include <tadah/models/eam/eam_refit.h>
#include <tadah/models/eam/eam_lattice.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <tadah/models/descriptors/basis_dispatch.h>

#include <algorithm>
#include <cmath>
#include <functional>

namespace tadah {
namespace models {

EAMRefit::EAMRefit(const EAMSetflIO &eam, const EAMRefitConfig &cfg)
  : eam_(eam), cfg_(cfg) {
  // Shared, single-source-of-truth parser (Blip/Blip5/Gaussian/BP/Knot5 and the
  // Knot3/KnotMix/Mixed/M6 aliases). The Knot5 family's per-atom degree arrays
  // then select cubic vs quintic knot atoms.
  family_ = basis_kind_from_string(cfg.basis);
  rcut_   = eam_.rcut();
  r_lo_   = std::max(cfg_.r_lo_frac * rcut_, eam_.dr());
  rhomax_ = eam_.rhomax();
  rho_join_ = rhomax_;  // dense default: F_setfl == F_fit on the table range
}

// Hard cutoff (Cut_Dummy): f_c=1 inside r_cut, 0 beyond — matches LAMMPS
// eam/alloy so the native fit reproduces the original tabulated potential.
double EAMRefit::cutoff(double r) const {
  return (r <= rcut_) ? 1.0 : 0.0;
}

double EAMRefit::basisval(double x, double center, double width,
                          bool emb, int deg) const {
  // The SAME shared dispatcher the DM_REAM / D2 descriptors use, so a refitted
  // potential evaluates the identical basis at fit time and at run time (no
  // parallel hand-copy of the basis math to drift). `width` is the eta gauge
  // (1.0 for the knot family, whose knot arms ignore it); `emb` picks the
  // RIGHT-active embedding orientation vs the LEFT-active phi/rho orientation.
  return emb ? embedding_value(family_, x, width, center, deg)
             : radial_value_bare(family_, x, width, center, deg);
}

void EAMRefit::make_grid(double lo, double hi, int n,
                         tadah::core::v_type &centers,
                         tadah::core::v_type &widths) const {
  centers.assign(n, 0.0);
  widths.assign(n, 0.0);
  double dx = (n > 1) ? (hi - lo) / (n - 1) : (hi - lo);
  if (dx <= 0.0) dx = (hi > 0.0 ? hi : 1.0);
  for (int k = 0; k < n; ++k) {
    centers[k] = lo + dx * k;
    // Blip: eta = 1/dx  -> support |r-rs| < 2*dx, i.e. each blip overlaps its
    //   two nearest neighbours on each side (cubic-B-spline-like coverage).
    // Gaussian: eta = 0.5/dx^2 (sigma ~ dx).
    // Knot5: eta is a pure gauge — pinned to 1 (the knot mu is the shape).
    widths[k] = (family_ == BasisKind::Gaussian) ? 0.5 / (dx * dx)
              : (family_ == BasisKind::Knot5)    ? 1.0
                                              : 1.0 / dx;  // Blip + Blip5
  }
}

// Two uniform sub-grids: the bulk of the basis on the physically REACHABLE
// window [win_lo,win_hi] and a small uniform allocation on the unreachable
// tail(s), so the tabulated output stays defined everywhere.
//
// Why uniform-per-region and not curvature-equidistributed: nonuniformly
// stretched blips are NOT a B-spline basis on nonuniform knots — irregular
// spacing destroys the approximation power of uniform translates. Measured
// (Ti, accurate): curvature placement gave 49 meV/atom ambient error where
// uniform gave 0.35 meV. Region-wise uniform keeps uniform-grid quality while
// concentrating capacity where a crystal can actually go.
void EAMRefit::make_grid_windowed(double lo, double hi, int n,
                                  double win_lo, double win_hi,
                                  tadah::core::v_type &centers,
                                  tadah::core::v_type &widths) const {
  win_lo = std::max(win_lo, lo);
  win_hi = std::min(win_hi, hi);
  const double span = hi - lo, wspan = win_hi - win_lo;
  if (n < 8 || wspan <= 0.0 || wspan >= 0.999 * span) {
    make_grid(lo, hi, n, centers, widths);
    return;
  }
  // tail blip count: proportional to tail length but capped so the window
  // keeps >= 80% of the basis; at least 4 blips per non-empty tail
  const double tail_lo_len = win_lo - lo, tail_hi_len = hi - win_hi;
  auto tail_n = [&](double len) {
    if (len <= 0.0) return 0;
    return std::max(4, static_cast<int>(std::floor(n * std::min(0.2, len / span))));
  };
  int n_lo = tail_n(tail_lo_len), n_hi = tail_n(tail_hi_len);
  if (n_lo + n_hi > n / 3) {  // keep the window dominant
    const double f = (n / 3.0) / (n_lo + n_hi);
    n_lo = static_cast<int>(n_lo * f);
    n_hi = static_cast<int>(n_hi * f);
  }
  const int n_win = n - n_lo - n_hi;

  centers.clear(); widths.clear();
  centers.reserve(n); widths.reserve(n);
  auto add_uniform = [&](double a, double b, int m, bool incl_b) {
    if (m <= 0) return;
    const double dx = (m > 1) ? (b - a) / (incl_b ? (m - 1) : m)
                              : std::max(b - a, 1e-12);
    for (int k = 0; k < m; ++k) {
      centers.push_back(a + dx * k);
      widths.push_back((family_ == BasisKind::Gaussian) ? 0.5 / (dx * dx)
                     : (family_ == BasisKind::Knot5)    ? 1.0
                                                     : 1.0 / dx);
    }
  };
  add_uniform(lo, win_lo, n_lo, false);             // lower tail (excl. win_lo)
  add_uniform(win_lo, win_hi, n_win, n_hi == 0);    // the window
  if (n_hi > 0) {
    const double dx = (win_hi - win_lo) / std::max(n_win - 1, 1);
    add_uniform(win_hi + 0.0 * dx, hi, n_hi, true); // upper tail
  }
}

FitStats EAMRefit::lsq_fit(const std::vector<double> &xs,
                           const std::vector<double> &ys,
                           const std::vector<double> &wts,
                           const tadah::core::v_type &centers,
                           const tadah::core::v_type &widths,
                           bool apply_cutoff,
                           tadah::core::v_type &coeffs) const {
  const std::size_t Nd = xs.size();
  const std::size_t Nb = centers.size();

  // Allocate Nd + Nb rows: the extra Nb rows stay zero for LAMBDA=0 (OLS) and
  // receive the sqrt(lambda) ridge augmentation when LAMBDA>0 (see
  // LinearRegressor::train).
  tadah::core::phi_type Phi(Nd + Nb, Nb);
  Phi.set_zero();
  tadah::core::t_type T(Nd + Nb);
  for (std::size_t i = 0; i < Nd; ++i) {
    const double fc = apply_cutoff ? cutoff(xs[i]) : 1.0;
    const double sw = wts.empty() ? 1.0 : std::sqrt(wts[i]);
    T(i) = ys[i] * sw;
    for (std::size_t k = 0; k < Nb; ++k)
      Phi(i, k) = basisval(xs[i], centers[k], widths[k]) * fc * sw;
  }

  tadah::core::Context lc;
  lc.add("OALGO", 2);                 // GELS, default rcond 1e-8
  lc.add("LAMBDA", cfg_.lambda);
  tadah::core::t_type w;
  tadah::core::Matrix Sigma;
  LinearRegressor::train(lc, Phi, T, w, Sigma);

  coeffs.assign(Nb, 0.0);
  for (std::size_t k = 0; k < Nb; ++k) coeffs[k] = w(k);

  // Stats over the RELEVANT window (rows with weight >= 0.5): out-of-window
  // regions are deliberately fit crudely and would drown the meaningful
  // numbers. Falls back to all rows when unweighted.
  FitStats st;
  double se = 0.0, ae = 0.0, mx = 0.0;
  std::size_t nrel = 0;
  for (std::size_t i = 0; i < Nd; ++i) {
    if (!wts.empty() && wts[i] < 0.5) continue;
    const double fc = apply_cutoff ? cutoff(xs[i]) : 1.0;
    double pred = 0.0;
    for (std::size_t k = 0; k < Nb; ++k)
      pred += coeffs[k] * basisval(xs[i], centers[k], widths[k]) * fc;
    const double res = pred - ys[i];
    se += res * res;
    ae += std::fabs(res);
    mx = std::max(mx, std::fabs(res));
    ++nrel;
  }
  st.n = nrel;
  if (nrel > 0) { st.rmse = std::sqrt(se / nrel); st.mae = ae / nrel; }
  st.maxabs = mx;
  return st;
}

FitStats EAMRefit::lsq_fit_windowed(const std::vector<double> &xs,
                                    const std::vector<double> &ys,
                                    double win_lo, double win_hi,
                                    const tadah::core::v_type &centers,
                                    const tadah::core::v_type &widths,
                                    bool apply_cutoff,
                                    tadah::core::v_type &coeffs,
                                    bool pin, double pin_x, double pin_y) const {
  const std::size_t Nd = xs.size();
  const std::size_t Nb = centers.size();
  coeffs.assign(Nb, 0.0);

  // Partition blips by support overlap with the window (blip support: |x-c| <
  // 2/eta for blips, ~3 sigma for gaussians) and samples by location.
  // NOTE the dense windowed solve is NOT used for the one-sided Knot5 family
  // (its support is one-sided and unbounded toward the active side, so the
  // two-stage block separation does not apply; the engine restricts Knot5 to
  // compact mode, which injects externally-optimised parameters instead).
  std::vector<std::size_t> bA, bB, sA, sB;
  for (std::size_t k = 0; k < Nb; ++k) {
    const double half = (family_ == BasisKind::Gaussian)
                            ? 3.0 / std::sqrt(2.0 * widths[k])
                        : (family_ == BasisKind::Blip5)
                            ? 3.0 / widths[k]
                            : 2.0 / widths[k];
    if (centers[k] + half >= win_lo && centers[k] - half <= win_hi)
      bA.push_back(k);
    else
      bB.push_back(k);
  }
  for (std::size_t i = 0; i < Nd; ++i)
    (xs[i] >= win_lo && xs[i] <= win_hi ? sA : sB).push_back(i);

  auto solve_block = [&](const std::vector<std::size_t> &rows,
                         const std::vector<std::size_t> &cols,
                         bool with_pin) {
    if (cols.empty()) return;
    const std::size_t extra = with_pin ? 1 : 0;
    if (rows.size() + extra == 0) return;
    tadah::core::phi_type Phi(rows.size() + extra + cols.size(), cols.size());
    Phi.set_zero();
    tadah::core::t_type T(rows.size() + extra + cols.size());
    for (std::size_t r = 0; r < rows.size(); ++r) {
      const std::size_t i = rows[r];
      const double fc = apply_cutoff ? cutoff(xs[i]) : 1.0;
      // subtract the contribution of already-solved (other-block) blips
      double held = 0.0;
      for (std::size_t k = 0; k < Nb; ++k)
        if (coeffs[k] != 0.0)
          held += coeffs[k] * basisval(xs[i], centers[k], widths[k]) * fc;
      T(r) = ys[i] - held;
      for (std::size_t q = 0; q < cols.size(); ++q)
        Phi(r, q) = basisval(xs[i], centers[cols[q]], widths[cols[q]]) * fc;
    }
    if (with_pin) {  // heavily-weighted single constraint row
      const double sw = 1e4;
      const std::size_t r = rows.size();
      T(r) = pin_y * sw;
      for (std::size_t q = 0; q < cols.size(); ++q)
        Phi(r, q) = basisval(pin_x, centers[cols[q]], widths[cols[q]]) * sw;
    }
    tadah::core::Context lc;
    // GELSD (SVD, rank-revealing): a tail block can be legitimately
    // rank-deficient (a blip whose support holds no samples); GELS (QR)
    // aborts on that (DGELS info=1).
    lc.add("OALGO", 1);
    lc.add("LAMBDA", cfg_.lambda);
    tadah::core::t_type w;
    tadah::core::Matrix Sigma;
    LinearRegressor::train(lc, Phi, T, w, Sigma);
    for (std::size_t q = 0; q < cols.size(); ++q) coeffs[cols[q]] = w(q);
  };

  const bool pinA = pin && pin_x >= win_lo && pin_x <= win_hi;
  solve_block(sA, bA, pinA);
  solve_block(sB, bB, pin && !pinA);

  // stats over the WINDOW samples (the physically relevant region)
  FitStats st;
  double se = 0.0, ae = 0.0, mx = 0.0;
  for (std::size_t i : sA) {
    const double fc = apply_cutoff ? cutoff(xs[i]) : 1.0;
    double pred = 0.0;
    for (std::size_t k = 0; k < Nb; ++k)
      pred += coeffs[k] * basisval(xs[i], centers[k], widths[k]) * fc;
    const double res = pred - ys[i];
    se += res * res; ae += std::fabs(res); mx = std::max(mx, std::fabs(res));
  }
  st.n = sA.size();
  if (!sA.empty()) { st.rmse = std::sqrt(se / sA.size()); st.mae = ae / sA.size(); }
  st.maxabs = mx;
  return st;
}

// Reachability analysis on the native lattice (header a0/lattice): what pair
// distances / embedding densities can a crystal actually visit between strong
// compression (a = reach_smin*a0, V/V0 ~ smin^3) and mild expansion? Basis
// capacity concentrates there. Called by fit(); compact mode calls it
// directly because it needs the windows BEFORE any fitting.
void EAMRefit::analyze_reachability() {
  lattice_known_ = eam_lattice::supported(eam_.lattice());
  const std::string lat = lattice_known_ ? eam_.lattice() : "bcc";
  const double a0 = eam_.lattice_param();
  r_reach_ = 0.0; rho_reach_ = 0.0;
  if (a0 > 0.0) {
    const double r1nn = eam_lattice::nn_distance(lat, a0);
    r_reach_ = cfg_.reach_smin * r1nn;
    // densest reachable embedding density over the compression sweep
    const auto dists = eam_lattice::neighbour_distances(
        lat, a0, rcut_ / std::max(cfg_.reach_smin, 0.1));
    double rho_top = 0.0;
    for (int is = 0; is <= 28; ++is) {
      const double s = cfg_.reach_smin + (1.10 - cfg_.reach_smin) * is / 28.0;
      double acc = 0.0;
      for (double d : dists) {
        const double rs = s * d;
        if (rs < rcut_) acc += eam_.rho_tab(rs);
      }
      rho_top = std::max(rho_top, acc);
    }
    rho_reach_ = std::min(cfg_.reach_margin * rho_top, rhomax_);
    if (rho_reach_ <= 0.0) rho_reach_ = rhomax_;
  }
}

void EAMRefit::fit() {
  const bool adaptive = (cfg_.placement != "uniform");
  analyze_reachability();
  if (!adaptive) { r_reach_ = 0.0; rho_reach_ = 0.0; }  // legacy uniform mode

  // ---- radial dense sample grid (shared by phi and rho) ----
  const int nr_samp = std::max(
      {2048, cfg_.n2b * cfg_.oversample, cfg_.nrho * cfg_.oversample});
  std::vector<double> xr, yphi, yrho;
  xr.reserve(nr_samp); yphi.reserve(nr_samp); yrho.reserve(nr_samp);
  for (int i = 0; i < nr_samp; ++i) {
    const double r = r_lo_ + (rcut_ - r_lo_) * i / (nr_samp - 1);
    if (r <= 0.0) continue;
    xr.push_back(r);
    yphi.push_back(eam_.phi_tab(r));
    yrho.push_back(eam_.rho_tab(r));
  }

  // Radial fit window. The window's lower edge introduces a boundary
  // error-decay zone of a few blip spacings (measured phi err 5.8e-3 at the
  // edge vs 4.5e-5 interior), and that zone is FATAL when it touches visited
  // shells of a material with steep F' (Ti: drho ~ 0.01 x F'=+14.7 -> 172
  // meV/atom). So the window shrinks below the full range ONLY when a genuine
  // unfittable wall is detected (|phi| or rho exceeding K x its reachable-
  // window scale): tame files keep the full uniform grid and no interior
  // boundary at all.
  double rwin_lo = r_lo_;
  if (adaptive && r_reach_ > 0.0) {
    const double Kfit = 30.0;
    auto p2p = [&](const std::vector<double> &ys, double xmin) {
      double lo = 1e300, hi = -1e300;
      for (std::size_t i = 0; i < xr.size(); ++i)
        if (xr[i] >= xmin) { lo = std::min(lo, ys[i]); hi = std::max(hi, ys[i]); }
      return std::max(hi - lo, 1e-300);
    };
    const double phi_sc = p2p(yphi, r_reach_), rho_sc = p2p(yrho, r_reach_);
    double r_K = r_lo_;
    for (std::size_t i = 0; i < xr.size(); ++i) {
      if (std::fabs(yphi[i]) > Kfit * phi_sc || yrho[i] > Kfit * rho_sc)
        r_K = xr[i];  // xr ascending: ends at the LARGEST offending r
    }
    // never let the wall edge encroach on the reachable zone
    rwin_lo = std::min(std::max(r_lo_, r_K), 0.8 * r_reach_);
  }

  // ---- 1. pair phi(r) -> 2-body linear weights ----
  if (adaptive) make_grid_windowed(r_lo_, rcut_, cfg_.n2b, rwin_lo, rcut_,
                                   cgrid2b_, sgrid2b_);
  else          make_grid(r_lo_, rcut_, cfg_.n2b, cgrid2b_, sgrid2b_);
  phi_stats_ = lsq_fit_windowed(xr, yphi, rwin_lo, rcut_, cgrid2b_, sgrid2b_,
                                true, c2b_, false, 0.0, 0.0);

  // ---- 2. density rho(r) -> baked coefficients ----
  if (adaptive) make_grid_windowed(r_lo_, rcut_, cfg_.nrho, rwin_lo, rcut_,
                                   crho_, srho_);
  else          make_grid(r_lo_, rcut_, cfg_.nrho, crho_, srho_);
  rho_stats_ = lsq_fit_windowed(xr, yrho, rwin_lo, rcut_, crho_, srho_,
                                true, arho_, false, 0.0, 0.0);

  // ---- 3. embedding F(rho) over [0, rhomax] -> many-body linear weights ----
  const int nrho_samp = std::max(2048, cfg_.nF * cfg_.oversample);
  std::vector<double> xrho, yF;
  xrho.reserve(nrho_samp); yF.reserve(nrho_samp);
  for (int i = 0; i < nrho_samp; ++i) {
    const double rho = rhomax_ * i / (nrho_samp - 1);
    xrho.push_back(rho);
    yF.push_back(eam_.F_tab(rho));
  }
  const double rho_win_hi = (adaptive && rho_reach_ > 0.0) ? rho_reach_
                                                           : rhomax_;
  if (adaptive) make_grid_windowed(0.0, rhomax_, cfg_.nF, 0.0, rho_win_hi,
                                   cF_, sF_);
  else          make_grid(0.0, rhomax_, cfg_.nF, cF_, sF_);
  // pin = the isolated-atom reference F(0) (the table's own value)
  F_stats_ = lsq_fit_windowed(xrho, yF, 0.0, rho_win_hi, cF_, sF_, false,
                              cmb_, cfg_.constrain_F0, 0.0, eam_.F_tab(0.0));

  // ---- concatenated linear weights [c2b | cmb] ----
  weights_.clear();
  weights_.reserve(c2b_.size() + cmb_.size());
  for (double v : c2b_) weights_.push_back(v);
  for (double v : cmb_) weights_.push_back(v);
}

namespace {
inline int deg_at(const std::vector<int> &d, std::size_t k) {
  return (k < d.size()) ? d[k] : 5;
}
}  // namespace

double EAMRefit::phi_fit(double r) const {
  const double fc = cutoff(r);
  if (fc == 0.0) return 0.0;
  double acc = 0.0;
  for (std::size_t k = 0; k < cgrid2b_.size(); ++k)
    acc += c2b_[k] * basisval(r, cgrid2b_[k], sgrid2b_[k], false,
                              deg_at(deg2b_, k));
  return acc * fc;
}

double EAMRefit::rho_fit(double r) const {
  const double fc = cutoff(r);
  if (fc == 0.0) return 0.0;
  double acc = 0.0;
  for (std::size_t n = 0; n < crho_.size(); ++n)
    acc += arho_[n] * basisval(r, crho_[n], srho_[n], false,
                               deg_at(degrho_, n));
  return acc * fc;
}

double EAMRefit::F_fit(double rho) const {
  double acc = 0.0;
  for (std::size_t m = 0; m < cF_.size(); ++m)
    acc += cmb_[m] * basisval(rho, cF_[m], sF_[m], /*emb=*/true,
                              deg_at(degF_, m));
  return acc;
}

// Below r_lo the most faithful tabulation is the ORIGINAL table itself, shifted
// by the (small, ~fit-error-sized) offset that makes it C0-continuous with the
// fit at the join. No invented functional form, no hyper-parameters; sentinel
// walls and deep-core shapes are preserved exactly as the source potential has
// them. (A C2-matched exp continuation was tried first and ABANDONED: it
// faithfully extrapolates the blip-edge curvature artefact and collapses to
// zero — see cycle2 report.)
void EAMRefit::set_blip_fit(Fn which, const std::vector<double> &mu,
                            const std::vector<double> &eta,
                            const std::vector<double> &c,
                            const std::vector<int> &deg) {
  auto assign = [](tadah::core::v_type &dst, const std::vector<double> &src) {
    dst.assign(src.begin(), src.end());
  };
  switch (which) {
    case Fn::PHI: assign(cgrid2b_, mu); assign(sgrid2b_, eta);
                  assign(c2b_, c); deg2b_ = deg; break;
    case Fn::RHO: assign(crho_, mu); assign(srho_, eta);
                  assign(arho_, c); degrho_ = deg; break;
    case Fn::F:   assign(cF_, mu); assign(sF_, eta);
                  assign(cmb_, c); degF_ = deg; break;
  }
}

void EAMRefit::finalize_injected(double rwin_lo, double rho_win_hi) {
  weights_.clear();
  weights_.reserve(c2b_.size() + cmb_.size());
  for (double v : c2b_) weights_.push_back(v);
  for (double v : cmb_) weights_.push_back(v);

  auto stats_on = [&](const std::function<double(double)> &fitf,
                      const std::function<double(double)> &tabf, double lo,
                      double hi) {
    FitStats st;
    const int N = 1000;
    double se = 0, ae = 0, mx = 0;
    for (int i = 0; i < N; ++i) {
      const double x = lo + (hi - lo) * i / (N - 1);
      const double d = fitf(x) - tabf(x);
      se += d * d; ae += std::fabs(d); mx = std::max(mx, std::fabs(d));
    }
    st.n = N; st.rmse = std::sqrt(se / N); st.mae = ae / N; st.maxabs = mx;
    return st;
  };
  phi_stats_ = stats_on([&](double r) { return phi_fit(r); },
                        [&](double r) { return eam_.phi_tab(r); }, rwin_lo,
                        rcut_);
  rho_stats_ = stats_on([&](double r) { return rho_fit(r); },
                        [&](double r) { return eam_.rho_tab(r); }, rwin_lo,
                        rcut_);
  F_stats_ = stats_on([&](double x) { return F_fit(x); },
                      [&](double x) { return eam_.F_tab(x); }, 0.0,
                      rho_win_hi);
}

double EAMRefit::F_setfl(double rho) const {
  if (rho <= rho_join_) return F_fit(rho);
  return eam_.F_tab(rho) + (F_fit(rho_join_) - eam_.F_tab(rho_join_));
}

double EAMRefit::z2_setfl(double r) const {
  if (r >= r_lo_) return phi_fit(r) * r;
  // original z2r + (C0 offset in phi) * r — finite at r = 0
  return eam_.z2_tab(r) + (phi_fit(r_lo_) - eam_.phi_tab(r_lo_)) * r;
}

double EAMRefit::phi_setfl(double r) const {
  if (r >= r_lo_) return phi_fit(r);
  return r > 1e-12 ? z2_setfl(r) / r : 0.0;
}

double EAMRefit::rho_setfl(double r) const {
  if (r >= r_lo_) return rho_fit(r);
  return eam_.rho_tab(r) + (rho_fit(r_lo_) - eam_.rho_tab(r_lo_));
}

tadah::core::Context EAMRefit::build_pot_context() const {
  tadah::core::Context c;
  std::string el = eam_.element();
  if (el.empty()) el = "X";
  const int bcode = static_cast<int>(family_);  // DM_REAM TYPEMB basis code

  c.add("MODEL", "M_BLR");
  c.add("MODEL", "BF_Linear");
  c.add("ATOM", el);
  // Pin the density weight explicitly: when WATOM is ABSENT the descriptors
  // default it to 1.0, but the HPO pipeline fills missing keys from blueprint
  // defaults (atomic number!) — which scaled rho by Z and produced ~500 eV/atom
  // garbage in the hpo evaluation of an otherwise-exact pot. Measured 2026-06-10.
  c.add("WATOM", 1.0);
  c.add("BIAS", "false");
  c.add("NORM", "false");
  c.add("FORCE", "false");
  c.add("STRESS", "false");
  c.add("LAMBDA", 0);
  c.add("OALGO", 2);

  // ---- two-body descriptor: must match the family the coefficients were
  //      FIT with (basisval), else the pot evaluates a different basis than
  //      was trained — D2_BP is the Gaussian twin of D2_Blip (same
  //      CGRID2B/SGRID2B keys); D2_Knot5 is the one-sided knot descriptor
  //      (knot positions only — NO width key, eta is a pure gauge; optional
  //      KDEG2B data list carries per-atom degrees 3|5).
  c.add("INIT2B", "true");
  const bool knot2b = (family_ == BasisKind::Knot5);
  c.add("TYPE2B", family_ == BasisKind::Gaussian ? "D2_BP"
                : knot2b                      ? "D2_Knot5"
                : family_ == BasisKind::Blip5    ? "D2_Blip5"
                                              : "D2_Blip");
  if (knot2b) {
    // contract: D2_Knot5 <nC> <nKdeg> El El
    c.add("TYPE2B", static_cast<int>(cgrid2b_.size()));
    c.add("TYPE2B", static_cast<int>(deg2b_.size()));
  } else {
    // contract: <name> <nC> <nS> El El
    c.add("TYPE2B", static_cast<int>(cgrid2b_.size()));
    c.add("TYPE2B", static_cast<int>(sgrid2b_.size()));
  }
  c.add("TYPE2B", el);
  c.add("TYPE2B", el);
  for (double v : cgrid2b_) c.add("CGRID2B", v);
  if (!knot2b)
    for (double v : sgrid2b_) c.add("SGRID2B", v);
  // optional per-atom degrees (KnotMix/Knot3; absent = all quintic)
  for (int d : deg2b_) c.add("KDEG2B", static_cast<double>(d));
  c.add("RCUT2B", rcut_);
  c.add("RCTYPE2B", cfg_.rctype);

  // ---- many-body: DM_REAM ----
  // contract: DM_REAM <bcode> <nC> <nS> <nAmp> <nEc> <nEs> <nKd> <nKe> El El
  // (knot family: width sizes are 0 — eta is a pure gauge, no SGRID data)
  const bool knotmb = (family_ == BasisKind::Knot5);
  c.add("INITMB", "true");
  c.add("TYPEMB", "DM_REAM");
  c.add("TYPEMB", bcode);
  c.add("TYPEMB", static_cast<int>(crho_.size()));
  c.add("TYPEMB", knotmb ? 0 : static_cast<int>(srho_.size()));
  c.add("TYPEMB", static_cast<int>(arho_.size()));
  c.add("TYPEMB", static_cast<int>(cF_.size()));
  c.add("TYPEMB", knotmb ? 0 : static_cast<int>(sF_.size()));
  c.add("TYPEMB", static_cast<int>(degrho_.size()));
  c.add("TYPEMB", static_cast<int>(degF_.size()));
  c.add("TYPEMB", el);
  c.add("TYPEMB", el);
  for (double v : crho_) c.add("CGRIDMB", v);
  if (!knotmb)
    for (double v : srho_) c.add("SGRIDMB", v);
  for (double v : arho_) c.add("AMPGRIDMB", v);
  for (int d : degrho_) c.add("KDEGMB", static_cast<double>(d));
  for (double v : cF_) c.add("CEMBFUNC", v);
  if (!knotmb)
    for (double v : sF_) c.add("SEMBFUNC", v);
  for (int d : degF_) c.add("KDEGEMB", static_cast<double>(d));
  c.add("RCUTMB", rcut_);
  c.add("RCTYPEMB", cfg_.rctype);

  // ---- linear weights [c2b | cmb] ----
  for (double v : weights_) c.add("WEIGHTS", v);

  // ---- HPO coefficient-search seed ----
  // FIXINDEX + FIXWEIGHT freeze every linear weight at its fitted value so HPO
  // can *search* the coefficients via `OPTIM FIXWEIGHT` (the standard linear
  // solve would otherwise overwrite them each iteration). `tadah predict`
  // ignores both keys — it reads WEIGHTS — so this does not change predictions.
  // `OPTIM FIXWEIGHT (i)` indices name MODEL WEIGHTS; because FIXINDEX is 1-K
  // here, weight index == FIXWEIGHT slot and the emitted target is valid
  // under both readings.
  const std::size_t K = weights_.size();
  if (K > 0) {
    c.add("FIXINDEX", std::string("1-") + std::to_string(K));
    for (double v : weights_) c.add("FIXWEIGHT", v);
  }

  return c;
}

} // namespace models
} // namespace tadah
