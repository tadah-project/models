#include <tadah/models/descriptors/dm/dm_ream.h>
#include <tadah/models/descriptors/basis_dispatch.h>
#include <tadah/core/logger.h>
#include <tadah/core/stats/data_stats.h>
#include <algorithm>
#include <atomic>

#include <stdexcept>

namespace tadah {
namespace models {

namespace {
// One notice per process when a DM_REAM model is built from a context that
// does not declare EMBF0 (the default ZERO applies). A trained potential
// (WEIGHTS stored) gets a WARNING: if it was written by tadah 1.5.0 or
// earlier its energies change.
std::atomic<bool> embf0_noticed{false};
void notice_missing_embf0(const tadah::core::Context &ctx) {
  if (embf0_noticed.exchange(true)) return;
  auto &log = tadah::core::Logger::instance();
  if (!ctx.has_stored("WEIGHTS")) {
    log.log(tadah::core::LogLevel::INFO,
            "DM_REAM: EMBF0 not set; using the default EMBF0 ZERO "
            "(zero-referenced embedding, F(0) = 0 by construction).");
    return;
  }
  if (model_is_linear(ctx))
    log.log(tadah::core::LogLevel::WARNING,
            "DM_REAM: EMBF0 not set; this potential is evaluated with the "
            "default EMBF0 ZERO (F(0) = 0 by construction). A potential "
            "written by tadah 1.5.0 or earlier now gives every atomic energy "
            "shifted by -F(0) (forces and stresses are unchanged for this "
            "linear model); add 'EMBF0 FREE' to the potential to reproduce "
            "its energies.");
  else
    log.log(tadah::core::LogLevel::WARNING,
            "DM_REAM: EMBF0 not set; this potential is evaluated with the "
            "default EMBF0 ZERO (F(0) = 0 by construction). The model is "
            "non-linear in the descriptor, so for a potential written by "
            "tadah 1.5.0 or earlier its energies AND forces change; add "
            "'EMBF0 FREE' to the potential to reproduce them.");
}
}  // namespace

DM_REAM::DM_REAM() { init(); }

void DM_REAM::validate_grids() const {
  if (crho.size() != arho.size())
    throw std::runtime_error("DM_REAM: CGRIDMB / AMPGRIDMB sizes differ.");
  if (cemb.empty())
    throw std::runtime_error("DM_REAM: embedding basis (CEMBFUNC) is empty.");
  if (basis_traits(basis).eta_is_gauge) {
    // eta is a pure gauge for this family (knots): the width channels must be
    // ABSENT (size-0 segments in TYPEMB) — no silent constant-1 lists. Asking
    // the trait (not comparing to a specific basis) keeps this future-proof.
    if (!srho.empty() || !semb.empty())
      throw std::runtime_error(
          "DM_REAM (knot family): SGRIDMB/SEMBFUNC must be empty — eta is a "
          "pure gauge ((c,eta)->(c s^k, eta/s)); the knot position grids "
          "CGRIDMB/CEMBFUNC define the basis alone.");
  } else {
    if (crho.size() != srho.size())
      throw std::runtime_error("DM_REAM: CGRIDMB / SGRIDMB sizes differ.");
    if (cemb.size() != semb.size())
      throw std::runtime_error("DM_REAM: CEMBFUNC / SEMBFUNC sizes differ.");
    for (const auto &e : srho)
      if (std::abs(e) < std::numeric_limits<double>::min())
        throw std::runtime_error("DM_REAM: zero entry in SGRIDMB.");
    for (const auto &e : semb)
      if (std::abs(e) < std::numeric_limits<double>::min())
        throw std::runtime_error("DM_REAM: zero entry in SEMBFUNC.");
  }
  if (!drho_deg.empty() && drho_deg.size() != crho.size())
    throw std::runtime_error("DM_REAM: KDEGMB / CGRIDMB sizes differ.");
  if (!demb_deg.empty() && demb_deg.size() != cemb.size())
    throw std::runtime_error("DM_REAM: KDEGEMB / CEMBFUNC sizes differ.");
  for (const auto *arr : {&drho_deg, &demb_deg})
    for (const auto &d : *arr)
      if (d != 3.0 && d != 5.0)
        throw std::runtime_error(
            "DM_REAM: KDEG* entries must be 3 or 5 (the knot-degree lists "
            "are structure data, not optimisable parameters).");
}

DM_REAM::DM_REAM(const tadah::core::Context& ctx) : DM_Base(ctx) {
  init();
  if (!ctx.get<bool>("INITMB")) return;

  int bcode = 0;
  try { bcode = ctx.get<int>("TYPEMB", 1); } catch (...) { bcode = 0; }
  basis = basis_kind_from_code(bcode);

  get_grid(ctx, "CGRIDMB", crho);
  srho.clear();
  if (ctx.contains("SGRIDMB") && ctx.size("SGRIDMB") > 0)
    get_grid(ctx, "SGRIDMB", srho);     // absent for the knot family
  get_grid(ctx, "AMPGRIDMB", arho);
  get_grid(ctx, "CEMBFUNC", cemb);
  semb.clear();
  if (ctx.contains("SEMBFUNC") && ctx.size("SEMBFUNC") > 0)
    get_grid(ctx, "SEMBFUNC", semb);    // absent for the knot family
  if (ctx.contains("KDEGMB") && ctx.size("KDEGMB") > 0)
    get_grid(ctx, "KDEGMB", drho_deg);  // optional knot degrees (3|5)
  if (ctx.contains("KDEGEMB") && ctx.size("KDEGEMB") > 0)
    get_grid(ctx, "KDEGEMB", demb_deg);
  D_Base::printGrid(crho, "CGRIDMB");
  D_Base::printGrid(cemb, "CEMBFUNC");

  validate_grids();
  setup_embedding_reference(ctx);

  s         = cemb.size();   // M embedding outputs (linear weights live in M_BLR)
  rhoisize  = 1;             // scalar density rho_i
  rhoipsize = s;             // one stored g'_m(rho_i) per output (for forces)
  set_rfidx(0);
  set_drfidx(rhoisize);

  init_for_atoms(get_init_atoms(ctx));
  warn_if_knots_beyond_rcut();
  warn_if_bump_support_beyond_rcut("DM_REAM", "CGRIDMB", "RCUTMB", crho,
                                   srho, basis);
}

void DM_REAM::warn_if_knots_beyond_rcut() {
  if (warned_rcut_ || basis != BasisKind::Knot5 || crho.empty() || !fcut)
    return;
  const double mx = *std::max_element(crho.begin(), crho.end());
  const double rc = get_rcut();
  if (mx > rc + 1e-12) {
    warned_rcut_ = true;
    tadah::core::Logger::instance().log(
        tadah::core::LogLevel::WARNING,
        "DM_REAM: the largest CGRIDMB density knot (" + std::to_string(mx) +
        ") lies beyond RCUTMB (" + std::to_string(rc) +
        ") — its contribution is silently hard-truncated at the cutoff. "
        "The effective cutoff of a knot channel is its largest knot: set "
        "RCUTMB to the largest knot (refit/edit/hpo tooling derives it "
        "automatically).");
  }
}

// ---------------------------- basis helpers ----------------------------

double DM_REAM::psi(std::size_t n, double r, double fc) const {
  // Knot family: SGRIDMB is absent (eta is a pure gauge) -> s = 1 (ignored by
  // the knot arms of the dispatcher). One shared switch, no if(basis==...) ladder.
  const double s = (n < srho.size()) ? srho[n] : 1.0;
  return radial_value(basis, r, fc, s, crho[n], rho_deg(n));
}

double DM_REAM::dpsi(std::size_t n, double r, double fc, double fcp) const {
  const double s = (n < srho.size()) ? srho[n] : 1.0;
  return radial_deriv(basis, r, fc, fcp, s, crho[n], rho_deg(n));
}

double DM_REAM::rho_raw(double r, double fc) const {
  double acc = 0.0;
  for (std::size_t n = 0; n < crho.size(); ++n) acc += arho[n] * psi(n, r, fc);
  return acc;
}

double DM_REAM::rho_raw_deriv(double r, double fc, double fcp) const {
  double acc = 0.0;
  for (std::size_t n = 0; n < crho.size(); ++n) acc += arho[n] * dpsi(n, r, fc, fcp);
  return acc;
}

double DM_REAM::gemb(std::size_t m, double rho) const {
  // Embedding: RIGHT-active knots, no cutoff. SEMBFUNC absent for the knot
  // family (eta gauge) -> s = 1 (ignored by the knot arms).
  const double s = (m < semb.size()) ? semb[m] : 1.0;
  return embedding_value(basis, rho, s, cemb[m], emb_deg(m));
}

double DM_REAM::dgemb(std::size_t m, double rho) const {
  const double s = (m < semb.size()) ? semb[m] : 1.0;
  return embedding_deriv(basis, rho, s, cemb[m], emb_deg(m));
}

void DM_REAM::setup_embedding_reference(const tadah::core::Context& ctx) {
  const bool declared = ctx.has_stored("EMBF0");
  // absent -> default ZERO, materialised; throws on an invalid value before
  // any member changes
  const EmbF0Mode mode = embf0_from_context(ctx);
  if (!declared) notice_missing_embf0(ctx);
  tadah::core::v_type g0a, g0e;
  if (mode == EmbF0Mode::Zero) {   // FREE: raw basis, no reference, no FP work
    const std::size_t M = cemb.size();
    g0a.resize(M);
    g0e.resize(M);
    for (std::size_t m = 0; m < M; ++m) {
      g0a[m] = gemb(m, 0.0);                       // calc_aed's routine
      const double s = (m < semb.size()) ? semb[m] : 1.0;
      double v, d;
      embedding_value_deriv(basis, 0.0, s, cemb[m], v, d, emb_deg(m));
      g0e[m] = v;                                  // eam_embed's routine
    }
  }
  embf0_ = mode;
  g0_aed_.swap(g0a);
  g0_emb_.swap(g0e);
}

// ------------------------------ density --------------------------------

void DM_REAM::calc_rho(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const tadah::core::Vec3d &,
  tadah::core::rho_type& rho,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  const double fc = fcut->calc(rij);
  rho(rfidx) += scale * weights[Zi][Zj] * rho_raw(rij, fc);
}

// --------------------------- embedding energy --------------------------

void DM_REAM::calc_aed(
  tadah::core::rho_type& rho,
  tadah::core::aed_type &aed) const
{
  const double rho_i = rho(rfidx);
  const std::size_t M = cemb.size();
  if (embf0_ == EmbF0Mode::Free) {
    for (std::size_t m = 0; m < M; ++m) {
      aed(fidx + m) += gemb(m, rho_i);   // F(rho) = sum_m c_m g_m(rho)  (c_m in M_BLR)
      rho(drfidx + m) = dgemb(m, rho_i); // store g'_m(rho_i) for the force descriptor
    }
    return;
  }
  // EMBF0 ZERO: F(rho) = sum_m c_m (g_m(rho) - g_m(0)). At rho == 0 every
  // output is exactly zero, so nothing is added (exact whatever the FP
  // contraction); g'_m(0) is still stored for the force descriptor.
  if (rho_i == 0.0) {
    for (std::size_t m = 0; m < M; ++m) rho(drfidx + m) = dgemb(m, rho_i);
    return;
  }
  for (std::size_t m = 0; m < M; ++m) {
    aed(fidx + m) += gemb(m, rho_i) - g0_aed_[m];
    rho(drfidx + m) = dgemb(m, rho_i);
  }
}

// ------------------------------- forces --------------------------------

void DM_REAM::calc_dXijdri_dXjidri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rhoi,
  tadah::core::rho_type& rhoj,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  const double fc  = fcut->calc(rij);
  const double fcp = fcut->calc_prime(rij);
  const double rpr = rho_raw_deriv(rij, fc, fcp);      // unweighted drho/dr
  const double rinv = 1.0 / rij;
  const double pref = scale * weights[Zi][Zj] * rpr * rinv;
  const std::size_t M = cemb.size();
  for (std::size_t m = 0; m < M; ++m) {
    const double coef = pref * (rhoi(drfidx + m) + rhoj(drfidx + m));
    fd_ij(fidx + m, 0) = coef * vec_ij[0];
    fd_ij(fidx + m, 1) = coef * vec_ij[1];
    fd_ij(fidx + m, 2) = coef * vec_ij[2];
  }
}

void DM_REAM::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rhoi,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  const double fc  = fcut->calc(rij);
  const double fcp = fcut->calc_prime(rij);
  const double rpr = rho_raw_deriv(rij, fc, fcp);      // unweighted drho/dr
  const double rinv = 1.0 / rij;
  const double pref = scale * weights[Zi][Zj] * rpr * rinv;
  const std::size_t M = cemb.size();
  for (std::size_t m = 0; m < M; ++m) {
    const double coef = pref * rhoi(drfidx + m);       // g'_m(rho_i) from calc_aed
    fd_ij(fidx + m, 0) = coef * vec_ij[0];
    fd_ij(fidx + m, 1) = coef * vec_ij[1];
    fd_ij(fidx + m, 2) = coef * vec_ij[2];
  }
}

// --------------------------- EAM collapse ------------------------------

void DM_REAM::eam_density(const int Zi, const int Zj, const double rij,
                          double &rho, double &rhop) const {
  rho = 0.0; rhop = 0.0;
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  const double fc  = fcut->calc(rij);
  const double fcp = fcut->calc_prime(rij);
  double rho_raw_both = 0.0;
  double rhop_raw = 0.0;
  for (std::size_t n = 0; n < crho.size(); ++n) {
    const double s = (n < srho.size()) ? srho[n] : 1.0;
    double psi_value, psi_deriv;
    radial_value_deriv(basis, rij, fc, fcp, s, crho[n],
                       psi_value, psi_deriv, rho_deg(n));
    rho_raw_both += arho[n] * psi_value;
    rhop_raw += arho[n] * psi_deriv;
  }
  const double w   = weights[Zi][Zj];
  rho  = w * rho_raw_both;
  rhop = w * rhop_raw;
}

void DM_REAM::eam_embed(const tadah::core::t_type &w_full, const double rho_bar,
                        double &F, double &Fp) const {
  F = 0.0; Fp = 0.0;
  const std::size_t M = cemb.size();
  if (embf0_ == EmbF0Mode::Free) {
    for (std::size_t m = 0; m < M; ++m) {
      const double w = w_full(fidx + m);
      const double s = (m < semb.size()) ? semb[m] : 1.0;
      double emb_value, emb_deriv;
      embedding_value_deriv(basis, rho_bar, s, cemb[m],
                            emb_value, emb_deriv, emb_deg(m));
      F  += w * emb_value;
      Fp += w * emb_deriv;
    }
    return;
  }
  // EMBF0 ZERO: the same zero-referenced basis as calc_aed, subtracted per
  // term so this kernel mirrors w.aed; F' is untouched.
  for (std::size_t m = 0; m < M; ++m) {
    const double w = w_full(fidx + m);
    const double s = (m < semb.size()) ? semb[m] : 1.0;
    double emb_value, emb_deriv;
    embedding_value_deriv(basis, rho_bar, s, cemb[m],
                          emb_value, emb_deriv, emb_deg(m));
    F  += w * (emb_value - g0_emb_[m]);
    Fp += w * emb_deriv;
  }
  if (rho_bar == 0.0) F = 0.0;
}

std::optional<double>
DM_REAM::raw_embedding_f0(const tadah::core::t_type &w_full) const {
  const std::size_t M = cemb.size();
  if (w_full.size() < fidx + M) return std::nullopt;
  double acc = 0.0;
  for (std::size_t m = 0; m < M; ++m) acc += w_full(fidx + m) * gemb(m, 0.0);
  return acc;
}

// ------------------------------- misc ----------------------------------

std::string DM_REAM::label() const { return lab; }

void DM_REAM::init_rhoi(tadah::core::rho_type& rhoi) {
  rhoi.resize(rhoisize + rhoipsize);
  rhoi.set_zero();
}

void DM_REAM::init() {
  // contract (see d_base.h): TYPEMB carries 1 family param (the basis code)
  // followed by one SIZE token per keys entry, then per data_keys entry:
  //   DM_REAM <bcode> <nC> <nS> <nAmp> <nEc> <nEs> <nKd> <nKe> El El
  // Knot pots carry 0 for the width channels (eta is a pure gauge) and bump
  // pots carry 0 for the KDEG channels.
  keys = {"CGRIDMB", "SGRIDMB", "AMPGRIDMB", "CEMBFUNC", "SEMBFUNC"};
  data_keys = {"KDEGMB", "KDEGEMB"};
  nparams = 1 + 5 + 2;
}

std::optional<DescriptorMeta> DM_REAM::meta() const {
  using DataStats = tadah::core::stats::DataStats;
  static constexpr std::size_t I_P50 = tadah::core::stats::IDX_P50;
  static constexpr std::size_t I_P99 = tadah::core::stats::IDX_P99;

  DescriptorMeta dm;
  dm.label = "DM_REAM";

  // RCUTMB — real-space cutoff for the density neighbour sum.
  ParamMeta p_rcut;
  p_rcut.key     = "RCUTMB";
  p_rcut.arg_pos = -1;
  p_rcut.scale   = Scale::Distance;
  p_rcut.units   = "Å";
  p_rcut.role    = "Real-space cutoff for the EAM density neighbour sum.";
  p_rcut.hard.lo = 0.0;
  p_rcut.crit    = Criticality::Critical;
  p_rcut.scope   = PairScope::Global;
  p_rcut.suggest = [](const DataStats& d, const SeedView&, int)
      -> std::pair<double, double> {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    return std::make_pair(pp.r_percentiles[I_P50], pp.r_percentiles[I_P99]);
  };
  dm.params.push_back(std::move(p_rcut));

  // Density-expansion grids (r-space). Tunable, but `tadah refit` emits
  // explicit ±% OPTIM bounds, so the data-driven suggesters stay neutral.
  auto add_vec = [&dm](const std::string& key, Scale sc, const std::string& units,
                       const std::string& role, double lo, Criticality crit) {
    ParamMeta p;
    p.key     = key;
    p.arg_pos = 0;
    p.scale   = sc;
    p.units   = units;
    p.role    = role;
    p.hard.lo = lo;
    p.crit    = crit;
    p.scope   = PairScope::Global;
    p.suggest = [](const DataStats&, const SeedView&, int)
        -> std::pair<double, double> { return std::make_pair(0.0, 0.0); };
    dm.params.push_back(std::move(p));
  };
  add_vec("CGRIDMB", Scale::Distance, "Å",
          "Centre (knot families: knot position) of the i-th density radial "
          "basis function. Optimizable (OPTIM CGRIDMB).", 0.0,
          Criticality::Recommended);
  add_vec("SGRIDMB", Scale::InverseDistanceSq, "Å^-2",
          "Inverse width of the i-th density radial basis function — bump "
          "families only; ABSENT for knot families (eta is a pure gauge). "
          "Optimizable where present.", 1.0e-6, Criticality::Recommended);
  add_vec("AMPGRIDMB", Scale::Dimensionless, "",
          "Density-expansion amplitude a_n: rho(r) = sum_n a_n psi_n(r). "
          "Excluded from the LINEAR solve by bilinearity, but a first-class "
          "NONLINEAR HPO parameter (OPTIM AMPGRIDMB).", -1.0e30,
          Criticality::Recommended);
  add_vec("CEMBFUNC", Scale::Density, "rho",
          "Centre (in rho-space) of the i-th embedding basis function. "
          "Optimizable (OPTIM CEMBFUNC).", -1.0e30,
          Criticality::Recommended);
  add_vec("SEMBFUNC", Scale::InverseDensity, "rho^-1",
          "Inverse width (in rho-space) of the i-th embedding basis "
          "function — bump families only; ABSENT for knot families. "
          "Optimizable where present.", 1.0e-6, Criticality::Recommended);
  add_vec("KDEGMB", Scale::Dimensionless, "",
          "Knot-degree DATA (3|5) per density atom: a discrete structure "
          "choice, NOT optimizable (data_keys; values outside {3,5} are "
          "rejected at construction — the de-facto OPTIM guard).", 3.0,
          Criticality::Optional);
  add_vec("KDEGEMB", Scale::Dimensionless, "",
          "Knot-degree DATA (3|5) per embedding atom: NOT optimizable "
          "(data_keys).", 3.0, Criticality::Optional);

  dm.expected_columns = [](const SeedView& s) -> std::size_t {
    if (!s.contains("CEMBFUNC")) return 0;
    return s.size("CEMBFUNC");
  };

  return dm;
}

bool DM_REAM::reinit(const tadah::core::Context& ctx) {
  const std::size_t old_s = s;
  int bcode = 0;
  try { bcode = ctx.get<int>("TYPEMB", 1); } catch (...) { bcode = 0; }
  basis = basis_kind_from_code(bcode);
  get_grid(ctx, "CGRIDMB", crho);
  srho.clear();
  if (ctx.contains("SGRIDMB") && ctx.size("SGRIDMB") > 0)
    get_grid(ctx, "SGRIDMB", srho);     // absent for the knot family
  get_grid(ctx, "AMPGRIDMB", arho);
  get_grid(ctx, "CEMBFUNC", cemb);
  semb.clear();
  if (ctx.contains("SEMBFUNC") && ctx.size("SEMBFUNC") > 0)
    get_grid(ctx, "SEMBFUNC", semb);    // absent for the knot family
  drho_deg.clear();
  demb_deg.clear();
  if (ctx.contains("KDEGMB") && ctx.size("KDEGMB") > 0)
    get_grid(ctx, "KDEGMB", drho_deg);
  if (ctx.contains("KDEGEMB") && ctx.size("KDEGEMB") > 0)
    get_grid(ctx, "KDEGEMB", demb_deg);
  validate_grids();
  setup_embedding_reference(ctx);   // unconditional: CEMBFUNC/SEMBFUNC may move
  s         = cemb.size();
  rhoisize  = 1;
  rhoipsize = s;
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUTMBMAX"));
  warn_if_knots_beyond_rcut();
  warn_if_bump_support_beyond_rcut("DM_REAM", "CGRIDMB", "RCUTMB", crho,
                                   srho, basis);
  return s != old_s;
}

} // namespace models
} // namespace tadah
