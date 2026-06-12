#include <tadah/models/descriptors/dm/dm_ream.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <tadah/core/stats/data_stats.h>

#include <stdexcept>

namespace tadah {
namespace models {

DM_REAM::DM_REAM() { init(); }

namespace {
// Map the TYPEMB basis index to the family enum. Unknown codes ERROR rather
// than silently degrading to Blip (a pot written by a newer Tadah would
// otherwise be evaluated with the wrong basis without any warning).
DM_REAM::Basis basis_from_bcode(int bcode) {
  switch (bcode) {
    case 0: return DM_REAM::Basis::Blip;
    case 1: return DM_REAM::Basis::Gaussian;
    case 2: return DM_REAM::Basis::Knot5;
    case 3: return DM_REAM::Basis::Blip5;
    default:
      throw std::runtime_error(
          "DM_REAM: unknown basis code " + std::to_string(bcode) +
          " in TYPEMB (known: 0=Blip, 1=Gaussian, 2=Knot5, 3=Blip5).");
  }
}
}  // namespace

void DM_REAM::validate_grids() const {
  if (crho.size() != arho.size())
    throw std::runtime_error("DM_REAM: CGRIDMB / AMPGRIDMB sizes differ.");
  if (cemb.empty())
    throw std::runtime_error("DM_REAM: embedding basis (CEMBFUNC) is empty.");
  if (basis == Basis::Knot5) {
    // eta is a pure gauge for the knot family: the width channels must be
    // ABSENT (size-0 segments in TYPEMB) — no silent constant-1 lists.
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
  basis = basis_from_bcode(bcode);

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

  s         = cemb.size();   // M embedding outputs (linear weights live in M_BLR)
  rhoisize  = 1;             // scalar density rho_i
  rhoipsize = s;             // one stored g'_m(rho_i) per output (for forces)
  set_rfidx(0);
  set_drfidx(rhoisize);

  init_for_atoms(get_init_atoms(ctx));
}

// ---------------------------- basis helpers ----------------------------

double DM_REAM::psi(std::size_t n, double r, double fc) const {
  if (basis == Basis::Blip)
    return B(srho[n] * (r - crho[n]), fc);
  if (basis == Basis::Knot5)   // LEFT-active; eta == 1 (gauge, no SGRIDMB)
    return rho_deg(n) == 3 ? K3L(r - crho[n], fc) : K5L(r - crho[n], fc);
  if (basis == Basis::Blip5)
    return B5(srho[n] * (r - crho[n]), fc);
  return G(r, srho[n], crho[n], fc);
}

double DM_REAM::dpsi(std::size_t n, double r, double fc, double fcp) const {
  if (basis == Basis::Blip)
    return dB(srho[n] * (r - crho[n]), srho[n], fc, fcp);
  if (basis == Basis::Knot5)
    return rho_deg(n) == 3 ? dK3L(r - crho[n], 1.0, fc, fcp)
                           : dK5L(r - crho[n], 1.0, fc, fcp);
  if (basis == Basis::Blip5)
    return dB5(srho[n] * (r - crho[n]), srho[n], fc, fcp);
  return dG(r, srho[n], crho[n], fc, fcp);
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
  if (basis == Basis::Blip)
    return B(semb[m] * (rho - cemb[m]));
  if (basis == Basis::Knot5)   // RIGHT-active; eta == 1 (gauge, no SEMBFUNC)
    return emb_deg(m) == 3 ? K3R(rho - cemb[m]) : K5R(rho - cemb[m]);
  if (basis == Basis::Blip5)
    return B5(semb[m] * (rho - cemb[m]));
  return G(rho, semb[m], cemb[m]);
}

double DM_REAM::dgemb(std::size_t m, double rho) const {
  if (basis == Basis::Blip)
    return dB(semb[m] * (rho - cemb[m]), semb[m]);
  if (basis == Basis::Knot5)
    return emb_deg(m) == 3 ? dK3R(rho - cemb[m], 1.0)
                           : dK5R(rho - cemb[m], 1.0);
  if (basis == Basis::Blip5)
    return dB5(semb[m] * (rho - cemb[m]), semb[m]);
  return dG(rho, semb[m], cemb[m]);
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
  for (std::size_t m = 0; m < M; ++m) {
    aed(fidx + m) += gemb(m, rho_i);   // F(rho) = sum_m c_m g_m(rho)  (c_m in M_BLR)
    rho(drfidx + m) = dgemb(m, rho_i); // store g'_m(rho_i) for the force descriptor
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
  basis = basis_from_bcode(bcode);
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
  s         = cemb.size();
  rhoisize  = 1;
  rhoipsize = s;
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUTMBMAX"));
  return s != old_s;
}

} // namespace models
} // namespace tadah
