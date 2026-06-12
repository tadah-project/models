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

DM_REAM::DM_REAM(const tadah::core::Context& ctx) : DM_Base(ctx) {
  init();
  if (!ctx.get<bool>("INITMB")) return;

  int bcode = 0;
  try { bcode = ctx.get<int>("TYPEMB", 1); } catch (...) { bcode = 0; }
  basis = basis_from_bcode(bcode);

  get_grid(ctx, "CGRIDRHOMB", crho);
  get_grid(ctx, "SGRIDRHOMB", srho);
  get_grid(ctx, "AGRIDMB",    arho);
  get_grid(ctx, "CGRIDFMB",   cemb);
  get_grid(ctx, "SGRIDFMB",   semb);
  if (ctx.contains("DGRIDRHOMB") && ctx.size("DGRIDRHOMB") > 0)
    get_grid(ctx, "DGRIDRHOMB", drho_deg);  // optional knot degrees (3|5)
  if (ctx.contains("DGRIDFMB") && ctx.size("DGRIDFMB") > 0)
    get_grid(ctx, "DGRIDFMB", demb_deg);
  D_Base::printGrid(crho, "CGRIDRHOMB");
  D_Base::printGrid(cemb, "CGRIDFMB");

  if (crho.size() != srho.size() || crho.size() != arho.size())
    throw std::runtime_error(
        "DM_REAM: CGRIDRHOMB / SGRIDRHOMB / AGRIDMB sizes differ.");
  if (cemb.size() != semb.size())
    throw std::runtime_error("DM_REAM: CGRIDFMB / SGRIDFMB sizes differ.");
  if (cemb.empty())
    throw std::runtime_error("DM_REAM: embedding basis (CGRIDFMB) is empty.");
  if (!drho_deg.empty() && drho_deg.size() != crho.size())
    throw std::runtime_error("DM_REAM: DGRIDRHOMB / CGRIDRHOMB sizes differ.");
  if (!demb_deg.empty() && demb_deg.size() != cemb.size())
    throw std::runtime_error("DM_REAM: DGRIDFMB / CGRIDFMB sizes differ.");
  for (const auto &e : srho)
    if (std::abs(e) < std::numeric_limits<double>::min())
      throw std::runtime_error("DM_REAM: zero entry in SGRIDRHOMB.");
  for (const auto &e : semb)
    if (std::abs(e) < std::numeric_limits<double>::min())
      throw std::runtime_error("DM_REAM: zero entry in SGRIDFMB.");

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
  if (basis == Basis::Knot5)   // LEFT-active: vanishes smoothly at the knot
    return rho_deg(n) == 3 ? K3L(srho[n] * (r - crho[n]), fc)
                           : K5L(srho[n] * (r - crho[n]), fc);
  if (basis == Basis::Blip5)
    return B5(srho[n] * (r - crho[n]), fc);
  return G(r, srho[n], crho[n], fc);
}

double DM_REAM::dpsi(std::size_t n, double r, double fc, double fcp) const {
  if (basis == Basis::Blip)
    return dB(srho[n] * (r - crho[n]), srho[n], fc, fcp);
  if (basis == Basis::Knot5)
    return rho_deg(n) == 3 ? dK3L(srho[n] * (r - crho[n]), srho[n], fc, fcp)
                           : dK5L(srho[n] * (r - crho[n]), srho[n], fc, fcp);
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
  if (basis == Basis::Knot5)   // RIGHT-active: F free at high density
    return emb_deg(m) == 3 ? K3R(semb[m] * (rho - cemb[m]))
                           : K5R(semb[m] * (rho - cemb[m]));
  if (basis == Basis::Blip5)
    return B5(semb[m] * (rho - cemb[m]));
  return G(rho, semb[m], cemb[m]);
}

double DM_REAM::dgemb(std::size_t m, double rho) const {
  if (basis == Basis::Blip)
    return dB(semb[m] * (rho - cemb[m]), semb[m]);
  if (basis == Basis::Knot5)
    return emb_deg(m) == 3 ? dK3R(semb[m] * (rho - cemb[m]), semb[m])
                           : dK5R(semb[m] * (rho - cemb[m]), semb[m]);
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
  keys = {"CGRIDRHOMB", "SGRIDRHOMB", "AGRIDMB", "CGRIDFMB", "SGRIDFMB"};
  nparams = 1;   // TYPEMB: DM_REAM <basis> <El> <El>
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
  add_vec("CGRIDRHOMB", Scale::Distance, "Å",
          "Centre of the i-th density radial basis function.", 0.0,
          Criticality::Recommended);
  add_vec("SGRIDRHOMB", Scale::InverseDistanceSq, "Å^-2",
          "Inverse width of the i-th density radial basis function.", 1.0e-6,
          Criticality::Recommended);
  add_vec("AGRIDMB", Scale::Dimensionless, "",
          "Baked density-expansion coefficient a_n (fitted, not a weight).",
          -1.0e30, Criticality::Optional);
  add_vec("CGRIDFMB", Scale::Density, "rho",
          "Centre (in rho-space) of the i-th embedding basis function.",
          -1.0e30, Criticality::Recommended);
  add_vec("SGRIDFMB", Scale::InverseDensity, "rho^-1",
          "Inverse width (in rho-space) of the i-th embedding basis function.",
          1.0e-6, Criticality::Recommended);

  dm.expected_columns = [](const SeedView& s) -> std::size_t {
    if (!s.contains("CGRIDFMB")) return 0;
    return s.size("CGRIDFMB");
  };

  return dm;
}

bool DM_REAM::reinit(const tadah::core::Context& ctx) {
  const std::size_t old_s = s;
  int bcode = 0;
  try { bcode = ctx.get<int>("TYPEMB", 1); } catch (...) { bcode = 0; }
  basis = basis_from_bcode(bcode);
  get_grid(ctx, "CGRIDRHOMB", crho);
  get_grid(ctx, "SGRIDRHOMB", srho);
  get_grid(ctx, "AGRIDMB",    arho);
  get_grid(ctx, "CGRIDFMB",   cemb);
  get_grid(ctx, "SGRIDFMB",   semb);
  drho_deg.clear();
  demb_deg.clear();
  if (ctx.contains("DGRIDRHOMB") && ctx.size("DGRIDRHOMB") > 0)
    get_grid(ctx, "DGRIDRHOMB", drho_deg);
  if (ctx.contains("DGRIDFMB") && ctx.size("DGRIDFMB") > 0)
    get_grid(ctx, "DGRIDFMB", demb_deg);
  s         = cemb.size();
  rhoisize  = 1;
  rhoipsize = s;
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUTMBMAX"));
  return s != old_s;
}

} // namespace models
} // namespace tadah
