#include "tadah/models/descriptors/d_base.h"
#include <tadah/models/descriptors/d2/d2_knot5.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <tadah/models/descriptors/param_meta.h>
#include <tadah/core/stats/data_stats.h>
#include <cstdlib>
#include <stdexcept>
namespace tadah {
namespace models {


D2_Knot5::D2_Knot5() {
  init();
}
D2_Knot5::D2_Knot5(const tadah::core::Context &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  init();

  get_grid(c,"CGRID2B",cgrid);
  get_grid(c,"SGRID2B",sgrid);
  if (c.contains("DGRID2B") && c.size("DGRID2B") > 0)
    get_grid(c,"DGRID2B",dgrid);   // optional per-atom degree (3|5)
  D_Base::printGrid(cgrid,"CGRID2B");
  D_Base::printGrid(sgrid,"SGRID2B");

  if (cgrid.size()!=sgrid.size()) {
    throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
  }
  if (!dgrid.empty()) {
    if (dgrid.size() != cgrid.size())
      throw std::runtime_error("DGRID2B and CGRID2B arrays differ in size.\n");
    for (const auto &d : dgrid)
      if (d != 3.0 && d != 5.0)
        throw std::runtime_error("DGRID2B entries must be 3 or 5.\n");
  }

  for (const auto &e:sgrid) {
    if(std::abs(e)<std::numeric_limits<double>::min())
      throw std::runtime_error("At least one of SGRID2B values is zero.\n");
  }

  s=cgrid.size();

  auto init_atoms = get_init_atoms(c);
  init_for_atoms(init_atoms);
}

void D2_Knot5::calc_aed(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  tadah::core::aed_type &aed,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  for (size_t c=0; c<cgrid.size(); c++) {
    double t = kval(c, sgrid[c]*(rij-cgrid[c]), fc_ij);
    aed(i++) += scale*weights[Zi][Zj]*t;
  }
}
void D2_Knot5::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double ,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  for (size_t c=0; c<cgrid.size(); c++) {
    fd_ij(i++,0) = scale*weights[Zi][Zj]*kder(c, sgrid[c]*(rij-cgrid[c]),sgrid[c],fc_ij,fcp_ij);
  }
}
void D2_Knot5::calc_all(
  const int Zi,
  const int Zj,
  const double rij,
  const double ,
  tadah::core::aed_type &aed,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  for (size_t c=0; c<cgrid.size(); c++) {
    aed(i) += scale*weights[Zi][Zj]*kval(c, sgrid[c]*(rij-cgrid[c]),fc_ij);
    fd_ij(i,0) = scale*weights[Zi][Zj]*kder(c, sgrid[c]*(rij-cgrid[c]),sgrid[c],fc_ij,fcp_ij);
    ++i;
  }
}

// Closed-form gradient of D2_Knot5 w.r.t. its hyperparameters — the same
// shift-scale chain rule as D2_Blip (mirrors calc_aed's loop shape so a
// parity test against central-difference of calc_aed works component-wise).
// dK5L(x_c, 1.0) is the bare dK5L/dx_c.
bool D2_Knot5::calc_dDdtheta(
  const std::string& theta_name,
  int                arg_pos,
  const int          Zi,
  const int          Zj,
  const double       rij,
  const double       /*rij_sq*/,
  tadah::core::aed_type& dD_buf,
  const double       scale) const
{
  if (!is_init_for_atoms(Zi, Zj) || rij > get_rcut()) return false;

  if (theta_name == "CGRID2B") {
    if (arg_pos < 1 || static_cast<std::size_t>(arg_pos) > cgrid.size()) return false;
    const std::size_t k   = static_cast<std::size_t>(arg_pos - 1);
    const double      mu  = cgrid[k];
    const double      eta = sgrid[k];
    const double      xc  = eta * (rij - mu);
    const double      fc  = fcut->calc(rij);
    const double      dK_dx = (deg(k) == 3) ? dK3L(xc, 1.0)
                                            : dK5L(xc, 1.0);  // bare dK/dx_c
    // d(x_c)/d(mu) = -eta
    const double      dD_dmu = -eta * dK_dx * fc;
    dD_buf(fidx + k) += scale * weights[Zi][Zj] * dD_dmu;
    return true;
  }

  if (theta_name == "SGRID2B") {
    if (arg_pos < 1 || static_cast<std::size_t>(arg_pos) > sgrid.size()) return false;
    const std::size_t k   = static_cast<std::size_t>(arg_pos - 1);
    const double      mu  = cgrid[k];
    const double      eta = sgrid[k];
    const double      xc  = eta * (rij - mu);
    const double      fc  = fcut->calc(rij);
    const double      dK_dx = (deg(k) == 3) ? dK3L(xc, 1.0)
                                            : dK5L(xc, 1.0);
    // d(x_c)/d(eta) = (r - mu)
    const double      dD_deta = (rij - mu) * dK_dx * fc;
    dD_buf(fidx + k) += scale * weights[Zi][Zj] * dD_deta;
    return true;
  }

  if (theta_name == "RCUT2B") {
    // RCUT2B enters only through f_c(r).
    const double dfc = fcut->calc_drcut(rij);
    if (dfc == 0.0) return false;  // discontinuous cutoff (Cut_Dummy);
                                   // caller falls back to FD under HYBRID.
    for (std::size_t c = 0; c < cgrid.size(); ++c) {
      const double xc  = sgrid[c] * (rij - cgrid[c]);
      dD_buf(fidx + c) += scale * weights[Zi][Zj] * kval(c, xc, 1.0) * dfc;
    }
    return true;
  }

  if (theta_name == "WATOM") {
    if (!has_watom_arg_pos(arg_pos)) return false;
    const double chain = dweight_dw_n(Zi, Zj, arg_pos);
    if (chain != 0.0) {
      const double fc_ij = fcut->calc(rij);
      for (std::size_t c = 0; c < cgrid.size(); ++c) {
        const double xc = sgrid[c] * (rij - cgrid[c]);
        dD_buf(fidx + c) += scale * chain * kval(c, xc, fc_ij);
      }
    }
    return true;
  }

  return false;  // theta_name not owned by this descriptor.
}

// Closed-form mixed derivative d(fd_ij_radial)/d theta — same math as
// D2_Blip with (b, b1, b2) -> (K5L, K5L', K5L''). The family is C^4 so
// b2 = d2K5L is CONTINUOUS everywhere (no FD-probe hazards at knots).
bool D2_Knot5::calc_dXijdri_dtheta(
  const std::string& theta_name,
  int                arg_pos,
  const int          Zi,
  const int          Zj,
  const double       rij,
  const double       /*rij_sq*/,
  tadah::core::fd_type& dfd_ij,
  const double       scale) const
{
  if (!is_init_for_atoms(Zi, Zj) || rij > get_rcut()) return false;

  const double fc  = fcut->calc(rij);
  const double fcp = fcut->calc_prime(rij);

  if (theta_name == "CGRID2B") {
    if (arg_pos < 1 || static_cast<std::size_t>(arg_pos) > cgrid.size()) return false;
    const std::size_t k   = static_cast<std::size_t>(arg_pos - 1);
    const double      mu  = cgrid[k];
    const double      eta = sgrid[k];
    const double      xc  = eta * (rij - mu);
    const double      b1  = (deg(k) == 3) ? dK3L(xc, 1.0) : dK5L(xc, 1.0);
    const double      b2  = (deg(k) == 3) ? d2K3L(xc) : d2K5L(xc);
    const double      val = -eta * (fc * eta * b2 + fcp * b1);
    dfd_ij(fidx + k, 0) += scale * weights[Zi][Zj] * val;
    return true;
  }

  if (theta_name == "SGRID2B") {
    if (arg_pos < 1 || static_cast<std::size_t>(arg_pos) > sgrid.size()) return false;
    const std::size_t k   = static_cast<std::size_t>(arg_pos - 1);
    const double      mu  = cgrid[k];
    const double      eta = sgrid[k];
    const double      xc  = eta * (rij - mu);
    const double      b1  = (deg(k) == 3) ? dK3L(xc, 1.0) : dK5L(xc, 1.0);
    const double      b2  = (deg(k) == 3) ? d2K3L(xc) : d2K5L(xc);
    const double      val = fc * b1 + (rij - mu) * (fc * eta * b2 + fcp * b1);
    dfd_ij(fidx + k, 0) += scale * weights[Zi][Zj] * val;
    return true;
  }

  if (theta_name == "RCUT2B") {
    const double dfc  = fcut->calc_drcut(rij);
    const double dfcp = fcut->calc_prime_drcut(rij);
    if (dfc == 0.0 && dfcp == 0.0) return false;   // Cut_Dummy fallback
    for (std::size_t c = 0; c < cgrid.size(); ++c) {
      const double eta = sgrid[c];
      const double xc  = eta * (rij - cgrid[c]);
      const double b   = kval(c, xc, 1.0);
      const double b1  = (deg(c) == 3) ? dK3L(xc, 1.0) : dK5L(xc, 1.0);
      const double val = dfc * eta * b1 + dfcp * b;
      dfd_ij(fidx + c, 0) += scale * weights[Zi][Zj] * val;
    }
    return true;
  }

  if (theta_name == "WATOM") {
    if (!has_watom_arg_pos(arg_pos)) return false;
    const double chain = dweight_dw_n(Zi, Zj, arg_pos);
    if (chain != 0.0) {
      for (std::size_t c = 0; c < cgrid.size(); ++c) {
        const double xc = sgrid[c] * (rij - cgrid[c]);
        dfd_ij(fidx + c, 0) +=
          scale * chain * kder(c, xc, sgrid[c], fc, fcp);
      }
    }
    return true;
  }

  return false;
}

std::string D2_Knot5::label() const {
  return lab;
}
void D2_Knot5::init() {
  keys.push_back("CGRID2B");
  keys.push_back("SGRID2B");
  nparams=2;
}
bool D2_Knot5::reinit(const tadah::core::Context& ctx) {
  const size_t old_s = s;
  get_grid(ctx, "CGRID2B", cgrid);
  get_grid(ctx, "SGRID2B", sgrid);
  dgrid.clear();
  if (ctx.contains("DGRID2B") && ctx.size("DGRID2B") > 0)
    get_grid(ctx, "DGRID2B", dgrid);
  s = cgrid.size();
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUT2BMAX"));
  return s != old_s;
}

std::optional<DescriptorMeta> D2_Knot5::meta() const {
  using DataStats = tadah::core::stats::DataStats;
  static constexpr std::size_t I_P99 = tadah::core::stats::IDX_P99;

  DescriptorMeta dm;
  dm.label = "D2_Knot5";

  // RCUT2B — outer cutoff. Knots at or below rcut vanish smoothly there.
  ParamMeta p_rcut;
  p_rcut.key      = "RCUT2B";
  p_rcut.arg_pos  = -1;
  p_rcut.scale    = Scale::Distance;
  p_rcut.units    = "Å";
  p_rcut.role     = "Outer cutoff for the 2-body neighbour sum.";
  p_rcut.hard.lo  = 0.0;
  p_rcut.crit     = Criticality::Critical;
  p_rcut.scope    = PairScope::Global;
  dm.params.push_back(std::move(p_rcut));

  // CGRID2B — vector: knot positions (the only true shape parameters).
  ParamMeta p_cgrid;
  p_cgrid.key      = "CGRID2B";
  p_cgrid.arg_pos  = 0;
  p_cgrid.scale    = Scale::Distance;
  p_cgrid.units    = "Å";
  p_cgrid.role     = "Knot position r_s of the i-th truncated power "
                     "(r_s - r)_+^5; keep r_s <= RCUT2B.";
  p_cgrid.hard.lo  = 0.0;
  p_cgrid.crit     = Criticality::Recommended;
  p_cgrid.scope    = PairScope::PerElementPair;
  p_cgrid.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    return std::make_pair(0.0, pp.r_percentiles[I_P99]);
  };
  dm.params.push_back(std::move(p_cgrid));

  // SGRID2B — vector: a pure gauge for this family (conventionally 1.0).
  ParamMeta p_sgrid;
  p_sgrid.key      = "SGRID2B";
  p_sgrid.arg_pos  = 0;
  p_sgrid.scale    = Scale::Dimensionless;  // gauge factor, conventionally 1
  p_sgrid.units    = "";
  p_sgrid.role     = "Gauge scale of the i-th truncated power "
                     "((c,eta)->(c s^5, eta/s) is an exact invariance); "
                     "conventionally 1.0 — optimise CGRID2B and the "
                     "coefficients instead.";
  p_sgrid.hard.lo  = 1.0e-6;
  p_sgrid.crit     = Criticality::Optional;
  p_sgrid.scope    = PairScope::Global;
  dm.params.push_back(std::move(p_sgrid));

  dm.cross_param_checks.push_back(
    [](const DataStats&, const SeedView& s) -> HealthVerdict {
      if (!s.contains("CGRID2B") || !s.contains("SGRID2B")) return HealthOk{};
      if (s.size("CGRID2B") != s.size("SGRID2B"))
        return HealthFail{"CGRID2B and SGRID2B sizes differ."};
      return HealthOk{};
    });

  dm.expected_columns = [](const SeedView& s) -> std::size_t {
    if (!s.contains("CGRID2B")) return 0;
    return s.size("CGRID2B");
  };

  return dm;
}

} // namespace models
} // namespace tadah
