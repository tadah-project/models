#include "tadah/models/descriptors/d_base.h"
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <cstdlib>
#include <stdexcept>
namespace tadah {
namespace models {


D2_Blip::D2_Blip() {
  init();
}
D2_Blip::D2_Blip(const tadah::core::Context &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  init();

  get_grid(c,"CGRID2B",cgrid);
  get_grid(c,"SGRID2B",sgrid);
  D_Base::printGrid(cgrid,"CGRIDMB");
  D_Base::printGrid(sgrid,"SGRIDMB");

  if (cgrid.size()!=sgrid.size()) {
    throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
  }

  for (const auto &e:sgrid) {
    if(std::abs(e)<std::numeric_limits<double>::min())
      throw std::runtime_error("At least one of SGRID2B values is zero.\n");
  }

  s=cgrid.size();

  auto init_atoms = get_init_atoms(c);
  init_for_atoms(init_atoms);
}

void D2_Blip::calc_aed(
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
    double t = B(sgrid[c]*(rij-cgrid[c]),fc_ij);
    aed(i++) += scale*weights[Zj]*t;
  }

}
void D2_Blip::calc_dXijdri(
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
    fd_ij(i++,0) = scale*weights[Zj]*dB(sgrid[c]*(rij-cgrid[c]),sgrid[c],fc_ij,fcp_ij);
  }
}
void D2_Blip::calc_all(
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
    aed(i) += scale*weights[Zj]*B(sgrid[c]*(rij-cgrid[c]),fc_ij);
    fd_ij(i,0) = scale*weights[Zj]*dB(sgrid[c]*(rij-cgrid[c]),sgrid[c],fc_ij,fcp_ij);
    ++i;
  }
}

std::string D2_Blip::label() const {
  return lab;
}
void D2_Blip::init() {
  keys.push_back("CGRID2B");
  keys.push_back("SGRID2B");
  nparams=2;
}
bool D2_Blip::reinit(const tadah::core::Context& ctx) {
  const size_t old_s = s;
  get_grid(ctx, "CGRID2B", cgrid);
  get_grid(ctx, "SGRID2B", sgrid);
  s = cgrid.size();
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUT2BMAX"));
  return s != old_s;
}

std::optional<DescriptorMeta> D2_Blip::meta() const {
  using DataStats = tadah::core::stats::DataStats;
  static constexpr std::size_t I_P10 = tadah::core::stats::IDX_P10;
  static constexpr std::size_t I_P50 = tadah::core::stats::IDX_P50;
  static constexpr std::size_t I_P90 = tadah::core::stats::IDX_P90;
  static constexpr std::size_t I_P99 = tadah::core::stats::IDX_P99;

  DescriptorMeta dm;
  dm.label = "D2_Blip";

  // RCUT2B — outer cutoff in r for the 2-body design matrix.
  ParamMeta p_rcut;
  p_rcut.key      = "RCUT2B";
  p_rcut.arg_pos  = -1;
  p_rcut.scale    = Scale::Distance;
  p_rcut.units    = "Å";
  p_rcut.role     = "Outer cutoff for the 2-body neighbour sum.";
  p_rcut.hard.lo  = 0.0;            // r >= 0
  p_rcut.crit     = Criticality::Critical;
  p_rcut.scope    = PairScope::Global;
  p_rcut.suggest  = [](const DataStats& d, const SeedView&, int) {
    // Lower: include the second peak of the per-pair r_ij KDE if
    // present, else fall back to the median of r_ij. Upper: just past
    // p99 to preserve the long tail.
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    double low = pp.r_percentiles[I_P50];
    if (pp.r_peaks.size() >= 2) low = pp.r_peaks[1].position * 1.05;
    return std::make_pair(low, pp.r_percentiles[I_P99]);
  };
  p_rcut.check    = [](const DataStats& d, const SeedView& s) -> HealthVerdict {
    if (!s.contains("RCUT2B")) return HealthOk{};
    const double rc = s.get<double>("RCUT2B");
    if (d.per_pair.empty() && d.n_atoms > 0) {
      return HealthFail{"RCUT2B=" + std::to_string(rc) +
                        " is too small — no neighbours were found in any "
                        "structure (per_pair stats empty)."};
    }
    if (d.per_pair.empty()) return HealthOk{};
    if (rc < d.per_pair[0].r_min)
      return HealthFail{"RCUT2B=" + std::to_string(rc) +
                        " is below the dataset r_min=" +
                        std::to_string(d.per_pair[0].r_min) +
                        " — no 2-body neighbours."};
    return HealthOk{};
  };
  dm.params.push_back(std::move(p_rcut));

  // CGRID2B — vector: centre of i-th radial Gaussian.
  ParamMeta p_cgrid;
  p_cgrid.key      = "CGRID2B";
  p_cgrid.arg_pos  = 0;
  p_cgrid.scale    = Scale::Distance;
  p_cgrid.units    = "Å";
  p_cgrid.role     = "Centre of the i-th radial Gaussian (CGRID2B[i]).";
  p_cgrid.hard.lo  = 0.0;
  p_cgrid.crit     = Criticality::Recommended;
  p_cgrid.scope    = PairScope::PerElementPair;
  p_cgrid.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    // Spread CGRID2B[i] uniformly between 0 and the data's r_ij p99.
    return std::make_pair(0.0, pp.r_percentiles[I_P99]);
  };
  dm.params.push_back(std::move(p_cgrid));

  // SGRID2B — vector: inverse Gaussian width.
  ParamMeta p_sgrid;
  p_sgrid.key      = "SGRID2B";
  p_sgrid.arg_pos  = 0;
  p_sgrid.scale    = Scale::InverseDistanceSq;
  p_sgrid.units    = "Å^-2";
  p_sgrid.role     = "Inverse Gaussian width of the i-th radial.";
  // hard: SGRID2B = 0 makes the existing init() throw.
  p_sgrid.hard.lo  = 1.0e-6;
  p_sgrid.crit     = Criticality::Recommended;
  p_sgrid.scope    = PairScope::Global;
  p_sgrid.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    const double spread = pp.r_percentiles[I_P90] - pp.r_percentiles[I_P10];
    if (spread <= 0.0) return std::make_pair(0.1, 1.0);
    // η ∈ [0.25 / spread^2 , 25 / spread^2] — the spread^2 anchor
    // covers Gaussians from "broad enough to cover r_ij p10..p90" to
    // "thermal-scale narrow".
    return std::make_pair(0.25 / (spread * spread),
                          25.0 / (spread * spread));
  };
  dm.params.push_back(std::move(p_sgrid));

  // Cross-param: CGRID and SGRID must have matching sizes (mirrors the
  // existing init() runtime check at d2_blip.cpp:24).
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

