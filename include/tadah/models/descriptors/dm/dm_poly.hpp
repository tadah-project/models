#ifndef DM_POLY_HPP
#define DM_POLY_HPP
#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_poly.h>
#include <tadah/models/descriptors/dm/dm_impl_utils.h>
#include <tadah/models/descriptors/warn_once.h>
#include <tadah/models/descriptors/dm/radial_policies.h>

#include <stdexcept>
#include <type_traits>

constexpr std::size_t DM_POLY_MAX_GRID = 64;
static_assert(DM_POLY_MAX_GRID >= 1, "DM_POLY_MAX_GRID must be > 0");

namespace {
using namespace tadah::models::detail;
} // namespace

namespace tadah {
namespace models {

template <typename ChannelFn, typename RadialBasis, typename AggMode>
DM_Poly<ChannelFn, RadialBasis, AggMode>::DM_Poly() { init(); }

template <typename ChannelFn, typename RadialBasis, typename AggMode>
DM_Poly<ChannelFn, RadialBasis, AggMode>::DM_Poly(const tadah::core::Context& ctx)
: DM_Base(ctx)
, rc(ctx.get<double>("RCUTMBMAX"))
{
  if (!ctx.get<bool>("INITMB")) return;
  init();

  get_grid(ctx, "CGRIDMB", cgrid);
  get_grid(ctx, "SGRIDMB", sgrid);
  D_Base::printGrid(cgrid, "CGRIDMB");
  D_Base::printGrid(sgrid, "SGRIDMB");

  if (cgrid.size() != sgrid.size())
    throw std::runtime_error("DM_Poly: CGRIDMB and SGRIDMB sizes differ.");

  emb = ChannelFn(ctx, cgrid.size(), sgrid.size());

  Lmax = ctx.get<int>("TYPEMB", 1);
  if (Lmax > 6)
    throw std::runtime_error("DM_Poly: Lmax > 6 is not supported (requested "
        + std::to_string(Lmax) + ").");

  /* OrbitalWise + F≠F_SQ + L>0 is not rotationally invariant; warned of
   * once per process (log_warning_once): tadah hpo builds the descriptors
   * of every evaluation. */
  if constexpr (std::is_same_v<AggMode, OrbitalWise>) {
    if (!ChannelFn::rotationally_invariant && Lmax > 0)
      log_warning_once(
        label() + " (OrbitalWise): channel not rotationally invariant for L>0",
        label() + " (OrbitalWise): channel function is not F_SQ — descriptor"
        " is NOT rotationally invariant for L>0. Use PowerSpectrum mode.");
  }
  odim     = emb.output_dim();
  const std::size_t ng = cgrid.size();
  s        = odim * (Lmax + 1);
  const std::size_t norbs = gen_atomic_orbitals(Lmax);
  rhoisize  = norbs * ng;      // ρ storage
  rhoipsize = norbs * odim;    // dF/dρ storage
  set_rfidx(0);
  set_drfidx(rhoisize);

  init_for_atoms(get_init_atoms(ctx));

}

/* ================================================================== *
 *  calc_aed                                                           *
 *                                                                     *
 *  PowerSpectrum mode:                                                *
 *    ρ_L[c] = Σ_o fac · ρ[o,c]²   then  V = F(ρ_L[c])              *
 *    derivative: ∂V/∂ρ[o,c] = F'(ρ_L) · 2·fac · ρ[o,c]             *
 *                                                                     *
 *  OrbitalWise mode:                                                  *
 *    V[L,c] = Σ_o fac · F(ρ[o,c])                                   *
 *    derivative: ∂V/∂ρ[o,c] = fac · F'(ρ[o,c])                      *
 *                                                                     *
 *  rhoi layout (each half of size rhoisize):                          *
 *    [0, rhoisize)          — orbital densities ρ[L,o,c]             *
 *    [rhoisize, 2*rhoisize) — per-orbital derivatives (see above)     *
 * ================================================================== */
template <typename ChannelFn, typename RadialBasis, typename AggMode>
void DM_Poly<ChannelFn, RadialBasis, AggMode>::calc_aed(
    tadah::core::rho_type& rhoi,
    tadah::core::aed_type& aed) const
{
  const std::size_t ng = cgrid.size();
  std::size_t Lshift  = 0;
  std::size_t ii_base = rfidx;
  std::size_t dii_base = drfidx;

  for (int L = 0; L <= Lmax; ++L) {
    const std::size_t N_orb = orbitals[L].size();

    if constexpr (std::is_same_v<AggMode, PowerSpectrum>) {

      if constexpr (ChannelFn::expands_dimensions) {
        /* PowerSpectrum + dimension-expanding (ng==1):
         * ρ_L = Σ_o fac · ρ[o]²  (single ρ per orbital)
         * then V[p] = F(ρ_L, p) for each embedding point p */
        double rho_L = 0.0;
        for (std::size_t o = 0; o < N_orb; ++o) {
          const double rho_o = rhoi(ii_base + o);
          rho_L += fac[L][o] * rho_o * rho_o;
        }
        for (std::size_t p = 0; p < odim; ++p) {
          const double FF  = emb.calc_F (rho_L, p);
          const double dFF = emb.calc_dF(rho_L, p);
          aed(p + Lshift + fidx) += FF;
          for (std::size_t o = 0; o < N_orb; ++o) {
            const double rho_o = rhoi(ii_base + o);
            rhoi(dii_base + o*odim + p) = dFF * 2.0 * fac[L][o] * rho_o;
          }
        }
      } else {
        /* PowerSpectrum + standard (odim == ng) */
        for (std::size_t c = 0; c < ng; ++c) {
          double rho_L = 0.0;
          for (std::size_t o = 0; o < N_orb; ++o) {
            const double rho_oc = rhoi(ii_base + o*ng + c);
            rho_L += fac[L][o] * rho_oc * rho_oc;
          }
          const double FF  = emb.calc_F (rho_L, c);
          const double dFF = emb.calc_dF(rho_L, c);
          aed(c + Lshift + fidx) += FF;
          for (std::size_t o = 0; o < N_orb; ++o) {
            const double rho_oc = rhoi(ii_base + o*ng + c);
            rhoi(dii_base + o*ng + c) = dFF * 2.0 * fac[L][o] * rho_oc;
          }
        }
      }

    } else {
      /* ---- OrbitalWise ---- */

      if constexpr (ChannelFn::expands_dimensions) {
        /* OrbitalWise + dimension-expanding (ng==1):
         * V[p] = Σ_o fac · F(ρ[o], p) */
        std::size_t ii = ii_base;
        std::size_t dii = dii_base;
        for (std::size_t o = 0; o < N_orb; ++o, ++ii) {
          const double f   = fac[L][o];
          const double rho = rhoi(ii);
          for (std::size_t p = 0; p < odim; ++p, ++dii) {
            const double FF  = emb.calc_F (rho, p);
            const double dFF = emb.calc_dF(rho, p);
            aed(p + Lshift + fidx) += f * FF;
            rhoi(dii) = f * dFF;
          }
        }
      } else {
        /* OrbitalWise + standard (odim == ng) */
        std::size_t ii = ii_base;
        std::size_t dii = dii_base;
        for (std::size_t o = 0; o < N_orb; ++o) {
          const double f = fac[L][o];
          for (std::size_t c = 0; c < ng; ++c, ++ii, ++dii) {
            const double rho  = rhoi(ii);
            const double FF   = emb.calc_F (rho, c);
            const double dFF  = emb.calc_dF(rho, c);
            aed(c + Lshift + fidx) += f * FF;
            rhoi(dii) = f * dFF;
          }
        }
      }
    }

    Lshift  += odim;
    ii_base += N_orb * ng;
    dii_base += N_orb * (ChannelFn::expands_dimensions ? odim : ng);
  }
}

/* ================================================================== *
 *  calc_dXijdri_dXjidri — half-list, symmetric pair contribution     *
 *  (identical for both AggModes — derivative layout is the same)     *
 * ================================================================== */
template <typename ChannelFn, typename RadialBasis, typename AggMode>
void DM_Poly<ChannelFn, RadialBasis, AggMode>::calc_dXijdri_dXjidri(
    int Zi, int Zj, double rij, double,
    const tadah::core::Vec3d& v,
    tadah::core::rho_type& rhoi,
    tadah::core::rho_type& rhoj,
    tadah::core::fd_type& fd_ij,
    double scale) const
{
  if (!is_init_for_atoms(Zi, Zj) || rij > get_rcut()) return;

  const double x = v[0], y = v[1], z = v[2];
  const double w_pair = weights[Zi][Zj];

  const double fc  = fcut->calc       (rij);
  const double fcp = fcut->calc_prime (rij);

  const std::size_t ng = cgrid.size();
  StackOrHeapBuffer<DM_POLY_MAX_GRID> buf(ng);
  RadialBasis::eval(rij, fc, fcp, sgrid, cgrid, buf.g(), buf.dg(), ng);

  double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
  precompute_powers(x, Lmax, px, dpx);
  precompute_powers(y, Lmax, py, dpy);
  precompute_powers(z, Lmax, pz, dpz);

  const double   rinv   = 1.0/rij;
  std::size_t    Lshift = 0;
  std::size_t    ii     = drfidx;

  for (int L = 0; L <= Lmax; ++L) {
    const double s_L = (L & 1) ? 1.0 : -1.0;

    for (std::size_t o = 0; o < orbitals[L].size(); ++o) {
      const int lx = orbitals[L][o][0],
                ly = orbitals[L][o][1],
                lz = orbitals[L][o][2];

      const double pxij  = px[lx], pyij = py[ly], pzij = pz[lz];
      const double pqrij = pxij*pyij*pzij;
      const double txij  = dpx[lx]*pyij*pzij;
      const double tyij  = dpy[ly]*pxij*pzij;
      const double tzij  = dpz[lz]*pyij*pxij;

      if constexpr (ChannelFn::expands_dimensions) {
        const double g  = buf.g ()[0];
        const double dg = buf.dg()[0];
        const double Bx = g*txij + pqrij*dg*x*rinv;
        const double By = g*tyij + pqrij*dg*y*rinv;
        const double Bz = g*tzij + pqrij*dg*z*rinv;

        for (std::size_t p = 0; p < odim; ++p, ++ii) {
          const double coef = scale * w_pair *
              ( rhoi(ii) - s_L * rhoj(ii) );
          const std::size_t row = p + Lshift + fidx;
          fd_ij(row, 0) += coef * Bx;
          fd_ij(row, 1) += coef * By;
          fd_ij(row, 2) += coef * Bz;
        }
      } else {
        for (std::size_t c = 0; c < ng; ++c, ++ii) {
          const double g  = buf.g ()[c];
          const double dg = buf.dg()[c];

          const double Bx = g*txij + pqrij*dg*x*rinv;
          const double By = g*tyij + pqrij*dg*y*rinv;
          const double Bz = g*tzij + pqrij*dg*z*rinv;

          const double coef = scale * w_pair *
              ( rhoi(ii) - s_L * rhoj(ii) );

          const std::size_t row = c + Lshift + fidx;
          fd_ij(row, 0) += coef * Bx;
          fd_ij(row, 1) += coef * By;
          fd_ij(row, 2) += coef * Bz;
        }
      }
    }
    Lshift += odim;
  }
}

/* ================================================================== *
 *  calc_dXijdri — full-list, one-sided                               *
 * ================================================================== */
template <typename ChannelFn, typename RadialBasis, typename AggMode>
void DM_Poly<ChannelFn, RadialBasis, AggMode>::calc_dXijdri(
    int Zi, int Zj, double rij, double,
    const tadah::core::Vec3d& v,
    tadah::core::rho_type& rhoi,
    tadah::core::fd_type& fd_ij,
    double scale) const
{
  if (!is_init_for_atoms(Zi, Zj) || rij > get_rcut()) return;

  const double x = v[0], y = v[1], z = v[2];
  const double w_pair = weights[Zi][Zj];

  const double fc  = scale*w_pair * fcut->calc       (rij);
  const double fcp = scale*w_pair * fcut->calc_prime (rij);

  const double      rinv = 1.0/rij;
  const std::size_t ng   = cgrid.size();

  StackOrHeapBuffer<DM_POLY_MAX_GRID> buf(ng);
  RadialBasis::eval(rij, fc, fcp, sgrid, cgrid, buf.g(), buf.dg(), ng);

  double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
  precompute_powers(x, Lmax, px, dpx);
  precompute_powers(y, Lmax, py, dpy);
  precompute_powers(z, Lmax, pz, dpz);

  std::size_t Lshift = 0;
  std::size_t ii     = drfidx;

  for (int L = 0; L <= Lmax; ++L) {
    for (std::size_t o = 0; o < orbitals[L].size(); ++o) {
      const int lx = orbitals[L][o][0],
                ly = orbitals[L][o][1],
                lz = orbitals[L][o][2];

      const double pxij  = px[lx], pyij = py[ly], pzij = pz[lz];
      const double pqrij = pxij*pyij*pzij;
      const double txij  = dpx[lx]*pyij*pzij;
      const double tyij  = dpy[ly]*pxij*pzij;
      const double tzij  = dpz[lz]*pyij*pxij;

      if constexpr (ChannelFn::expands_dimensions) {
        const double g  = buf.g ()[0];
        const double dg = buf.dg()[0];
        const double termx = g*txij + dg*pqrij*rinv*x;
        const double termy = g*tyij + dg*pqrij*rinv*y;
        const double termz = g*tzij + dg*pqrij*rinv*z;

        for (std::size_t p = 0; p < odim; ++p, ++ii) {
          const double f_r = rhoi(ii);
          const std::size_t row = p + Lshift + fidx;
          fd_ij(row, 0) += f_r * termx;
          fd_ij(row, 1) += f_r * termy;
          fd_ij(row, 2) += f_r * termz;
        }
      } else {
        const double t2x = pqrij * rinv * x;
        const double t2y = pqrij * rinv * y;
        const double t2z = pqrij * rinv * z;

#pragma omp simd
        for (std::size_t c = 0; c < ng; ++c) {
          const double g  = buf.g ()[c];
          const double dg = buf.dg()[c];

          const double termx = g*txij + dg*t2x;
          const double termy = g*tyij + dg*t2y;
          const double termz = g*tzij + dg*t2z;

          const double f_r = rhoi(ii + c);

          const std::size_t row = c + Lshift + fidx;
          fd_ij(row, 0) += f_r * termx;
          fd_ij(row, 1) += f_r * termy;
          fd_ij(row, 2) += f_r * termz;
        }
        ii += ng;
      }
    }
    Lshift += odim;
  }
}

/* ================================================================== *
 *  calc_rho — density accumulation (same for both AggModes)          *
 * ================================================================== */
template <typename ChannelFn, typename RadialBasis, typename AggMode>
void DM_Poly<ChannelFn, RadialBasis, AggMode>::calc_rho(
    int Zi, int Zj, double rij, double,
    const tadah::core::Vec3d& v,
    tadah::core::rho_type& rhoi,
    double scale) const
{
  if (!is_init_for_atoms(Zi, Zj) || rij > get_rcut()) return;

  const double fc = scale * weights[Zi][Zj] * fcut->calc(rij);

  const std::size_t ng = cgrid.size();
  StackOrHeapBuffer<DM_POLY_MAX_GRID> buf(ng);
  RadialBasis::eval(rij, fc, sgrid, cgrid, buf.g(), ng);

  double px[7], py[7], pz[7];
  precompute_powers(v[0], Lmax, px);
  precompute_powers(v[1], Lmax, py);
  precompute_powers(v[2], Lmax, pz);

  std::size_t ii = rfidx;

  for (int L = 0; L <= Lmax; ++L) {
    for (std::size_t o = 0; o < orbitals[L].size(); ++o) {
      const int lx = orbitals[L][o][0],
                ly = orbitals[L][o][1],
                lz = orbitals[L][o][2];

      const double poly = px[lx]*py[ly]*pz[lz];

      for (std::size_t c = 0; c < ng; ++c)
        rhoi(ii++) += buf.g()[c] * poly;
    }
  }
}

/* ================================================================== *
 *  Utilities                                                          *
 * ================================================================== */
template <typename ChannelFn, typename RadialBasis, typename AggMode>
void DM_Poly<ChannelFn, RadialBasis, AggMode>::init_rhoi(tadah::core::rho_type& rhoi)
{
  rhoi.resize(rhoisize+rhoipsize);
  rhoi.set_zero();
}

template <typename ChannelFn, typename RadialBasis, typename AggMode>
std::string DM_Poly<ChannelFn, RadialBasis, AggMode>::label() const { return lab; }

template <typename ChannelFn, typename RadialBasis, typename AggMode>
void DM_Poly<ChannelFn, RadialBasis, AggMode>::init()
{
  keys = {"CGRIDMB", "SGRIDMB"};
  context_keys = ChannelFn::context_keys();
  nparams = 2 + 1;
}

template <typename ChannelFn, typename RadialBasis, typename AggMode>
bool DM_Poly<ChannelFn, RadialBasis, AggMode>::reinit(const tadah::core::Context& ctx)
{
  const size_t old_s = s;
  get_grid(ctx, "CGRIDMB", cgrid);
  get_grid(ctx, "SGRIDMB", sgrid);
  Lmax = ctx.get<int>("TYPEMB", 1);
  // rc must be updated before gen_atomic_orbitals: fac[][] depends on rc
  rc = ctx.get<double>("RCUTMBMAX");
  if (fcut) fcut->set_rcut(rc);
  // Re-construct embedding function (re-reads its keys, e.g. CEMBFUNC/SEMBFUNC)
  emb = ChannelFn(ctx, cgrid.size(), sgrid.size());
  odim     = emb.output_dim();
  const std::size_t ng = cgrid.size();
  const std::size_t norbs = gen_atomic_orbitals(Lmax);
  s         = odim * static_cast<size_t>(Lmax + 1);
  rhoisize  = norbs * ng;
  rhoipsize = norbs * odim;
  set_drfidx(rhoisize);
  return s != old_s;
}

template <typename ChannelFn, typename RadialBasis, typename AggMode>
size_t DM_Poly<ChannelFn, RadialBasis, AggMode>::gen_atomic_orbitals(int Lm)
{
  orbitals.assign(Lm + 1, {});   // clear inner vectors so re-calls don't accumulate
  fac     .assign(Lm + 1, {});

  size_t cnt = 0;
  for (int L = 0; L <= Lm; ++L) {
    for (int lx = 0; lx <= L; ++lx) {
      for (int ly = 0; ly <= L-lx; ++ly) {
        const int lz = L - lx - ly;
        orbitals[L].push_back({lx, ly, lz});
        const double f = detail::fact_cached(L) /
            (detail::fact_cached(lx) * detail::fact_cached(ly) * detail::fact_cached(lz)) /
            detail::pow_i06(4*rc, L);
        fac[L].push_back(f);
        ++cnt;
      }
    }
  }
  return cnt;
}

/* ================================================================== *
 *  meta() — DM_Poly<F, R, M> shares the same parameter shape as
 *  DM_mEAD<F>: RCUTMB scalar, CGRIDMB / SGRIDMB paired vectors. The
 *  per-channel `F` may add CEMBFUNC/SEMBFUNC; those are picked up by
 *  the audit only when the seed declares them.
 * ================================================================== */
template <typename ChannelFn, typename RadialBasis, typename AggMode>
std::optional<DescriptorMeta>
DM_Poly<ChannelFn, RadialBasis, AggMode>::meta() const {
  using DataStats = tadah::core::stats::DataStats;
  static constexpr std::size_t I_P10 = tadah::core::stats::IDX_P10;
  static constexpr std::size_t I_P50 = tadah::core::stats::IDX_P50;
  static constexpr std::size_t I_P90 = tadah::core::stats::IDX_P90;
  static constexpr std::size_t I_P99 = tadah::core::stats::IDX_P99;

  DescriptorMeta dm;
  dm.label = label();   // runtime descriptor label

  ParamMeta p_rcut;
  p_rcut.key      = "RCUTMB";
  p_rcut.arg_pos  = -1;
  p_rcut.scale    = Scale::Distance;
  p_rcut.units    = "Å";
  p_rcut.role     = "Real-space cutoff for the MB neighbour sum.";
  p_rcut.hard.lo  = 0.0;
  p_rcut.crit     = Criticality::Critical;
  p_rcut.scope    = PairScope::Global;
  p_rcut.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    double low = pp.r_percentiles[I_P50];
    if (pp.r_peaks.size() >= 2) low = pp.r_peaks[1].position * 1.05;
    return std::make_pair(low, pp.r_percentiles[I_P99]);
  };
  p_rcut.check    = [](const DataStats& d, const SeedView& s) -> HealthVerdict {
    if (!s.contains("RCUTMB")) return HealthOk{};
    const double rc = s.get<double>("RCUTMB");
    if (d.per_pair.empty() && d.n_atoms > 0)
      return HealthFail{"RCUTMB=" + std::to_string(rc) +
                        " is too small — no MB neighbours found."};
    if (d.per_pair.empty()) return HealthOk{};
    if (rc < d.per_pair[0].r_min)
      return HealthFail{"RCUTMB=" + std::to_string(rc) +
                        " is below dataset r_min=" +
                        std::to_string(d.per_pair[0].r_min) + "."};
    return HealthOk{};
  };
  dm.params.push_back(std::move(p_rcut));

  ParamMeta p_cgrid;
  p_cgrid.key      = "CGRIDMB";
  p_cgrid.arg_pos  = 0;
  p_cgrid.scale    = Scale::Distance;
  p_cgrid.units    = "Å";
  p_cgrid.role     = "Radial centre. Must lie within dataset r_ij range.";
  p_cgrid.hard.lo  = 0.0;
  p_cgrid.crit     = Criticality::Critical;
  p_cgrid.scope    = PairScope::PerElementPair;
  p_cgrid.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    if (!pp.r_peaks.empty()) {
      const double mu = pp.r_peaks[0].position;
      const double w  = pp.r_peaks[0].fwhm > 0 ? pp.r_peaks[0].fwhm
                                               : 0.2 * mu;
      return std::make_pair(std::max(0.0, mu - w), mu + w);
    }
    return std::make_pair(pp.r_percentiles[I_P10],
                          pp.r_percentiles[I_P90]);
  };
  p_cgrid.check    = [](const DataStats& d, const SeedView& s) -> HealthVerdict {
    if (!s.contains("CGRIDMB")) return HealthOk{};
    if (d.per_pair.empty()) return HealthOk{};
    const std::size_t n = s.size("CGRIDMB");
    for (std::size_t i = 0; i < n; ++i) {
      double mu;
      try { mu = s.get<double>("CGRIDMB", i); } catch (...) { continue; }
      const auto& pp = d.per_pair[0];
      if (mu < pp.r_min || mu > pp.r_max)
        return HealthFail{"CGRIDMB[" + std::to_string(i) +
                          "]=" + std::to_string(mu) +
                          " is outside dataset r_ij range — ρ will collapse."};
    }
    return HealthOk{};
  };
  dm.params.push_back(std::move(p_cgrid));

  ParamMeta p_sgrid;
  p_sgrid.key      = "SGRIDMB";
  p_sgrid.arg_pos  = 0;
  p_sgrid.scale    = Scale::InverseDistanceSq;
  p_sgrid.units    = "Å^-2";
  p_sgrid.role     = "Inverse Gaussian width η.";
  p_sgrid.hard.lo  = 1.0e-6;
  p_sgrid.crit     = Criticality::Critical;
  p_sgrid.scope    = PairScope::Global;
  p_sgrid.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    const double w = pp.r_peaks.empty()
                       ? (pp.r_percentiles[I_P90] - pp.r_percentiles[I_P10])
                       : (pp.r_peaks[0].fwhm > 0 ? pp.r_peaks[0].fwhm : 0.5);
    return std::make_pair(0.25 / (w * w), 9.0 / (w * w));
  };
  dm.params.push_back(std::move(p_sgrid));

  dm.cross_param_checks.push_back(
    [](const DataStats&, const SeedView& s) -> HealthVerdict {
      if (!s.contains("CGRIDMB") || !s.contains("SGRIDMB")) return HealthOk{};
      if (s.size("CGRIDMB") != s.size("SGRIDMB"))
        return HealthFail{"CGRIDMB and SGRIDMB sizes differ."};
      return HealthOk{};
    });

  dm.expected_columns = [](const SeedView& s) -> std::size_t {
    if (!s.contains("CGRIDMB")) return 0;
    return s.size("CGRIDMB");
  };

  return dm;
}

} // namespace models
} // namespace tadah

#endif
