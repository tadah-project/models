#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_poly.h>
#include <tadah/models/descriptors/dm/dm_impl_utils.h>
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

  /* OrbitalWise + F≠F_SQ + L>0 is not rotationally invariant */
  if constexpr (std::is_same_v<AggMode, OrbitalWise>) {
    if (!ChannelFn::rotationally_invariant && Lmax > 0)
      core::Logger::instance().log(core::LogLevel::WARNING,
        label() + " (OrbitalWise): channel function is not F_SQ — descriptor"
        " is NOT rotationally invariant for L>0. Use PowerSpectrum mode.");
  }

  s        = cgrid.size() * (Lmax + 1);
  rhoisize = gen_atomic_orbitals(Lmax) * cgrid.size();
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
      /* ---- PowerSpectrum: F applied after invariant quadratic sum ---- */
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
    } else {
      /* ---- OrbitalWise: F applied per orbital before sum ------------ */
      std::size_t ii = ii_base;
      std::size_t dii = dii_base;
      for (std::size_t o = 0; o < N_orb; ++o) {
        const double f = fac[L][o];
        for (std::size_t c = 0; c < ng; ++c, ++ii, ++dii) {
          const double rho  = rhoi(ii);
          const double FF   = emb.calc_F (rho, c);
          const double dFF  = emb.calc_dF(rho, c);
          aed(c + Lshift + fidx) += f * FF;
          rhoi(dii)     = f * dFF;
        }
      }
    }

    Lshift  += ng;
    ii_base += N_orb * ng;
    dii_base += N_orb * ng;
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
  const double wi = weights[Zi];
  const double wj = weights[Zj];

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

      for (std::size_t c = 0; c < ng; ++c, ++ii) {
        const double g  = buf.g ()[c];
        const double dg = buf.dg()[c];

        const double Bx = g*txij + pqrij*dg*x*rinv;
        const double By = g*tyij + pqrij*dg*y*rinv;
        const double Bz = g*tzij + pqrij*dg*z*rinv;

        const double coef = scale *
            ( rhoi(ii)*wj - s_L * rhoj(ii)*wi );

        const std::size_t row = c + Lshift + fidx;
        fd_ij(row, 0) += coef * Bx;
        fd_ij(row, 1) += coef * By;
        fd_ij(row, 2) += coef * Bz;
      }
    }
    Lshift += ng;
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
  const double wj = weights[Zj];

  const double fc  = scale*wj * fcut->calc       (rij);
  const double fcp = scale*wj * fcut->calc_prime (rij);

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
    Lshift += ng;
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

  const double fc = scale * weights[Zj] * fcut->calc(rij);

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
  rhoi.resize(2*rhoisize);
  rhoi.set_zero();
}

template <typename ChannelFn, typename RadialBasis, typename AggMode>
std::string DM_Poly<ChannelFn, RadialBasis, AggMode>::label() const { return lab; }

template <typename ChannelFn, typename RadialBasis, typename AggMode>
void DM_Poly<ChannelFn, RadialBasis, AggMode>::init()
{
  keys = {"CGRIDMB", "SGRIDMB"};
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
  // TODO Probablky we do not need to do that:
  s        = cgrid.size() * static_cast<size_t>(Lmax + 1);
  rhoisize = gen_atomic_orbitals(Lmax) * cgrid.size();
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

} // namespace models
} // namespace tadah
