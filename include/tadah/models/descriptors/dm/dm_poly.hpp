#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_poly.h>
#include <tadah/models/descriptors/dm/dm_impl_utils.h>
#include <tadah/models/descriptors/d_basis_functions.h>

#include <stdexcept>

constexpr std::size_t DM_POLY_MAX_GRID = 64;
static_assert(DM_POLY_MAX_GRID >= 1, "DM_POLY_MAX_GRID must be > 0");

namespace {
using namespace tadah::models;
using namespace tadah::models::detail;

/* Blip radial pre-computation (will become a policy in Stage 2) */
inline void poly_precompute_radials(double r, double fc,
    const std::vector<double>& sg, const std::vector<double>& cg,
    double* g, std::size_t n)
{
  for (std::size_t i = 0; i < n; ++i) {
    const double x = sg[i]*(r - cg[i]);
    g[i] = B(x, fc);
  }
}
inline void poly_precompute_radials(double r, double fc, double fcp,
    const std::vector<double>& sg, const std::vector<double>& cg,
    double* g, double* dg, std::size_t n)
{
  for (std::size_t i = 0; i < n; ++i) {
    const double s = sg[i];
    const double x = s*(r - cg[i]);
    g [i] = B (x, fc);
    dg[i] = dB(x, s, fc, fcp);
  }
}
} // namespace

namespace tadah {
namespace models {

template <typename ChannelFn>
DM_Poly<ChannelFn>::DM_Poly() { init(); }

template <typename ChannelFn>
DM_Poly<ChannelFn>::DM_Poly(const tadah::core::Context& ctx)
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

  s        = cgrid.size() * (Lmax + 1);
  rhoisize = gen_atomic_orbitals(Lmax) * cgrid.size();

  init_for_atoms(get_init_atoms(ctx));
}

/* ================================================================== *
 *  calc_aed — power spectrum then channel function                   *
 *                                                                    *
 *  rhoi layout (each half of size rhoisize):                         *
 *    [0, rhoisize)          — orbital densities ρ[L,o,c]            *
 *    [rhoisize, 2*rhoisize) — per-orbital derivatives                *
 *                             ∂aed/∂ρ[L,o,c] = F'(ρ_L) · 2·fac·ρ  *
 * ================================================================== */
template <typename ChannelFn>
void DM_Poly<ChannelFn>::calc_aed(
    tadah::core::rho_type& rhoi,
    tadah::core::aed_type& aed) const
{
  const std::size_t ng = cgrid.size();
  std::size_t Lshift  = 0;
  std::size_t ii_base = rfidx;   // start of density part for current L

  for (int L = 0; L <= Lmax; ++L) {
    const std::size_t N_orb = orbitals[L].size();

    for (std::size_t c = 0; c < ng; ++c) {
      /* pass 1: power spectrum  ρ_L[c] = Σ_o fac[o] · ρ[o,c]² */
      double rho_L = 0.0;
      for (std::size_t o = 0; o < N_orb; ++o) {
        const double rho_oc = rhoi(ii_base + o*ng + c);
        rho_L += fac[L][o] * rho_oc * rho_oc;
      }

      /* apply channel function */
      const double FF  = emb.calc_F (rho_L, c);
      const double dFF = emb.calc_dF(rho_L, c);
      aed(c + Lshift + fidx) += FF;

      /* store per-orbital derivatives ∂aed/∂ρ[o,c] = dFF · 2·fac·ρ[o,c] */
      for (std::size_t o = 0; o < N_orb; ++o) {
        const double rho_oc = rhoi(ii_base + o*ng + c);
        rhoi(ii_base + o*ng + c + rhoisize) = dFF * 2.0 * fac[L][o] * rho_oc;
      }
    }

    Lshift  += ng;
    ii_base += N_orb * ng;
  }
}

/* ================================================================== *
 *  calc_dXijdri_dXjidri — half-list, symmetric pair contribution     *
 * ================================================================== */
template <typename ChannelFn>
void DM_Poly<ChannelFn>::calc_dXijdri_dXjidri(
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
  poly_precompute_radials(rij, fc, fcp, sgrid, cgrid, buf.g(), buf.dg(), ng);

  double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
  precompute_powers(x, Lmax, px, dpx);
  precompute_powers(y, Lmax, py, dpy);
  precompute_powers(z, Lmax, pz, dpz);

  const double   rinv   = 1.0/rij;
  std::size_t    Lshift = 0;
  std::size_t    ii     = rhoisize + rfidx;   // derivative part

  for (int L = 0; L <= Lmax; ++L) {
    const double s_L = (L & 1) ? 1.0 : -1.0; // (-1)^{L+1}

    for (std::size_t o = 0; o < orbitals[L].size(); ++o) {
      const int lx = orbitals[L][o][0],
                ly = orbitals[L][o][1],
                lz = orbitals[L][o][2];

      const double pxij  = px[lx], pyij = py[ly], pzij = pz[lz];
      const double pqrij = pxij*pyij*pzij;

      const double txij = dpx[lx]*pyij*pzij;
      const double tyij = dpy[ly]*pxij*pzij;
      const double tzij = dpz[lz]*pyij*pxij;

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
template <typename ChannelFn>
void DM_Poly<ChannelFn>::calc_dXijdri(
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

  const double       rinv = 1.0/rij;
  const std::size_t  ng   = cgrid.size();

  StackOrHeapBuffer<DM_POLY_MAX_GRID> buf(ng);
  poly_precompute_radials(rij, fc, fcp, sgrid, cgrid, buf.g(), buf.dg(), ng);

  double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
  precompute_powers(x, Lmax, px, dpx);
  precompute_powers(y, Lmax, py, dpy);
  precompute_powers(z, Lmax, pz, dpz);

  std::size_t Lshift = 0;
  std::size_t ii     = rhoisize + rfidx;

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
 *  calc_rho — density accumulation (identical to DM_Blip)            *
 * ================================================================== */
template <typename ChannelFn>
void DM_Poly<ChannelFn>::calc_rho(
    int Zi, int Zj, double rij, double,
    const tadah::core::Vec3d& v,
    tadah::core::rho_type& rhoi,
    double scale) const
{
  if (!is_init_for_atoms(Zi, Zj) || rij > get_rcut()) return;

  const double fc = scale * weights[Zj] * fcut->calc(rij);

  const std::size_t ng = cgrid.size();
  StackOrHeapBuffer<DM_POLY_MAX_GRID> buf(ng);
  poly_precompute_radials(rij, fc, sgrid, cgrid, buf.g(), ng);

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
template <typename ChannelFn>
void DM_Poly<ChannelFn>::init_rhoi(tadah::core::rho_type& rhoi)
{
  rhoi.resize(2*rhoisize);
  rhoi.set_zero();
}

template <typename ChannelFn>
std::string DM_Poly<ChannelFn>::label() const { return lab; }

template <typename ChannelFn>
void DM_Poly<ChannelFn>::init()
{
  keys = {"CGRIDMB", "SGRIDMB"};
  nparams = 2 + 1;
}

template <typename ChannelFn>
size_t DM_Poly<ChannelFn>::gen_atomic_orbitals(int Lm)
{
  orbitals.resize(Lm + 1);
  fac     .resize(Lm + 1);

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
