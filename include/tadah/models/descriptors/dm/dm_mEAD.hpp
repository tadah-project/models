#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_mEAD.h>
#include <tadah/models/descriptors/dm/dm_impl_utils.h>
#include <tadah/models/descriptors/d_basis_functions.h>

constexpr std::size_t MAX_GRID = 64;
static_assert(MAX_GRID >= 1, "MAX_GRID must be > 0");

namespace {
using namespace tadah::models;
using namespace tadah::models::detail;

inline void precompute_radials(double rij,
      double fc,
      const std::vector<double>& sgrid,
      const std::vector<double>& cgrid,
      double* g, size_t ng)
  {
    for(std::size_t c=0; c<ng; ++c){
      g [c] = G (rij, sgrid[c], cgrid[c], fc );
    }
  }
inline void precompute_radials(double rij,
      double fc, double fcp,
      const std::vector<double>& sgrid,
      const std::vector<double>& cgrid,
      double* g, double* dg, size_t ng)
  {
    for(std::size_t c=0; c<ng; ++c){
      g [c] = G (rij, sgrid[c], cgrid[c], fc );
      dg[c] = dG(rij, sgrid[c], cgrid[c], fc , fcp);
    }
  }
}   // unnamed namespace

namespace tadah {
namespace models {

template <typename F>
  DM_mEAD<F>::DM_mEAD() {
    init();
  }

/* Index ii iterates rhoi array
 * rhoi array is split 50/50 into density and derivative
 * of the embedding function
 */
template <typename F>
  DM_mEAD<F>::DM_mEAD(const tadah::core::Context &context): DM_Base(context),
  rc(context.get<double>("RCUTMBMAX"))
{

  if (!context.get<bool>("INITMB")) return;
  init();


  get_grid(context,"CGRIDMB",cgrid);
  get_grid(context,"SGRIDMB",sgrid);

  D_Base::printGrid(cgrid,"CGRIDMB");
  D_Base::printGrid(sgrid,"SGRIDMB");


  core::Logger::instance().log(core::LogLevel::INFO,
      label()+" - No valid task found.");

  if (cgrid.size()!=sgrid.size()) {
    throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
  }

  emb = F(context, cgrid.size(), sgrid.size());

  Lmax = context.get<int>("TYPEMB",1);
  if (Lmax > 6)
    throw std::runtime_error(
        "DM_mEAD: Lmax > 6 is not supported (requested "
        + std::to_string(Lmax) + ").");

  if constexpr (!F::rotationally_invariant) {
    if (Lmax > 0)
      core::Logger::instance().log(core::LogLevel::WARNING,
        label() + ": channel function is not F_SQ — descriptor is NOT"
        " rotationally invariant for L>0. Use DM_PS_* descriptors instead.");
  }

  s=sgrid.size()*(Lmax+1);
  rhoisize = gen_atomic_orbitals(Lmax);
  rhoisize *= sgrid.size();

  auto init_atoms = get_init_atoms(context);
  init_for_atoms(init_atoms);
}
template <typename F>
  void DM_mEAD<F>::calc_aed(core::rho_type& rhoi, core::aed_type& aed) const
  {
    std::size_t Lshift = 0;
    std::size_t ii     = rfidx;

    for(int L=0; L<=Lmax; ++L){
      const std::size_t ng = cgrid.size();
      for(std::size_t o=0; o<orbitals[L].size(); ++o){
        const double f = fac[L][o];

        for(std::size_t c=0; c<ng; ++c, ++ii){
          const double rho = rhoi(ii);                // load once
          const double FF  = emb.calc_F (rho, c);
          const double dFF = emb.calc_dF(rho, c);

          aed (c+Lshift+fidx) += f * FF;
          rhoi(ii+rhoisize )  = f * dFF;
        }
      }
      Lshift += ng;
    }
  }
/* ------------------------------------------------------------------ *
 *  dX_ij/dr_i  and  dX_ji/dr_i   (pair contribution, symmetric)
 *  – corrected and de-duplicated version –                           *
 * ------------------------------------------------------------------ */
template <typename F>
void DM_mEAD<F>::calc_dXijdri_dXjidri(
        int Zi, int Zj, double rij, double,
        const core::Vec3d& v,
        core::rho_type&   rhoi,
        core::rho_type&   rhoj,
        core::fd_type&    fd_ij,
        double            scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

  /* --- geometry & constants --------------------------------------- */
  const double x=v[0],  y=v[1],  z=v[2];
  const double wi = weights[Zi];
  const double wj = weights[Zj];

  const double fc  = fcut->calc       (rij);
  const double fcp = fcut->calc_prime (rij);

  const double rinv = 1.0/rij;
  const std::size_t ng = cgrid.size();

  /* --- radial Gaussians ------------------------------------------- */
  StackOrHeapBuffer<MAX_GRID> buf(ng);
  precompute_radials(rij,fc,fcp,sgrid,cgrid,buf.g(),buf.dg(),ng);

  /* --- Cartesian powers and derivatives (up to L=6) ---------------- */
  double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
  precompute_powers(x,Lmax,px ,dpx);
  precompute_powers(y,Lmax,py ,dpy);
  precompute_powers(z,Lmax,pz ,dpz);

  /* --------------------------------------------------------------- */
  std::size_t Lshift = 0;                 // offset into fd_ij
  std::size_t ii     = rhoisize + rfidx;  // start of derivative part

  for (int L = 0; L <= Lmax; ++L)
  {
    const double s_L = (L & 1) ?  1.0   /* odd  L →  +1  */
                               : -1.0;  /* even L →  -1  = (-1)^{L+1} */

    for (std::size_t o = 0; o < orbitals[L].size(); ++o)
    {
      const int lx = orbitals[L][o][0],
                ly = orbitals[L][o][1],
                lz = orbitals[L][o][2];

      /* -------- i-partner Cartesian polynomial & derivatives ------- */
      const double pxij = px[lx],   pyij = py[ly],   pzij = pz[lz];
      const double pqrij = pxij*pyij*pzij;

      const double txij = dpx[lx]*pyij*pzij;
      const double tyij = dpy[ly]*pxij*pzij;
      const double tzij = dpz[lz]*pyij*pxij;

      for (std::size_t c = 0; c < ng; ++c, ++ii)
      {
        const double g  = buf.g ()[c];
        const double dg = buf.dg()[c];

        /* --- common Cartesian bracket (evaluated once) ------------- */
        const double Bx = g*txij + pqrij*dg*x*rinv;
        const double By = g*tyij + pqrij*dg*y*rinv;
        const double Bz = g*tzij + pqrij*dg*z*rinv;

        /* --- i/j prefactor ---------------------------------------- */
        const double coef = scale *
                            (  rhoi(ii)*wj          /* i-side        */
                            - s_L * rhoj(ii)*wi );  /* j-side, sign */

        const std::size_t row = c + Lshift + fidx;
        fd_ij(row,0) += coef * Bx;
        fd_ij(row,1) += coef * By;
        fd_ij(row,2) += coef * Bz;
      }
    }
    Lshift += ng;
  }
}

/* ================================================================== *
 *  one–sided descriptor derivative     dX_ij / dr_i                  *
 *  – streamlined, SIMD–friendly version                              *
 * ================================================================== */
template <typename F>
void DM_mEAD<F>::calc_dXijdri(
        int Zi,int Zj,double rij,double,
        const core::Vec3d& v,
        core::rho_type&    rhoi,
        core::fd_type&     fd_ij,
        double             scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

  /* ---- geometry & cut-off ---------------------------------------- */
  const double x = v[0],  y = v[1],  z = v[2];
  const double wj = weights[Zj];

  const double fc  = scale*wj * fcut->calc       (rij);
  const double fcp = scale*wj * fcut->calc_prime (rij);

  const double rinv = 1.0/rij;
  const std::size_t ng = cgrid.size();

  /* ---- radial Gaussians ------------------------------------------ */
  StackOrHeapBuffer<MAX_GRID> buf(ng);
  precompute_radials(rij,fc,fcp,sgrid,cgrid,buf.g(),buf.dg(),ng);

  /* ---- Cartesian powers & derivatives ---------------------------- */
  double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
  precompute_powers(x,Lmax,px ,dpx);
  precompute_powers(y,Lmax,py ,dpy);
  precompute_powers(z,Lmax,pz ,dpz);

  /* ---- main loops ------------------------------------------------ */
  std::size_t Lshift = 0;                 // row offset in fd_ij
  std::size_t ii     = rhoisize + rfidx;  // derivative part of ρ_i

  for (int L = 0; L <= Lmax; ++L)
  {
    const auto& orbs = orbitals[L];

    for (std::size_t o = 0; o < orbs.size(); ++o)
    {
      const int lx = orbs[o][0],
                ly = orbs[o][1],
                lz = orbs[o][2];

      /* --- polynomial and Cartesian derivatives (constant for grid) */
      const double pxij = px[lx],  pyij = py[ly],  pzij = pz[lz];

      const double pqrij = pxij*pyij*pzij;          //  x^{lx} y^{ly} z^{lz}

      const double txij = dpx[lx]*pyij*pzij;        // ∂P/∂x
      const double tyij = dpy[ly]*pxij*pzij;        // ∂P/∂y
      const double tzij = dpz[lz]*pyij*pxij;        // ∂P/∂z

      /* --- parts of the radial derivative that are constant in c --- */
      const double t2x_pref = pqrij * rinv * x;
      const double t2y_pref = pqrij * rinv * y;
      const double t2z_pref = pqrij * rinv * z;

      /* --- inner loop over radial grid – SIMD friendly ------------- */
#pragma omp simd
      for (std::size_t c = 0; c < ng; ++c)
      {
        const double g  = buf.g ()[c];       //   G(r)
        const double dg = buf.dg()[c];       // dG/dr

        const double termx = g*txij + dg*t2x_pref;
        const double termy = g*tyij + dg*t2y_pref;
        const double termz = g*tzij + dg*t2z_pref;

        const double f_r = rhoi(ii + c);     // ρ_i′ for that (o,c)

        const std::size_t row = c + Lshift + fidx;
        fd_ij(row,0) += f_r * termx;
        fd_ij(row,1) += f_r * termy;
        fd_ij(row,2) += f_r * termz;
      }
      ii += ng;                              // advance derivative index
    }
    Lshift += ng;
  }
}


template <typename F>
  std::string DM_mEAD<F>::label() const {
    return lab;
  }

template <typename F>
  void DM_mEAD<F>::init_rhoi(tadah::core::rho_type& rhoi)
  {
    rhoi.resize(2*rhoisize);
    rhoi.set_zero();
  }
template <typename F>
  void DM_mEAD<F>::calc_rho(
      int Zi,int Zj,double rij,double,const core::Vec3d& v,
      core::rho_type& rhoi,double scale) const
  {
    if(!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

    const double fc = scale*weights[Zj]*fcut->calc(rij);

    const std::size_t ng = cgrid.size();
    StackOrHeapBuffer<MAX_GRID> buf(ng);
    precompute_radials(rij, fc, sgrid, cgrid, buf.g(), ng);
    double px[7], py[7], pz[7];
    precompute_powers(v[0],Lmax,px);
    precompute_powers(v[1],Lmax,py);
    precompute_powers(v[2],Lmax,pz);

    std::size_t ii = rfidx;

    for(int L=0; L<=Lmax; ++L){
      const auto& orbs = orbitals[L];
      for(std::size_t o=0;o<orbs.size();++o){
        const int lx=orbs[o][0], ly=orbs[o][1], lz=orbs[o][2];
        const double poly = px[lx]*py[ly]*pz[lz];
        for(std::size_t c=0;c<ng;++c) {
          rhoi(ii++) += buf.g()[c]*poly;
        }
      }
    }
  }

template <typename F>
  size_t DM_mEAD<F>::gen_atomic_orbitals(int Lmax) {
    orbitals.assign(Lmax+1, {});   // clear inner vectors so re-calls don't accumulate
    fac.assign(Lmax+1, {});

    size_t count = 0;
    for (int L = 0; L <= Lmax; ++L) {
      for (int lx = 0; lx <= L; ++lx) {
        for (int ly = 0; ly <= L-lx; ++ly) {
          int lz = L - lx - ly;
          if (lx+ly+lz == L) {
            std::vector<int> o = {lx, ly, lz};
            orbitals[L].push_back(o);
            const double f = fact_cached(L)/(fact_cached(lx)*fact_cached(ly)*fact_cached(lz))/pow_i06(4*rc,L);
            fac[L].push_back(f);
            count++;
          }
        }
      }
    }
    return count;
  }
template <typename F>
  void DM_mEAD<F>::init() {
    keys.push_back("CGRIDMB");
    keys.push_back("SGRIDMB");
    keys.push_back("CEMBFUNC");
    keys.push_back("SEMBFUNC");
    nparams=5;
  }
template <typename F>
  bool DM_mEAD<F>::reinit(const tadah::core::Context& ctx) {
    const size_t old_s = s;
    get_grid(ctx, "CGRIDMB", cgrid);
    get_grid(ctx, "SGRIDMB", sgrid);
    Lmax = ctx.get<int>("TYPEMB", 1);
    // rc must be updated before gen_atomic_orbitals: fac[][] depends on rc
    rc = ctx.get<double>("RCUTMBMAX");
    if (fcut) fcut->set_rcut(rc);
    // Re-construct embedding function (re-reads CEMBFUNC/SEMBFUNC from ctx)
    emb = F(ctx, cgrid.size(), sgrid.size());
    s        = cgrid.size() * static_cast<size_t>(Lmax + 1);
    rhoisize = gen_atomic_orbitals(Lmax) * cgrid.size();
    return s != old_s;
  }

} // namespace models
} // namespace tadah
