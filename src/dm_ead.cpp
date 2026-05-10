#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_ead.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <tadah/models/descriptors/param_meta.h>
#include <tadah/core/stats/data_stats.h>

#include <array>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

constexpr std::size_t MAX_GRID = 64;      // tune to your input deck
static_assert(MAX_GRID>=1,"MAX_GRID must be > 0");

namespace {
using namespace tadah::models;

/* ------------------------------------------------------------------ *
 *  Fast, stack/heap scratch buffer (for G and dG/dr)                 *
 * ------------------------------------------------------------------ */
template<std::size_t N = 64>
class StackOrHeapBuffer {
public:
  explicit StackOrHeapBuffer(std::size_t n)
  : sz(n), use_heap(n>N), heap(use_heap?new double[2*n]:nullptr) {}

  double* g ()       { return use_heap?  heap.get()      : stack;     }
  double* dg()       { return use_heap? &heap.get()[sz] : &stack[N]; }
  std::size_t size() const { return sz; }

private:
  std::size_t              sz;
  bool                     use_heap;
  std::unique_ptr<double[]> heap;
  alignas(64) double       stack[2*N];
};

/* -------- factorials & small integer powers (0…6) ----------------- */
constexpr std::array<double,7> FACT_TBL = {1.,1.,2.,6.,24.,120.,720.};

inline double fact_cached(int n)      { return FACT_TBL[n]; }

inline double pow_i06(double x,int n)
{
  switch(n){
    case 0: return 1.;
    case 1: return x;
    case 2: return x*x;
    case 3: return x*x*x;
    case 4: { double x2=x*x; return x2*x2; }
    case 5: { double x2=x*x; return x2*x2*x; }
    case 6: { double x2=x*x; return x2*x2*x2; }
    default:return std::pow(x,n);
  }
}

/* -------- Cartesian powers and their derivatives ------------------ */
inline void precompute_powers(double x,int L,double p[7])
{
  p[0]=1.; for(int n=1;n<=L;++n) p[n]=p[n-1]*x;
#if defined(__SANITIZE_ADDRESS__)||defined(ADDRESS_SANITIZER)
  for(int n=L+1;n<=6;++n) p[n]=0.;
#endif
}
inline void precompute_powers(double x,int L,double p[7],double dp[7])
{
  p[0]=1.; dp[0]=0.;
  for(int n=1;n<=L;++n){ p[n]=p[n-1]*x; dp[n]=n*p[n-1]; }
#if defined(__SANITIZE_ADDRESS__)||defined(ADDRESS_SANITIZER)
  for(int n=L+1;n<=6;++n){ p[n]=0.; dp[n]=0.; }
#endif
}

/* -------- pre-compute Gaussian (and derivative) for a whole grid --- */
inline void precompute_radials(double r,double fc,
                               const std::vector<double>& sg,
                               const std::vector<double>& cg,
                               double *g,std::size_t n)
{
  for(std::size_t i=0;i<n;++i) g[i]=G(r,sg[i],cg[i],fc);
}
inline void precompute_radials(double r,double fc,double fcp,
                               const std::vector<double>& sg,
                               const std::vector<double>& cg,
                               double *g,double* dg,std::size_t n)
{
  for(std::size_t i=0;i<n;++i){
    g [i]=G (r,sg[i],cg[i],fc );
    dg[i]=dG(r,sg[i],cg[i],fc,fcp);
  }
}

/* sign-flip without run-time branch */
inline double neg(double v,int n){ return (n&1)?-v:v; }

} // unnamed namespace


/* ================================================================== *
 *  DM_EAD  (ρ² embedding)                                           *
 * ================================================================== */
namespace tadah {
namespace models {
    

DM_EAD::DM_EAD() { init(); }

/* ------------------------------------------------------------------ */
DM_EAD::DM_EAD(const core::Context& ctx)
: DM_Base(ctx)
, rc(ctx.get<double>("RCUTMBMAX"))
{
  if(!ctx.get<bool>("INITMB"))
    throw std::runtime_error("INITMB=false; cannot initialise "+label()+"\n");

  init();

  get_grid(ctx,"CGRIDMB",cgrid);
  get_grid(ctx,"SGRIDMB",sgrid);
  D_Base::printGrid(cgrid,"CGRIDMB");
  D_Base::printGrid(sgrid,"SGRIDMB");

  if(cgrid.size()!=sgrid.size())
    throw std::runtime_error("SGRIDMB and CGRIDMB array sizes differ!");

  Lmax = ctx.get<int>("TYPEMB",1);
  if(Lmax>6)
    throw std::runtime_error("DM_EAD: Lmax>6 not supported (requested "+
                             std::to_string(Lmax)+").");

  s         = cgrid.size()*(Lmax+1);
  rhoisize  = gen_atomic_orbitals(Lmax)*cgrid.size();
  set_rfidx(0);
  set_drfidx(rhoisize);

  init_for_atoms(get_init_atoms(ctx));
}

/* ------------------------------------------------------------------ *
 *  embedding energy + derivative store                               *
 * ------------------------------------------------------------------ */
void DM_EAD::calc_aed(core::rho_type& rhoi, core::aed_type& aed) const
{
  std::size_t Lshift = 0;
  std::size_t ii     = rfidx;          // start of density part
  std::size_t dii    = drfidx;         // start of derivative part

  for(int L=0; L<=Lmax; ++L){
    const std::size_t ng = cgrid.size();
    for(std::size_t o=0;o<orbitals[L].size();++o){
      const double f = fac[L][o];

      for(std::size_t c=0;c<ng;++c,++ii,++dii){
        const double rho =  rhoi(ii);
        aed (c+Lshift+fidx) += f * rho * rho;   // F(ρ)=ρ²
        rhoi(dii)  = 2.0 * f * rho;   // F′(ρ)=2ρ
      }
    }
    Lshift += ng;
  }
}

/* ------------------------------------------------------------------ *
 *  dX_ij / dr_i   &   dX_ji / dr_i  (pair contribution, symmetric)
 *  --- single-evaluation version (no duplicate j-branch) ------------
 * ------------------------------------------------------------------ */
void DM_EAD::calc_dXijdri_dXjidri(
        int Zi,int Zj,double rij,double,
        const core::Vec3d& v,
        core::rho_type&  rhoi,
        core::rho_type&  rhoj,
        core::fd_type&   fd_ij,
        double           scale) const
{
  if(!is_init_for_atoms(Zi,Zj) || rij>get_rcut()) return;

  const double x=v[0], y=v[1], z=v[2];
  const double wi = weights[Zi];
  const double wj = weights[Zj];

  /* --- radial part ------------------------------------------------- */
  const double fc  = fcut->calc       (rij);
  const double fcp = fcut->calc_prime (rij);

  const std::size_t ng = cgrid.size();
  StackOrHeapBuffer<MAX_GRID> buf(ng);
  precompute_radials(rij,fc,fcp,sgrid,cgrid,buf.g(),buf.dg(),ng);

  /* --- Cartesian powers (and derivatives) -------------------------- */
  double px[7],py[7],pz[7],dpx[7],dpy[7],dpz[7];
  precompute_powers(x,Lmax,px,dpx);
  precompute_powers(y,Lmax,py,dpy);
  precompute_powers(z,Lmax,pz,dpz);

  const double rinv = 1.0/rij;
  std::size_t   Lshift = 0;
  std::size_t   ii     = drfidx;      // derivative part

  for(int L=0; L<=Lmax; ++L)
  {
    const double s_L = (L & 1) ? 1.0 : -1.0;    //  (-1)^{L+1}

    for(std::size_t o=0; o<orbitals[L].size(); ++o)
    {
      const int lx = orbitals[L][o][0],
                ly = orbitals[L][o][1],
                lz = orbitals[L][o][2];

      /* polynomials and their Cartesian derivatives (partner i) ----- */
      const double pxij = px[lx],   pyij = py[ly],   pzij = pz[lz];
      const double pqrij = pxij*pyij*pzij;

      const double txij = dpx[lx]*pyij*pzij;
      const double tyij = dpy[ly]*pxij*pzij;
      const double tzij = dpz[lz]*pyij*pxij;

      for(std::size_t c=0; c<ng; ++c, ++ii)
      {
        const double g  = buf.g ()[c];
        const double dg = buf.dg()[c];

        /* common Cartesian brackets --------------------------------- */
        const double Bx = g*txij + pqrij*dg*x*rinv;
        const double By = g*tyij + pqrij*dg*y*rinv;
        const double Bz = g*tzij + pqrij*dg*z*rinv;

        /* i/j prefactor --------------------------------------------- */
        const double coef = scale *
                            ( rhoi(ii)*wj      /* i-side  */
                            - s_L * rhoj(ii)*wi ); /* j-side */

        const std::size_t row = c + Lshift + fidx;
        fd_ij(row,0) += coef * Bx;
        fd_ij(row,1) += coef * By;
        fd_ij(row,2) += coef * Bz;
      }
    }
    Lshift += ng;
  }
}

/* ------------------------------------------------------------------ *
 *  dX_ij/dr_i  (only one-sided, used e.g. for forces on atom i)      *
 * ------------------------------------------------------------------ */
void DM_EAD::calc_dXijdri(
        int Zi,int Zj,double rij,double,
        const core::Vec3d& v,
        core::rho_type&   rhoi,
        core::fd_type&    fd_ij,
        double            scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

  const double x = v[0],  y = v[1],  z = v[2];
  const double wj = weights[Zj];

  const double fc  = scale*wj * fcut->calc       (rij);
  const double fcp = scale*wj * fcut->calc_prime (rij);

  const double rinv = 1.0/rij;
  const std::size_t ng = cgrid.size();

  /* radial arrays -------------------------------------------------- */
  StackOrHeapBuffer<MAX_GRID> buf(ng);
  precompute_radials(rij,fc,fcp,sgrid,cgrid,buf.g(),buf.dg(),ng);

  /* Cartesian powers & derivatives --------------------------------- */
  double px[7],py[7],pz[7],dpx[7],dpy[7],dpz[7];
  precompute_powers(x,Lmax,px ,dpx);
  precompute_powers(y,Lmax,py ,dpy);
  precompute_powers(z,Lmax,pz ,dpz);

  std::size_t Lshift = 0;                 // row offset in fd
  std::size_t ii     = drfidx;  // derivative part of rho_i

  /* ========================== main loops ========================== */
  for (int L = 0; L <= Lmax; ++L)
  {
    for (std::size_t o = 0; o < orbitals[L].size(); ++o)
    {
      const int lx = orbitals[L][o][0],
                ly = orbitals[L][o][1],
                lz = orbitals[L][o][2];

      const double pxij = px[lx],  pyij = py[ly],  pzij = pz[lz];

      const double pqrij = pxij*pyij*pzij;          // x^lx y^ly z^lz

      const double txij = dpx[lx]*pyij*pzij;        // ∂/∂x
      const double tyij = dpy[ly]*pxij*pzij;        // ∂/∂y
      const double tzij = dpz[lz]*pyij*pxij;        // ∂/∂z

      /* coord-dependent part of radial derivative (constant for c) -- */
      const double t2x_pref = pqrij * rinv * x;
      const double t2y_pref = pqrij * rinv * y;
      const double t2z_pref = pqrij * rinv * z;

      /* ----------------------- grid loop --------------------------- */
#pragma omp simd
      for (std::size_t c = 0; c < ng; ++c)
      {
        const double g  = buf.g ()[c];
        const double dg = buf.dg()[c];

        const double termx = g*txij + dg*t2x_pref;
        const double termy = g*tyij + dg*t2y_pref;
        const double termz = g*tzij + dg*t2z_pref;

        const double f_r = rhoi(ii + c);            // ρ′_i for (o,c)

        const std::size_t row = c + Lshift + fidx;
        fd_ij(row,0) += f_r * termx;
        fd_ij(row,1) += f_r * termy;
        fd_ij(row,2) += f_r * termz;
      }
      ii += ng;                                     // advance derivative index
    }
    Lshift += ng;
  }
}

/* ------------------------------------------------------------------ *
 *  neighbour-density accumulation                                    *
 * ------------------------------------------------------------------ */
void DM_EAD::calc_rho(int Zi,int Zj,double rij,double,
                      const core::Vec3d& v,
                      core::rho_type&    rhoi,
                      double             scale) const
{
  if(!is_init_for_atoms(Zi,Zj) || rij>get_rcut()) return;

  const double fc = scale*weights[Zj]*fcut->calc(rij);

  const std::size_t ng = cgrid.size();
  StackOrHeapBuffer<MAX_GRID> buf(ng);
  precompute_radials(rij,fc,sgrid,cgrid,buf.g(),ng);

  double px[7],py[7],pz[7];
  precompute_powers(v[0],Lmax,px);
  precompute_powers(v[1],Lmax,py);
  precompute_powers(v[2],Lmax,pz);

  std::size_t ii = rfidx;

  for(int L=0; L<=Lmax; ++L){
    for(std::size_t o=0;o<orbitals[L].size();++o){
      const int lx=orbitals[L][o][0],
                ly=orbitals[L][o][1],
                lz=orbitals[L][o][2];

      const double poly = px[lx]*py[ly]*pz[lz];

      for(std::size_t c=0;c<ng;++c)
        rhoi(ii++) += buf.g()[c]*poly;
    }
  }
}

/* ------------------------------------------------------------------ *
 *  allocate & zero rho_i                                             *
 * ------------------------------------------------------------------ */
void DM_EAD::init_rhoi(core::rho_type& rhoi)
{
  rhoi.resize(2*rhoisize);
  rhoi.set_zero();
}

/* ------------------------------------------------------------------ *
 *  multinomial prefactors & orbital list (unique, L≤6)               *
 * ------------------------------------------------------------------ */
std::size_t DM_EAD::gen_atomic_orbitals(int Lmax)
{
  orbitals.assign(Lmax+1, {});   // clear inner vectors so re-calls don't accumulate
  fac     .assign(Lmax+1, {});

  std::size_t count=0;
  for(int L=0;L<=Lmax;++L){
    for(int lx=0;lx<=L;++lx){
      for(int ly=0;ly<=L-lx;++ly){
        int lz=L-lx-ly;
        std::vector<int> orb{lx,ly,lz};
        orbitals[L].push_back(orb);

        const double f = fact_cached(L) /
                         (fact_cached(lx)*fact_cached(ly)*fact_cached(lz)) /
                         pow_i06(4*rc,L);
        fac[L].push_back(f);
        ++count;
      }
    }
  }
  return count;
}

/* ------------------------------------------------------------------ */
std::string DM_EAD::label() const { return lab; }

/* ------------------------------------------------------------------ */
void DM_EAD::init()
{
  keys = {"CGRIDMB","SGRIDMB"};
  nparams = 2+1;  // (two grids + r_c)
}
bool DM_EAD::reinit(const tadah::core::Context& ctx) {
  const size_t old_s = s;
  get_grid(ctx, "CGRIDMB", cgrid);
  get_grid(ctx, "SGRIDMB", sgrid);
  Lmax = ctx.get<int>("TYPEMB", 1);
  // rc must be updated before gen_atomic_orbitals: fac[][] depends on rc
  rc = ctx.get<double>("RCUTMBMAX");
  if (fcut) fcut->set_rcut(rc);
  // TODO Probablky we do not need to do that:
  s        = cgrid.size() * static_cast<size_t>(Lmax + 1);
  rhoisize = gen_atomic_orbitals(Lmax) * cgrid.size();
  return s != old_s;
}

// ----------------------------------------------------------------
// meta() — DM_EAD shares the (RCUTMB, CGRIDMB, SGRIDMB) parameter
// shape with the BASINS-fix DM_mBlipEnv: per-pair radial GTO with
// centre/width vectors. The ρ-collapse cross-param check applies to
// the same per-descriptor key (`DM_EAD`) — the suggester anchors all
// numbers to the dataset's r_ij percentiles.
// ----------------------------------------------------------------
std::optional<DescriptorMeta> DM_EAD::meta() const {
  using DataStats = tadah::core::stats::DataStats;
  static constexpr std::size_t I_P10 = tadah::core::stats::IDX_P10;
  static constexpr std::size_t I_P50 = tadah::core::stats::IDX_P50;
  static constexpr std::size_t I_P90 = tadah::core::stats::IDX_P90;
  static constexpr std::size_t I_P99 = tadah::core::stats::IDX_P99;

  DescriptorMeta dm;
  dm.label = "DM_EAD";

  // RCUTMB — real-space cutoff for the MB neighbour sum.
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
    if (d.per_pair.empty()) return HealthOk{};
    const double rc = s.get<double>("RCUTMB");
    if (rc < d.per_pair[0].r_min)
      return HealthFail{"RCUTMB is below the smallest r_ij — no MB neighbours."};
    return HealthOk{};
  };
  dm.params.push_back(std::move(p_rcut));

  // CGRIDMB — vector of GTO centres (r_s in the EAD formula).
  ParamMeta p_cgrid;
  p_cgrid.key      = "CGRIDMB";
  p_cgrid.arg_pos  = 0;
  p_cgrid.scale    = Scale::Distance;
  p_cgrid.units    = "Å";
  p_cgrid.role     = "GTO centre r_s. Must lie within the dataset r_ij range or ρ collapses.";
  p_cgrid.hard.lo  = 0.0;
  p_cgrid.crit     = Criticality::Critical;
  p_cgrid.scope    = PairScope::PerElementPair;
  p_cgrid.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    if (!pp.r_peaks.empty()) {
      const double mu = pp.r_peaks[0].position;
      const double w  = pp.r_peaks[0].fwhm > 0
                          ? pp.r_peaks[0].fwhm
                          : 0.2 * mu;
      return std::make_pair(std::max(0.0, mu - w), mu + w);
    }
    return std::make_pair(pp.r_percentiles[I_P10],
                          pp.r_percentiles[I_P90]);
  };
  p_cgrid.check    = [](const DataStats& d, const SeedView& s) -> HealthVerdict {
    if (!s.contains("CGRIDMB")) return HealthOk{};
    if (d.per_pair.empty()) return HealthOk{};
    // Per-arg-pos check: any vector entry outside the data's r_ij range collapses ρ.
    const std::size_t n = s.size("CGRIDMB");
    for (std::size_t i = 0; i < n; ++i) {
      double mu;
      try { mu = s.get<double>("CGRIDMB", i); } catch (...) { continue; }
      const auto& pp = d.per_pair[0];
      if (mu < pp.r_min || mu > pp.r_max)
        return HealthFail{"CGRIDMB[" + std::to_string(i) +
                          "]=" + std::to_string(mu) +
                          " is outside the dataset's r_ij range — ρ will be ~0."};
    }
    return HealthOk{};
  };
  dm.params.push_back(std::move(p_cgrid));

  // SGRIDMB — vector of GTO inverse widths (η).
  ParamMeta p_sgrid;
  p_sgrid.key      = "SGRIDMB";
  p_sgrid.arg_pos  = 0;
  p_sgrid.scale    = Scale::InverseDistanceSq;
  p_sgrid.units    = "Å^-2";
  p_sgrid.role     = "Inverse Gaussian width η. Narrow GTOs combined with off-shell centres collapse ρ to 0.";
  p_sgrid.hard.lo  = 1.0e-6;
  p_sgrid.crit     = Criticality::Critical;
  p_sgrid.scope    = PairScope::Global;
  p_sgrid.suggest  = [](const DataStats& d, const SeedView&, int) {
    if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
    const auto& pp = d.per_pair[0];
    const double w = pp.r_peaks.empty()
                       ? (pp.r_percentiles[I_P90] - pp.r_percentiles[I_P10])
                       : (pp.r_peaks[0].fwhm > 0 ? pp.r_peaks[0].fwhm
                                                 : 0.5);
    return std::make_pair(0.25 / (w * w), 9.0 / (w * w));
  };
  dm.params.push_back(std::move(p_sgrid));

  // Cross-param: CGRIDMB and SGRIDMB sizes must agree.
  dm.cross_param_checks.push_back(
    [](const DataStats&, const SeedView& s) -> HealthVerdict {
      if (!s.contains("CGRIDMB") || !s.contains("SGRIDMB")) return HealthOk{};
      if (s.size("CGRIDMB") != s.size("SGRIDMB"))
        return HealthFail{"CGRIDMB and SGRIDMB sizes differ."};
      return HealthOk{};
    });

  // expected_columns: |CGRIDMB| × (Lmax + 1).
  dm.expected_columns = [](const SeedView& s) -> std::size_t {
    if (!s.contains("CGRIDMB")) return 0;
    const std::size_t n = s.size("CGRIDMB");
    int Lmax = 0;
    if (s.contains("TYPEMB")) {
      try { Lmax = s.context().get<int>("TYPEMB", 1); } catch (...) {}
    }
    return n * static_cast<std::size_t>(Lmax + 1);
  };

  return dm;
}

} // namespace models
} // namespace tadah

