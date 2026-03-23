#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_ead.h>
#include <tadah/models/descriptors/d_basis_functions.h>

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

  init_for_atoms(get_init_atoms(ctx));
}

/* ------------------------------------------------------------------ *
 *  embedding energy + derivative store                               *
 * ------------------------------------------------------------------ */
void DM_EAD::calc_aed(core::rho_type& rhoi, core::aed_type& aed) const
{
  std::size_t Lshift = 0;
  std::size_t ii     = rfidx;          // start of density part

  for(int L=0; L<=Lmax; ++L){
    const std::size_t ng = cgrid.size();
    for(std::size_t o=0;o<orbitals[L].size();++o){
      const double f = fac[L][o];

      for(std::size_t c=0;c<ng;++c,++ii){
        const double rho =  rhoi(ii);
        aed (c+Lshift+fidx) += f * rho * rho;   // F(ρ)=ρ²
        rhoi(ii+rhoisize )  = 2.0 * f * rho;   // F′(ρ)=2ρ
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
  std::size_t   ii     = rhoisize + rfidx;      // derivative part

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
  std::size_t ii     = rhoisize + rfidx;  // derivative part of rho_i

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
  s        = cgrid.size() * static_cast<size_t>(Lmax + 1);
  rhoisize = gen_atomic_orbitals(Lmax) * cgrid.size();
  return s != old_s;
}

} // namespace models
} // namespace tadah

