#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_blip.h>
#include <tadah/models/descriptors/dm/dm_impl_utils.h>
#include <tadah/models/descriptors/d_basis_functions.h>

#include <vector>
#include <stdexcept>

constexpr std::size_t MAX_GRID = 64;
static_assert(MAX_GRID>=1,"MAX_GRID must be > 0");

namespace {
using namespace tadah::models;
using namespace tadah::models::detail;

/* ------------ radial pre-computation (Blip) ----------------------- */
inline void precompute_radials(double r,double fc,
                               const std::vector<double>& sg,
                               const std::vector<double>& cg,
                               double* g,std::size_t n)
{
  for(std::size_t i=0;i<n;++i){
    const double x = sg[i]*(r-cg[i]);
    g[i] = B(x,fc);
  }
}
inline void precompute_radials(double r,double fc,double fcp,
                               const std::vector<double>& sg,
                               const std::vector<double>& cg,
                               double* g,double* dg,std::size_t n)
{
  for(std::size_t i=0;i<n;++i){
    const double s  = sg[i];
    const double x  = s*(r-cg[i]);
    g [i] = B (x,fc);
    dg[i] = dB(x,s,fc,fcp);          // derivative w.r.t. r
  }
}

} // namespace

/* ================================================================== *
 *  DM_Blip implementation                                            *
 * ================================================================== */
namespace tadah {
namespace models {

DM_Blip::DM_Blip() { init(); }

/* ------------------------------------------------------------------ */
DM_Blip::DM_Blip(const core::Context& ctx)
: DM_Base(ctx)
, rc(ctx.get<double>("RCUTMBMAX"))
{
  if(!ctx.get<bool>("INITMB")) return;
  init();

  get_grid(ctx,"CGRIDMB",cgrid);
  get_grid(ctx,"SGRIDMB",sgrid);
  D_Base::printGrid(cgrid,"CGRIDMB");
  D_Base::printGrid(sgrid,"SGRIDMB");

  if(cgrid.size()!=sgrid.size())
    throw std::runtime_error("CGRIDMB and SGRIDMB sizes differ.");

  Lmax = ctx.get<int>("TYPEMB",1);
  if(Lmax>6)
    throw std::runtime_error("DM_Blip: Lmax>6 not supported.");

  s         = cgrid.size()*(Lmax+1);
  rhoisize  = gen_atomic_orbitals(Lmax)*cgrid.size();

  init_for_atoms( get_init_atoms(ctx) );
}

/* ------------------ embedding energy / derivative ----------------- */
void DM_Blip::calc_aed(core::rho_type& rhoi, core::aed_type& aed) const
{
  std::size_t Lshift=0, ii=rfidx;
  const std::size_t ng=cgrid.size();

  for(int L=0;L<=Lmax;++L){
    for(std::size_t o=0;o<orbitals[L].size();++o){
      const double f=fac[L][o];
      for(std::size_t c=0;c<ng;++c,++ii){
        const double rho=rhoi(ii);
        aed (c+Lshift+fidx) += f*rho*rho;    // F(ρ)=ρ²
        rhoi(ii+rhoisize)   = 2.0*f*rho;     // F′(ρ)=2ρ
      }
    }
    Lshift+=ng;
  }
}
/* ------------------------------------------------------------------ *
 *  dX_ij / dr_i   &   dX_ji / dr_i  (pair contribution, symmetric)
 *  --- single-evaluation version (no duplicate j-branch) ------------
 * ------------------------------------------------------------------ */
void DM_Blip::calc_dXijdri_dXjidri(
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

/* -------------------- dXij/dr_i (one-sided) ----------------------- */
void DM_Blip::calc_dXijdri(
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

/* ----------------------- ρ_i accumulation ------------------------- */
void DM_Blip::calc_rho(int Zi,int Zj,double rij,double,
                       const core::Vec3d& v,
                       core::rho_type&    rhoi,
                       double             scale) const
{
  if(!is_init_for_atoms(Zi,Zj) || rij>get_rcut()) return;

  const double fc = scale*weights[Zj]*fcut->calc(rij);

  const std::size_t ng=cgrid.size();
  StackOrHeapBuffer<MAX_GRID> buf(ng);
  precompute_radials(rij,fc,sgrid,cgrid,buf.g(),ng);

  double px[7],py[7],pz[7];
  precompute_powers(v[0],Lmax,px);
  precompute_powers(v[1],Lmax,py);
  precompute_powers(v[2],Lmax,pz);

  std::size_t ii=rfidx;

  for(int L=0;L<=Lmax;++L){
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

/* ----------------------- rho_i memory ----------------------------- */
void DM_Blip::init_rhoi(core::rho_type& rhoi)
{
  rhoi.resize(2*rhoisize);
  rhoi.set_zero();
}

/* ------------------- orbital enumeration & factors ---------------- */
std::size_t DM_Blip::gen_atomic_orbitals(int Lm)
{
  orbitals.resize(Lm+1);
  fac     .resize(Lm+1);

  std::size_t cnt=0;
  for(int L=0;L<=Lm;++L){
    for(int lx=0;lx<=L;++lx){
      for(int ly=0;ly<=L-lx;++ly){
        int lz=L-lx-ly;
        orbitals[L].push_back({lx,ly,lz});
        const double f = fact_cached(L)/
                         (fact_cached(lx)*fact_cached(ly)*fact_cached(lz))/
                         pow_i06(4*rc,L);
        fac[L].push_back(f);
        ++cnt;
      }
    }
  }
  return cnt;
}

/* ------------------------------- misc ----------------------------- */
std::string DM_Blip::label() const { return lab; }

void DM_Blip::init()
{
  keys = {"CGRIDMB","SGRIDMB"};
  nparams = 2+1;
}

} // namespace models
} // namespace tadah

