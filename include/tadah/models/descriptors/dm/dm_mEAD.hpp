#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_mEAD.h>
#include <tadah/models/descriptors/d_basis_functions.h>

#include <array>
#include <memory>

constexpr std::size_t MAX_GRID = 64;          // <-- adjust to your input deck
static_assert(MAX_GRID >= 1, "MAX_GRID must be > 0");



namespace {
using namespace tadah::models;

template<std::size_t N = 64>
class StackOrHeapBuffer
{
  public:
    explicit StackOrHeapBuffer(std::size_t n)
      : sz(n)
        , use_heap(n > N)
        , heap(use_heap ? new double[2*n] : nullptr)   // ×2 for g and dg
        {}

    /* first half: G(r)   second half: dG/dr (if needed) */
    double* g ()       { return use_heap ?  heap.get()      : stack; }
    double* dg()       { return use_heap ? &heap.get()[sz] : &stack[N]; }
    std::size_t size() const { return sz; }

  private:
    std::size_t  sz;
    bool         use_heap;
    std::unique_ptr<double[]> heap;
    alignas(64) double stack[2*N];    // on-stack contiguous block
};

/* Factorials 0! … 6! at compile time */
constexpr std::array<double,7> FACT_TBL =
{1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0};

inline double fact_cached(int n)
{
  return FACT_TBL[n];          // n is guaranteed 0…6
}

inline double pow_i06(double x, int n)
{
  switch (n) {
    case 0:  return 1.0;
    case 1:  return x;
    case 2:  return x * x;
    case 3:  return x * x * x;
    case 4:  { double x2 = x * x;  return x2 * x2; }                // x⁴
    case 5:  { double x2 = x * x;  return x2 * x2 * x; }            // x⁵
    case 6:  { double x2 = x * x;  return x2 * x2 * x2; }           // x⁶
    default: return std::pow(x, n);  // should never hit (L ≤ 6)
  }
}
inline void precompute_powers(double x, int L, double powp[7])
{
  powp[0] = 1.0;
  for(int n=1; n<=L; ++n) powp[n] = powp[n-1]*x;
#if defined(__SANITIZE_ADDRESS__) || defined(ADDRESS_SANITIZER)
  for(int n=L+1; n<=6; ++n){ powp[n]=0.0; }
#endif
}
inline void precompute_powers(double x,           // coordinate
      int    L,           // 0 ≤ L ≤ 6
      double powp[7],     // x⁰ … xᴸ   (output)
      double dpow[7])     // ∂/∂x x⁰ … xᴸ (output)
  {
    powp[0] = 1.0;      dpow[0] = 0.0;
    for(int n = 1; n <= L; ++n){
      powp[n] = powp[n-1] * x;           // xⁿ
      dpow[n] = n * powp[n-1];           // n·xⁿ⁻¹
    }
#if defined(__SANITIZE_ADDRESS__) || defined(ADDRESS_SANITIZER)
    for(int n=L+1; n<=6; ++n){ powp[n]=0.0; dpow[n]=0.0; }
#endif
  }

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

  s=sgrid.size()*(Lmax+1);
  rhoisize = gen_atomic_orbitals(Lmax);
  rhoisize *= sgrid.size();

  auto init_atoms = get_init_atoms(context);
  init_for_atoms(init_atoms);
}
template <typename F>
  void DM_mEAD<F>::calc_aed(core::rho_type& rhoi, core::aed_type& aed)
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
template <typename F>
  void DM_mEAD<F>::calc_dXijdri_dXjidri(
      int Zi, int Zj, double rij, double,
      const core::Vec3d& v, core::rho_type& rhoi, core::rho_type& rhoj,
      core::fd_type& fd_ij, double scale)
  {
    if(!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
    /* -------- constants & scratch arrays --------------------------- */
    const double x=v[0], y=v[1], z=v[2];
    const double wi = weights[Zi], wj = weights[Zj];

    const double fc  = fcut->calc       (rij);
    const double fcp = fcut->calc_prime (rij);

    /* radial part ---------------------------------------------------- */
    const std::size_t ng = cgrid.size();
    StackOrHeapBuffer<MAX_GRID> buf(ng);
    precompute_radials(rij, fc, fcp, sgrid, cgrid, buf.g(), buf.dg(), ng);


    /* Cartesian powers up to 6 -------------------------------------- */
    double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
    precompute_powers(x, Lmax ,px ,dpx);
    precompute_powers(y, Lmax ,py ,dpy);
    precompute_powers(z, Lmax ,pz ,dpz);

    /* negative versions are sign flips (avoid extra mul) */
    auto neg = [](double v,int n){ return (n&1)? -v : v; };

    /* --------------------------------------------------------------- */
    std::size_t Lshift = 0, ii = rhoisize + rfidx;

    for(int L=0; L<=Lmax; ++L){
      const auto& orbs = orbitals[L];

      for(std::size_t o=0; o<orbs.size(); ++o){
        const int lx=orbs[o][0], ly=orbs[o][1], lz=orbs[o][2];

        const double pxij =  px[lx], pyij =  py[ly], pzij =  pz[lz];
        const double pxji = neg(pxij,lx), pyji = neg(pyij,ly), pzji = neg(pzij,lz);

        const double pqrij = pxij*pyij*pzij;
        const double pqrji = pxji*pyji*pzji;

        const double txij = dpx[lx]*pyij*pzij;
        const double tyij = dpy[ly]*pxij*pzij;
        const double tzij = dpz[lz]*pyij*pxij;

        const double txji = neg(txij,lx-1);
        const double tyji = neg(tyij,ly-1);
        const double tzji = neg(tzij,lz-1);

        for(std::size_t c = 0; c < ng; ++c, ++ii)
        {
          const double bf  = buf.g()[c];           //   G(rij,…)
          const double dbf = buf.dg()[c];           //  dG/dr

          const double rijInv = 1.0 / rij;

          /* --- term-2, radial derivative -------------------------------- */
          const double term2x = dbf * rijInv * x;
          const double term2y = dbf * rijInv * y;
          const double term2z = dbf * rijInv * z;

          /* --- term-1, angular derivative ------------------------------- */
          const double T1ijx = bf * txij;
          const double T1ijy = bf * tyij;
          const double T1ijz = bf * tzij;

          const double T1jix = bf * txji;
          const double T1jiy = bf * tyji;
          const double T1jiz = bf * tzji;

          const double f_ri = rhoi(ii);
          const double f_rj = rhoj(ii);

          const std::size_t row = c + Lshift + fidx;

          /* ---- accumulate ---------------------------------------------- */
          fd_ij(row,0) += scale * ( f_ri*(T1ijx + pqrij*term2x)*wj
              - f_rj*(T1jix + pqrji*term2x)*wi );

          fd_ij(row,1) += scale * ( f_ri*(T1ijy + pqrij*term2y)*wj
              - f_rj*(T1jiy + pqrji*term2y)*wi );

          fd_ij(row,2) += scale * ( f_ri*(T1ijz + pqrij*term2z)*wj
              - f_rj*(T1jiz + pqrji*term2z)*wi );
        }
      }
      Lshift += ng;
    }
  }
template <typename F>
  void DM_mEAD<F>::calc_dXijdri(
      int Zi,int Zj,double rij,double,
      const core::Vec3d& v, core::rho_type& rhoi,
      core::fd_type& fd_ij, double scale)
  {
    if(!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

    const double x=v[0], y=v[1], z=v[2];
    const double wj = weights[Zj];

    const double fc  = scale*wj*fcut->calc       (rij);
    const double fcp = scale*wj*fcut->calc_prime (rij);

    const std::size_t ng = cgrid.size();
    StackOrHeapBuffer<MAX_GRID> buf(ng);
    precompute_radials(rij, fc, fcp, sgrid, cgrid, buf.g(), buf.dg(), ng);

    double px[7], py[7], pz[7], dpx[7], dpy[7], dpz[7];
    precompute_powers(x, Lmax ,px ,dpx);
    precompute_powers(y, Lmax ,py ,dpy);
    precompute_powers(z, Lmax ,pz ,dpz);

    std::size_t Lshift=0, ii=rhoisize+rfidx;
    const double rijInv = 1.0/rij;

    for(int L=0; L<=Lmax; ++L){
      const auto& orbs = orbitals[L];

      for(std::size_t o=0; o<orbs.size(); ++o){
        const int lx=orbs[o][0], ly=orbs[o][1], lz=orbs[o][2];

        const double pxij = px[lx],  pyij = py[ly],  pzij = pz[lz];
        const double pqrij = pxij*pyij*pzij;

        const double txij = dpx[lx]*pyij*pzij;
        const double tyij = dpy[ly]*pxij*pzij;
        const double tzij = dpz[lz]*pyij*pxij;

        for(std::size_t c=0; c<ng; ++c, ++ii){
          const double bf  = buf.g() [c];
          const double dbf = buf.dg()[c];

          const double term2x = pqrij*dbf*x*rijInv;
          const double term2y = pqrij*dbf*y*rijInv;
          const double term2z = pqrij*dbf*z*rijInv;

          const double f_r = rhoi(ii);
          const std::size_t row = c + Lshift + fidx;

          fd_ij(row,0) += f_r*( bf*txij + term2x );
          fd_ij(row,1) += f_r*( bf*tyij + term2y );
          fd_ij(row,2) += f_r*( bf*tzij + term2z );
        }
      }
      Lshift += ng;
    }
  }
template <typename F>
  std::string DM_mEAD<F>::label() {
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
      core::rho_type& rhoi,double scale)
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
    orbitals.resize(Lmax+1);
    fac.resize(Lmax+1);

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

} // namespace models
} // namespace tadah
