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
 * rhoi array is split into density and derivative
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
    throw std::runtime_error("SGRIDMB and CGRIDMB arrays differ in size.\n");
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

  odim = emb.output_dim();
  const std::size_t ng = cgrid.size();
  const std::size_t norbs = gen_atomic_orbitals(Lmax);

  s = odim * (Lmax+1);
  rhoisize  = norbs * ng;      // ρ storage
  rhoipsize = norbs * odim;    // dF/dρ storage

  set_rfidx(0);
  set_drfidx(rhoisize);

  auto init_atoms = get_init_atoms(context);
  init_for_atoms(init_atoms);
}

/* ================================================================== *
 *  calc_aed – compute atomic energy descriptor from accumulated ρ    *
 * ================================================================== */
template <typename F>
  void DM_mEAD<F>::calc_aed(core::rho_type& rhoi, core::aed_type& aed) const
  {
    const std::size_t ng = cgrid.size();

    if constexpr (F::expands_dimensions) {
      /* --- dimension-expanding path (e.g. F_Blip, ng==1) ------------ */
      std::size_t Lshift = 0;           // offset into aed / fd_ij rows
      std::size_t ii_rho = rfidx;       // index into ρ   (one per orbital)
      std::size_t ii_drho = drfidx;     // index into dF/dρ storage

      for (int L = 0; L <= Lmax; ++L) {
        for (std::size_t o = 0; o < orbitals[L].size(); ++o) {
          const double f = fac[L][o];
          const double rho = rhoi(ii_rho);     // single ρ for this orbital
          ++ii_rho;

          for (std::size_t p = 0; p < odim; ++p, ++ii_drho) {
            const double FF  = emb.calc_F (rho, p);
            const double dFF = emb.calc_dF(rho, p);

            aed(p + Lshift + fidx) += f * FF;
            rhoi(ii_drho)           = f * dFF;
          }
        }
        Lshift += odim;
      }
    } else {
      /* --- standard path (odim == ng) ------------------------------- */
      std::size_t Lshift = 0;
      std::size_t ii     = rfidx;
      std::size_t dii    = drfidx;

      for (int L = 0; L <= Lmax; ++L) {
        for (std::size_t o = 0; o < orbitals[L].size(); ++o) {
          const double f = fac[L][o];

          for (std::size_t c = 0; c < ng; ++c, ++ii, ++dii) {
            const double rho = rhoi(ii);
            const double FF  = emb.calc_F (rho, c);
            const double dFF = emb.calc_dF(rho, c);

            aed(c + Lshift + fidx) += f * FF;
            rhoi(dii) = f * dFF;
          }
        }
        Lshift += ng;
      }
    }
  }

/* ------------------------------------------------------------------ *
 *  dX_ij/dr_i  and  dX_ji/dr_i   (pair contribution, symmetric)     *
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

  if constexpr (F::expands_dimensions) {
    /* --- dimension-expanding path --------------------------------- */
    std::size_t Lshift = 0;
    std::size_t ii_drho = drfidx;  // dF/dρ storage index

    for (int L = 0; L <= Lmax; ++L)
    {
      const double s_L = (L & 1) ?  1.0 : -1.0;

      for (std::size_t o = 0; o < orbitals[L].size(); ++o)
      {
        const int lx = orbitals[L][o][0],
                  ly = orbitals[L][o][1],
                  lz = orbitals[L][o][2];

        const double pxij = px[lx], pyij = py[ly], pzij = pz[lz];
        const double pqrij = pxij*pyij*pzij;

        const double txij = dpx[lx]*pyij*pzij;
        const double tyij = dpy[ly]*pxij*pzij;
        const double tzij = dpz[lz]*pyij*pxij;

        /* single radial basis (ng==1) */
        const double g  = buf.g ()[0];
        const double dg = buf.dg()[0];

        const double Bx = g*txij + pqrij*dg*x*rinv;
        const double By = g*tyij + pqrij*dg*y*rinv;
        const double Bz = g*tzij + pqrij*dg*z*rinv;

        /* loop over embedding output channels */
        for (std::size_t p = 0; p < odim; ++p, ++ii_drho)
        {
          const double coef = scale *
                              (  rhoi(ii_drho)*wj
                              - s_L * rhoj(ii_drho)*wi );

          const std::size_t row = p + Lshift + fidx;
          fd_ij(row,0) += coef * Bx;
          fd_ij(row,1) += coef * By;
          fd_ij(row,2) += coef * Bz;
        }
      }
      Lshift += odim;
    }
  } else {
    /* --- standard path -------------------------------------------- */
    std::size_t Lshift = 0;
    std::size_t ii     = drfidx;

    for (int L = 0; L <= Lmax; ++L)
    {
      const double s_L = (L & 1) ?  1.0 : -1.0;

      for (std::size_t o = 0; o < orbitals[L].size(); ++o)
      {
        const int lx = orbitals[L][o][0],
                  ly = orbitals[L][o][1],
                  lz = orbitals[L][o][2];

        const double pxij = px[lx], pyij = py[ly], pzij = pz[lz];
        const double pqrij = pxij*pyij*pzij;

        const double txij = dpx[lx]*pyij*pzij;
        const double tyij = dpy[ly]*pxij*pzij;
        const double tzij = dpz[lz]*pyij*pxij;

        for (std::size_t c = 0; c < ng; ++c, ++ii)
        {
          const double g  = buf.g ()[c];
          const double dg = buf.dg()[c];

          const double Bx = g*txij + pqrij*dg*x*rinv;
          const double By = g*tyij + pqrij*dg*y*rinv;
          const double Bz = g*tzij + pqrij*dg*z*rinv;

          const double coef = scale *
                              (  rhoi(ii)*wj
                              - s_L * rhoj(ii)*wi );

          const std::size_t row = c + Lshift + fidx;
          fd_ij(row,0) += coef * Bx;
          fd_ij(row,1) += coef * By;
          fd_ij(row,2) += coef * Bz;
        }
      }
      Lshift += ng;
    }
  }
}

/* ================================================================== *
 *  one–sided descriptor derivative     dX_ij / dr_i                  *
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

  if constexpr (F::expands_dimensions) {
    /* ---- dimension-expanding path -------------------------------- */
    std::size_t Lshift = 0;
    std::size_t ii_drho = drfidx;

    for (int L = 0; L <= Lmax; ++L)
    {
      const auto& orbs = orbitals[L];

      for (std::size_t o = 0; o < orbs.size(); ++o)
      {
        const int lx = orbs[o][0],
                  ly = orbs[o][1],
                  lz = orbs[o][2];

        const double pxij = px[lx], pyij = py[ly], pzij = pz[lz];
        const double pqrij = pxij*pyij*pzij;

        const double txij = dpx[lx]*pyij*pzij;
        const double tyij = dpy[ly]*pxij*pzij;
        const double tzij = dpz[lz]*pyij*pxij;

        /* single radial basis (ng==1) */
        const double g  = buf.g ()[0];
        const double dg = buf.dg()[0];

        const double termx = g*txij + dg*pqrij*rinv*x;
        const double termy = g*tyij + dg*pqrij*rinv*y;
        const double termz = g*tzij + dg*pqrij*rinv*z;

        /* loop over embedding output channels */
        for (std::size_t p = 0; p < odim; ++p, ++ii_drho)
        {
          const double f_r = rhoi(ii_drho);

          const std::size_t row = p + Lshift + fidx;
          fd_ij(row,0) += f_r * termx;
          fd_ij(row,1) += f_r * termy;
          fd_ij(row,2) += f_r * termz;
        }
      }
      Lshift += odim;
    }
  } else {
    /* ---- standard path ------------------------------------------- */
    std::size_t Lshift = 0;
    std::size_t ii     = drfidx;

    for (int L = 0; L <= Lmax; ++L)
    {
      const auto& orbs = orbitals[L];

      for (std::size_t o = 0; o < orbs.size(); ++o)
      {
        const int lx = orbs[o][0],
                  ly = orbs[o][1],
                  lz = orbs[o][2];

        const double pxij = px[lx],  pyij = py[ly],  pzij = pz[lz];

        const double pqrij = pxij*pyij*pzij;

        const double txij = dpx[lx]*pyij*pzij;
        const double tyij = dpy[ly]*pxij*pzij;
        const double tzij = dpz[lz]*pyij*pxij;

        const double t2x_pref = pqrij * rinv * x;
        const double t2y_pref = pqrij * rinv * y;
        const double t2z_pref = pqrij * rinv * z;

#pragma omp simd
        for (std::size_t c = 0; c < ng; ++c)
        {
          const double g  = buf.g ()[c];
          const double dg = buf.dg()[c];

          const double termx = g*txij + dg*t2x_pref;
          const double termy = g*tyij + dg*t2y_pref;
          const double termz = g*tzij + dg*t2z_pref;

          const double f_r = rhoi(ii + c);

          const std::size_t row = c + Lshift + fidx;
          fd_ij(row,0) += f_r * termx;
          fd_ij(row,1) += f_r * termy;
          fd_ij(row,2) += f_r * termz;
        }
        ii += ng;
      }
      Lshift += ng;
    }
  }
}


template <typename F>
  std::string DM_mEAD<F>::label() const {
    return lab;
  }

template <typename F>
  void DM_mEAD<F>::init_rhoi(tadah::core::rho_type& rhoi)
  {
    rhoi.resize(rhoisize+rhoipsize);
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
    odim     = emb.output_dim();
    const std::size_t ng = cgrid.size();
    const std::size_t norbs = gen_atomic_orbitals(Lmax);
    s         = odim * static_cast<size_t>(Lmax + 1);
    rhoisize  = norbs * ng;
    rhoipsize = norbs * odim;
    set_drfidx(rhoisize);
    return s != old_s;
  }

/* ================================================================== *
 *  meta() — default templated body: covers the parameter shape every
 *  DM_mEAD<F> instantiation shares (RCUTMB + CGRIDMB + SGRIDMB), plus
 *  the per-param ρ-collapse check on CGRIDMB. F-specific extras
 *  (CEMBFUNC/SEMBFUNC for F_RLR, RCUTENV/CEMBFUNC/SEMBFUNC for
 *  F_mEnv<Cut_*>) are added by explicit specialisations declared in
 *  dm_all.cpp; those specialisations call back into the templated
 *  body to share the common parameter declarations.
 * ================================================================== */
template <typename F>
std::optional<DescriptorMeta> DM_mEAD<F>::meta() const {
  using DataStats = tadah::core::stats::DataStats;
  static constexpr std::size_t I_P10 = tadah::core::stats::IDX_P10;
  static constexpr std::size_t I_P50 = tadah::core::stats::IDX_P50;
  static constexpr std::size_t I_P90 = tadah::core::stats::IDX_P90;
  static constexpr std::size_t I_P99 = tadah::core::stats::IDX_P99;

  DescriptorMeta dm;
  dm.label = label();   // runtime label string (e.g. "DM_mRLR")

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
    if (d.per_pair.empty() && d.n_atoms > 0) {
      return HealthFail{"RCUTMB=" + std::to_string(rc) +
                        " is too small — no MB neighbours found."};
    }
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
  p_cgrid.role     = "GTO centre r_s. Must lie within dataset r_ij range or ρ collapses.";
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
