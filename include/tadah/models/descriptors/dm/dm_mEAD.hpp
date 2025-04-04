#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_mEAD.h>
#include <tadah/models/descriptors/d_basis_functions.h>

template <typename F>
DM_mEAD<F>::DM_mEAD() {
  init();
}

/* Index ii iterates rhoi array
 * rhoi array is split 50/50 into density and derivative
 * of the embedding function
 */
template <typename F>
DM_mEAD<F>::DM_mEAD(Config &config): DM_Base(config),
  rc(config.get<double>("RCUTMBMAX")),
  emb(config)
{

  if (!config.get<bool>("INITMB")) return;
  init();

  get_grid(config,"CGRIDMB",cgrid);
  get_grid(config,"SGRIDMB",sgrid);
  // update config with generated grid
  config.remove("CGRIDMB");
  for (const auto &i: cgrid) config.add("CGRIDMB",i);
  config.remove("SGRIDMB");
  for (const auto &i: sgrid) config.add("SGRIDMB",i);

  if (cgrid.size()!=sgrid.size()) {
    throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
  }

  Lmax = config.get<int>("TYPEMB",1);
  s=sgrid.size()*(Lmax+1);
  rhoisize = gen_atomic_orbitals(Lmax);
  rhoisize *= sgrid.size();

  if (verbose) {
    std::cout << std::endl;
    std::cout << "SGRID (SGRIDMB): " << sgrid.size() << std::endl;
    for (auto e:sgrid) std::cout << e << "  ";
    std::cout << std::endl;

    std::cout << "CGRID (CGRIDMB): " << cgrid.size() << std::endl;
    for (auto L:cgrid) std::cout << L << "  ";
    std::cout << std::endl;
    std::cout << "rhoisize: " << rhoisize << std::endl;
  }
  auto init_atoms = get_init_atoms(config);
  init_for_atoms(init_atoms);
}

template <typename F>
void DM_mEAD<F>::calc_aed(
  rho_type& rhoi,
  aed_type &aed)
{
  size_t Lshift = 0;
  size_t ii=rfidx;
  for (int L=0; L<=Lmax; ++L) {
    for (size_t o=0; o<orbitals[L].size(); ++o) {
      const double f = fac[L][o];
      for (size_t c=0; c<cgrid.size(); c++) {
        aed(c+Lshift+fidx) += f*emb.calc_F(rhoi(ii),c);
        rhoi(ii+rhoisize) = f*emb.calc_dF(rhoi(ii),c);
        ii++;
      }
    }
    Lshift += sgrid.size();
  }
}
template <typename F>
void DM_mEAD<F>::calc_dXijdri_dXjidri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const Vec3d &vec_ij,
  rho_type& rhoi,
  rho_type& rhoj,
  fd_type &fd_ij,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  const double x = vec_ij[0];
  const double y = vec_ij[1];
  const double z = vec_ij[2];

  size_t Lshift = 0;

  double wi = weights[Zi];
  double wj = weights[Zj];
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  size_t ii=rhoisize+rfidx;
  for (int L=0; L<=Lmax; ++L) {
    for (size_t o=0; o<orbitals[L].size(); ++o) {
      const int lx = orbitals[L][o][0];
      const int ly = orbitals[L][o][1];
      const int lz = orbitals[L][o][2];

      const double powxlxij = my_pow(x,lx);
      const double powylyij = my_pow(y,ly);
      const double powzlzij = my_pow(z,lz);

      const double powxlxji = my_pow(-x,lx);
      const double powylyji = my_pow(-y,ly);
      const double powzlzji = my_pow(-z,lz);

      double txij=0.0,tyij=0.0,tzij=0.0;
      double txji=0.0,tyji=0.0,tzji=0.0;
      if (lx!=0) {
        txij = lx*my_pow(x,lx-1)*powylyij*powzlzij;
        txji = lx*my_pow(-x,lx-1)*powylyji*powzlzji;
      }
      if (ly!=0) {
        tyij = ly*my_pow(y,ly-1)*powxlxij*powzlzij;
        tyji = ly*my_pow(-y,ly-1)*powxlxji*powzlzji;
      }
      if (lz!=0) {
        tzij = lz*my_pow(z,lz-1)*powylyij*powxlxij;
        tzji = lz*my_pow(-z,lz-1)*powylyji*powxlxji;
      }
      const double pqrij = powxlxij*powylyij*powzlzij;
      const double pqrji = powxlxji*powylyji*powzlzji;

      for (size_t c=0; c<cgrid.size(); c++) {
        double bf= G(rij,sgrid[c],cgrid[c],fc_ij);
        double dbf= dG(rij,sgrid[c],cgrid[c],fc_ij,fcp_ij);

        const double term2ijx=pqrij*dbf*x/rij;
        const double term2ijy=pqrij*dbf*y/rij;
        const double term2ijz=pqrij*dbf*z/rij;

        const double term2jix=-pqrji*dbf*x/rij;
        const double term2jiy=-pqrji*dbf*y/rij;
        const double term2jiz=-pqrji*dbf*z/rij;

        const double term1ijx = bf*txij;
        const double term1ijy = bf*tyij;
        const double term1ijz = bf*tzij;

        const double term1jix = bf*txji;
        const double term1jiy = bf*tyji;
        const double term1jiz = bf*tzji;

        fd_ij(c+Lshift+fidx,0) += scale*(
          rhoi(ii)*(term1ijx + term2ijx)*wj
          -rhoj(ii)*(term1jix + term2jix)*wi)
        ;

        fd_ij(c+Lshift+fidx,1) += scale*(
          rhoi(ii)*(term1ijy + term2ijy)*wj
          -rhoj(ii)*(term1jiy + term2jiy)*wi)
        ;

        fd_ij(c+Lshift+fidx,2) += scale*(
          rhoi(ii)*(term1ijz + term2ijz)*wj
          -rhoj(ii)*(term1jiz + term2jiz)*wi)
        ;
        ii++;
      }
    }
    Lshift+=sgrid.size();
  }
}
template <typename F>
void DM_mEAD<F>::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const Vec3d &vec_ij,
  rho_type& rhoi,
  fd_type &fd_ij,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  const double x = vec_ij[0];
  const double y = vec_ij[1];
  const double z = vec_ij[2];

  double wj = weights[Zj];
  double fc_ij = scale*wj*fcut->calc(rij);
  double fcp_ij = scale*wj*fcut->calc_prime(rij);

  size_t Lshift = 0;

  size_t ii=rhoisize+rfidx;
  for (int L=0; L<=Lmax; ++L) {
    for (size_t o=0; o<orbitals[L].size(); ++o) {
      const int lx = orbitals[L][o][0];
      const int ly = orbitals[L][o][1];
      const int lz = orbitals[L][o][2];

      const double powxlxij = my_pow(x,lx);
      const double powylyij = my_pow(y,ly);
      const double powzlzij = my_pow(z,lz);

      double txij=0.0,tyij=0.0,tzij=0.0;
      if (lx!=0) {
        txij = lx*my_pow(x,lx-1)*powylyij*powzlzij;
      }
      if (ly!=0) {
        tyij = ly*my_pow(y,ly-1)*powxlxij*powzlzij;
      }
      if (lz!=0) {
        tzij = lz*my_pow(z,lz-1)*powylyij*powxlxij;
      }
      const double pqrij = powxlxij*powylyij*powzlzij;

      for (size_t c=0; c<cgrid.size(); c++) {
        double bf= G(rij,sgrid[c],cgrid[c],fc_ij);
        double dbf= dG(rij,sgrid[c],cgrid[c],fc_ij,fcp_ij);

        const double term2ijx=pqrij*dbf*x/rij;
        const double term2ijy=pqrij*dbf*y/rij;
        const double term2ijz=pqrij*dbf*z/rij;

        const double term1ijx = bf*txij;
        const double term1ijy = bf*tyij;
        const double term1ijz = bf*tzij;

        fd_ij(c+Lshift+fidx,0) +=
          rhoi(ii)*(term1ijx + term2ijx)
        ;

        fd_ij(c+Lshift+fidx,1) +=
          rhoi(ii)*(term1ijy + term2ijy)
        ;

        fd_ij(c+Lshift+fidx,2) += 
          rhoi(ii)*(term1ijz + term2ijz)
        ;
        ii++;
      }
    }
    Lshift+=sgrid.size();
  }
}
template <typename F>
std::string DM_mEAD<F>::label() {
  return lab;
}

template <typename F>
void DM_mEAD<F>::init_rhoi(rho_type& rhoi)
{
  rhoi.resize(2*rhoisize);
  rhoi.set_zero();
}
template <typename F>
void DM_mEAD<F>::calc_rho(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const Vec3d &vec_ij,
  rho_type& rhoi,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  size_t ii=rfidx;
  double fc_ij = scale*weights[Zj]*fcut->calc(rij);
  for (int L=0; L<=Lmax; ++L) {
    for (size_t o=0; o<orbitals[L].size(); ++o) {
      const size_t lx = orbitals[L][o][0];
      const size_t ly = orbitals[L][o][1];
      const size_t lz = orbitals[L][o][2];
      double t = (my_pow(vec_ij[0],lx)*my_pow(vec_ij[1],ly)*my_pow(vec_ij[2],lz));
      for (size_t c=0; c<cgrid.size(); c++) {
        rhoi(ii++) += G(rij,sgrid[c],cgrid[c],fc_ij)*t;
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
      for (int ly = 0; ly <= L; ++ly) {
        for (int lz = 0; lz <= L; ++lz) {
          if (lx+ly+lz == L) {
            std::vector<int> o = {lx, ly, lz};
            orbitals[L].push_back(o);
            const double f = fact(L)/(fact(lx)*fact(ly)*fact(lz))/my_pow(4*rc,L);
            fac[L].push_back(f);
            count++;
          }
        }
      }
    }
  }
  return count;
}
template <typename F>
double DM_mEAD<F>::fact(int n)
{
  double f=1.0;
  while (n) {
    f *= n;
    --n;
  }
  return f;
}
template <typename F>
double DM_mEAD<F>::my_pow(double x, int n){
  double r = 1.0;

  while(n > 0){
    r *= x;
    --n;
  }
  return r;
}
template <typename F>
void DM_mEAD<F>::init() {
  keys.push_back("CGRIDMB");
  keys.push_back("SGRIDMB");
  // keys.push_back("CEMBFUNC");
  // keys.push_back("SEMBFUNC");
  nparams=3;
}
