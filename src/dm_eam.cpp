#include "tadah/models/descriptors/dm/dm_base.h"
#include <tadah/models/descriptors/dm/dm_eam.h>
namespace tadah {
namespace models {


DM_EAM::DM_EAM() {
  init();
}
DM_EAM::DM_EAM(tadah::core::Config &c): DM_Base(c)
{
  init();
  if (!c.get<bool>("INITMB")) return;
  s=1;
  ef.file_path = c.get<std::string>("TYPEMB",1);
  read_setfl();

  if (std::abs(ef.rcut - c.get<double>("RCUTMB")) > std::numeric_limits<double>::min()) {
    if (verbose) {
      std::cout << "tadah::core::Config file cutoff and setfl file cutoff differ: "
        << c.get<double>("RCUTMB") << " " << ef.rcut << std::endl;
      std::cout << "Enforcing SETFL file cutoff: " << ef.rcut << std::endl;
    }
    c.remove("RCUTMB");
    c.add("RCUTMB", ef.rcut);
    c.postprocess();
  }

  frho_spline.resize(ef.nrho+1, std::vector<double>(7));
  rhor_spline.resize(ef.nr+1, std::vector<double>(7));
  z2r_spline.resize(ef.nr+1, std::vector<double>(7));

  gen_splines(ef.nrho, ef.drho, ef.frho, frho_spline);
  gen_splines(ef.nr, ef.dr, ef.rhor, rhor_spline);
  gen_splines(ef.nr, ef.dr, ef.z2r, z2r_spline);

  rhoisize = 1;
  auto init_atoms = get_init_atoms(c);
  init_for_atoms(init_atoms);
}

void DM_EAM::calc_aed(
  tadah::core::rho_type& rho,
  tadah::core::aed_type &aed)
{
  //std::cout << "rho[i]:   " << rho(0) << std::endl;
  double p = rho(rfidx)*ef.rdrho + 1.0;
  int m = static_cast<int> (p);
  m = std::max(1,std::min(m,ef.nrho-1));
  p -= m;
  p = std::min(p,1.0);
  //std::vector<double> coeff = frho_spline[m];
  double phi = ((frho_spline[m][3]*p + frho_spline[m][4])*p + frho_spline[m][5])*p + frho_spline[m][6];

  rho(1+rfidx) = (frho_spline[m][0]*p + frho_spline[m][1])*p + frho_spline[m][2];    // lammps fp[i]
  //std::cout << "fp[i]:   " << rho(1) << std::endl;

  if (rho(rfidx) > ef.rhomax) phi += rho(1+rfidx) * (rho(rfidx)-ef.rhomax);
  //phi *= scale[type[i]][type[i]];

  aed(fidx) += phi;


}
void DM_EAM::calc_dXijdri_dXjidri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rhoi,
  tadah::core::rho_type& rhoj,
  tadah::core::fd_type &fd_ij,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double r = rij;
  double rinv = 1/r;
  double p = r*ef.rdr + 1.0;
  int m = static_cast<int> (p);
  m = std::min(m,ef.nr-1);
  p -= m;
  p = std::min(p,1.0);
  double rhoip = (rhor_spline[m][0]*p + rhor_spline[m][1])*p + rhor_spline[m][2];
  double rhojp = (rhor_spline[m][0]*p + rhor_spline[m][1])*p + rhor_spline[m][2];

  // here rho_i[1] is a derivative of rho_i at local atom, similarly for rho_j
  //double psip = rho_i[1]*rhojp + rho_j[1]*rhoip;
  double psip = rhoi(rfidx+1)*rhojp + rhoj(rfidx+1)*rhoip;

  double f = scale*psip*rinv;
  fd_ij(fidx,0) = f*vec_ij[0]; 
  fd_ij(fidx,1) = f*vec_ij[1]; 
  fd_ij(fidx,2) = f*vec_ij[2]; 
}
void DM_EAM::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rhoi,
  tadah::core::fd_type &fd_ij,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double r = rij;
  double rinv = 1/r;
  double p = r*ef.rdr + 1.0;
  int m = static_cast<int> (p);
  m = std::min(m,ef.nr-1);
  p -= m;
  p = std::min(p,1.0);
  double rhojp = (rhor_spline[m][0]*p + rhor_spline[m][1])*p + rhor_spline[m][2];

  // here rho_i[1] is a derivative of rho_i at local atom
  double psip = rhoi(rfidx+1)*rhojp;

  double f = scale*psip*rinv;
  fd_ij(fidx,0) = f*vec_ij[0];
  fd_ij(fidx,1) = f*vec_ij[1];
  fd_ij(fidx,2) = f*vec_ij[2];
}
std::string DM_EAM::label() {
  return lab;
}

void DM_EAM::read_setfl()
{
  std::string line;
  std::ifstream in_file(ef.file_path);
  if (!in_file.good())
    throw std::runtime_error("SETFL file does not exists.\n");

  if (in_file.is_open()) {
    if (verbose)
      std::cout << "<DM_EAM> Reading setfl: " << ef.file_path << std::endl;
    // skip ficgridt 3 comment lines
    getline (in_file,line);
    getline (in_file,line);
    getline (in_file,line);
    // skip number of types and types
    getline (in_file,line);
    // read 5th line
    in_file >> ef.nrho >> ef.drho >> ef.nr >> ef.dr >> ef.rcut;
    in_file >> ef.atomic_number >> ef.atomic_mass >> ef.lattice_param >> ef.lattice;
    ef.rdr = 1.0/ef.dr;
    ef.rdrho = 1.0/ef.drho;
    // prepare arrays
    ef.frho.resize(ef.nrho);
    ef.rhor.resize(ef.nr);
    ef.z2r.resize(ef.nr);
    // read all data
    for (int i=0; i<ef.nrho; ++i) in_file >> ef.frho[i];
    for (int i=0; i<ef.nr; ++i) in_file >> ef.rhor[i];
    for (int i=0; i<ef.nr; ++i) in_file >> ef.z2r[i];
    in_file.close();
  }
  else {
    std::cout << "<DM_EAM> Unable to open file: " << ef.file_path << std::endl;;
  }
  // get max value of rho
  //ef.rhomax = *std::max_element(ef.rhor.begin(), ef.rhor.end());
  ef.rhomax = (ef.nrho-1) * ef.drho;

}
void DM_EAM::gen_splines(int &n, double &delta, std::vector<double> &f,
                         std::vector<std::vector<double>> &spline)
{
  // in lammps f is n+1, here is size n
  for (int m=1; m<=n; m++) spline[m][6] = f[m-1];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
      8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) - 2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] - 2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}
void DM_EAM::init_rhoi(tadah::core::rho_type& rhoi)
{
  rhoi.resize(2);
  rhoi.set_zero();
}
void DM_EAM::calc_rho(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  const tadah::core::Vec3d &,
  tadah::core::rho_type& rho,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double r = rij;
  double p = r*ef.rdr + 1.0;
  int m = static_cast<int> (p);
  m = std::min(m,ef.nr-1);
  p -= m;
  p = std::min(p,1.0);
  std::vector<double> coeff = rhor_spline[m];
  rho(rfidx) += scale*( ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);

}
void DM_EAM::init() {
  nparams=1;
}

} // namespace models
} // namespace tadah

