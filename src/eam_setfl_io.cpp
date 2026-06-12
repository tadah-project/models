#include <tadah/models/eam/eam_setfl_io.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace tadah {
namespace models {

void EAMSetflIO::read_setfl(const std::string &path) {
  file_path_ = path;
  std::ifstream in(path);
  if (!in.good())
    throw std::runtime_error("EAMSetflIO: setfl file does not exist: " + path);

  std::string line;
  std::getline(in, line); // comment 1
  std::getline(in, line); // comment 2
  std::getline(in, line); // comment 3
  std::string line2 = line;  // (setfl: comment; funcfl: "nrho drho nr dr rcut")
  std::getline(in, line); // line 4: <ntypes> <El1> [<El2> ...]
  {
    std::istringstream iss(line);
    int ntypes = 0;
    if (!(iss >> ntypes) || ntypes <= 0) {
      // DYNAMO funcfl heuristic: its line 3 is "nrho drho nr dr rcut"
      // (int double int double double) where setfl has a comment.
      std::istringstream f3(line2);
      int fn1 = 0, fn2 = 0; double fd1 = 0, fd2 = 0, fd3 = 0;
      if (f3 >> fn1 >> fd1 >> fn2 >> fd2 >> fd3)
        throw std::runtime_error(
            "EAMSetflIO: '" + path + "' looks like a DYNAMO funcfl file "
            "(pair_style eam). Only setfl files (eam/alloy or eam/fs) are "
            "supported in this version; convert the potential to setfl or "
            "obtain the setfl distribution.");
      throw std::runtime_error(
          "EAMSetflIO: malformed setfl header in '" + path +
          "' (line 4 must start with the number of element types).");
    }
    if (ntypes != 1)
      throw std::runtime_error(
          "EAMSetflIO: only single-element setfl (eam/alloy or eam/fs) is "
          "supported in this version; '" + path + "' declares " +
          std::to_string(ntypes) + " element types. Multi-element refit is a "
          "planned extension.");
    iss >> element_;
  }

  in >> nrho_ >> drho_ >> nr_ >> dr_ >> rcut_;
  in >> atomic_number_ >> atomic_mass_ >> lattice_param_ >> lattice_;
  if (!in.good())
    throw std::runtime_error("EAMSetflIO: malformed setfl header in: " + path);

  rdr_   = 1.0 / dr_;
  rdrho_ = 1.0 / drho_;

  frho_.resize(nrho_);
  rhor_.resize(nr_);
  z2r_.resize(nr_);
  for (int i = 0; i < nrho_; ++i) in >> frho_[i];
  for (int i = 0; i < nr_; ++i)   in >> rhor_[i];
  for (int i = 0; i < nr_; ++i)   in >> z2r_[i];
  if (in.fail())
    throw std::runtime_error("EAMSetflIO: truncated tabulated data in: " + path);

  // LAMMPS convention: the embedding spline is valid up to (nrho-1)*drho.
  rhomax_ = (nrho_ - 1) * drho_;
}

// Natural cubic spline, identical to the former D2_EAM/DM_EAM gen_splines.
void EAMSetflIO::gen_spline(int n, double delta,
                            const std::vector<double> &f,
                            std::vector<std::vector<double>> &spline) {
  // in lammps f is n+1, here is size n
  for (int m = 1; m <= n; m++) spline[m][6] = f[m - 1];

  spline[1][5]   = spline[2][6] - spline[1][6];
  spline[2][5]   = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5]   = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
    spline[m][5] = ((spline[m - 2][6] - spline[m + 2][6]) +
                    8.0 * (spline[m + 1][6] - spline[m - 1][6])) / 12.0;

  for (int m = 1; m <= n - 1; m++) {
    spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) - 2.0 * spline[m][5] - spline[m + 1][5];
    spline[m][3] = spline[m][5] + spline[m + 1][5] - 2.0 * (spline[m + 1][6] - spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5] / delta;
    spline[m][1] = 2.0 * spline[m][4] / delta;
    spline[m][0] = 3.0 * spline[m][3] / delta;
  }
}

void EAMSetflIO::build_splines() {
  frho_spline_.assign(nrho_ + 1, std::vector<double>(7, 0.0));
  rhor_spline_.assign(nr_ + 1, std::vector<double>(7, 0.0));
  z2r_spline_.assign(nr_ + 1, std::vector<double>(7, 0.0));
  gen_spline(nrho_, drho_, frho_, frho_spline_);
  gen_spline(nr_,   dr_,   rhor_, rhor_spline_);
  gen_spline(nr_,   dr_,   z2r_,  z2r_spline_);
}

double EAMSetflIO::eval_value(const std::vector<std::vector<double>> &spline,
                              double x, double rdx, int n) {
  double p = x * rdx + 1.0;
  int m = static_cast<int>(p);
  if (m < 1) m = 1;
  if (m > n - 1) m = n - 1;
  p -= m;
  if (p > 1.0) p = 1.0;
  const std::vector<double> &c = spline[m];
  return ((c[3] * p + c[4]) * p + c[5]) * p + c[6];
}

double EAMSetflIO::phi_tab(double r) const {
  return eval_value(z2r_spline_, r, rdr_, nr_) / r;
}

double EAMSetflIO::z2_tab(double r) const {
  return eval_value(z2r_spline_, r, rdr_, nr_);
}

double EAMSetflIO::rho_tab(double r) const {
  return eval_value(rhor_spline_, r, rdr_, nr_);
}

double EAMSetflIO::F_tab(double rho) const {
  return eval_value(frho_spline_, rho, rdrho_, nrho_);
}

EAMGrid EAMSetflIO::scaled_grid(double scale) const {
  EAMGrid g;
  g.nr   = std::max(2, static_cast<int>(std::lround(nr_ * scale)));
  g.nrho = std::max(2, static_cast<int>(std::lround(nrho_ * scale)));
  // Preserve the physical spans (nr-1)*dr and (nrho-1)*drho.
  g.dr   = dr_   * static_cast<double>(nr_ - 1)   / static_cast<double>(g.nr - 1);
  g.drho = drho_ * static_cast<double>(nrho_ - 1) / static_cast<double>(g.nrho - 1);
  g.rcut = rcut_;
  return g;
}

void EAMSetflIO::write_setfl(const std::string &path,
                             const EAMGrid &grid,
                             const std::function<double(double)> &frho_fn,
                             const std::function<double(double)> &rhor_fn,
                             const std::function<double(double)> &z2r_fn,
                             const std::string &title) const {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("EAMSetflIO: cannot open output setfl file: " + path);

  // 3 comment lines (LAMMPS skips exactly three) + element line.
  out << title << "\n";
  out << "Native tadah refit functions (NOT splines of the original table)\n";
  out << "source: " << (file_path_.empty() ? "-" : file_path_) << "\n";
  out << " 1 " << (element_.empty() ? "X" : element_) << "\n";

  out << std::scientific << std::setprecision(16);
  out << grid.nrho << "  " << grid.drho << "  " << grid.nr << "  " << grid.dr
      << "  " << grid.rcut << "\n";
  out << atomic_number_ << "  " << atomic_mass_ << "  " << lattice_param_ << "  "
      << (lattice_.empty() ? "unknown" : lattice_) << "\n";

  auto dump = [&out](int N, double spacing,
                     const std::function<double(double)> &fn) {
    for (int i = 0; i < N; ++i) {
      out << " " << fn(i * spacing);
      if ((i + 1) % 5 == 0) out << "\n";
    }
    if (N % 5 != 0) out << "\n";
  };

  dump(grid.nrho, grid.drho, frho_fn);  // F(rho)
  dump(grid.nr,   grid.dr,   rhor_fn);  // rho(r)
  dump(grid.nr,   grid.dr,   z2r_fn);   // z2r(r) = phi(r)*r
}

} // namespace models
} // namespace tadah
