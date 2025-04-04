#include <tadah/models/descriptors/dm/dm_func.h>

#include <tadah/core/config.h>

F_RLR::F_RLR(Config &conf, size_t cSize, size_t sSize):
  config(&conf)
{
  if (!config->get<bool>("INITMB")) return;
  D_Base::get_grid(*config,"SEMBFUNC",sfuncgrid);
  D_Base::get_grid(*config,"CEMBFUNC",cfuncgrid);
  // check if size is correct, throw otherwise
  if (sfuncgrid.size() != sSize) {
    throw std::runtime_error("F_RLR: size of SEMBFUNC does not match the number of [C/S]GRIDMB.");
  }
  if (cfuncgrid.size() != cSize) {
    throw std::runtime_error("F_RLR: size of CEMBFUNC does not match the number of [C/S]GRIDMB.");
  }
}
// F_RLR::F_RLR():
//   config(nullptr)
// {}
double F_RLR::calc_F(double rho,size_t p) {
  return sfuncgrid[p]*rho*log(cfuncgrid[p]*rho);
}
double F_RLR::calc_dF(double rho,size_t p) {
  return sfuncgrid[p]*(log(cfuncgrid[p]*rho)+1);
}

//////////////////////////
// F_SQRT
F_SQRT::F_SQRT(Config &conf, size_t, size_t sSize):
  config(&conf)
{
  if (!config->get<bool>("INITMB")) return;
  if (config->exist("SEMBFUNC")) {
    D_Base::get_grid(*config, "SEMBFUNC", sfuncgrid);
    // check if size is correct, throw otherwise
    if (sfuncgrid.size() != sSize) {
      throw std::runtime_error("F_SQRT: size of SEMBFUNCT does not match the number of [C/S]GRIDMB.");
    }
  } else {
    // size_t n = config->size("SGRIDMB");
    sfuncgrid.resize(sSize);
    for (size_t i=0; i<sSize; ++i) {
      sfuncgrid[i] = 1.0;
    }
  }
}
// F_SQRT::F_SQRT():
//   config(nullptr)
// {
//   sfuncgrid.resize(1);
//   sfuncgrid[0] = 1.0;
// }
double F_SQRT::calc_F(double rho,size_t p) {
  if (rho>0)
    return sfuncgrid[p]*std::sqrt(rho);
  return 0;
}
double F_SQRT::calc_dF(double rho,size_t p) {
  if (rho>0)
    return 0.5*sfuncgrid[p]/(std::sqrt(rho));
  return 0;
}
