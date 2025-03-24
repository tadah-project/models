#include <tadah/models/descriptors/dm/dm_func.h>

#include <tadah/core/config.h>

F_RLR::F_RLR(Config &conf):
  config(&conf)
{
  if (!config->get<bool>("INITMB")) return;
  D_Base::get_grid(*config,"SEMBFUNC",sgrid);
  D_Base::get_grid(*config,"CEMBFUNC",cgrid);

}
F_RLR::F_RLR():
  config(nullptr)
{}
double F_RLR::calc_F(double rho,size_t p) {
  if (rho>0)
    return sgrid[p]*rho*log(cgrid[p]*rho);
  return 0;
}
double F_RLR::calc_dF(double rho,size_t p) {
  if (rho>0)
    return sgrid[p]*(log(cgrid[p]*rho)+1);
  return sgrid[p];
}

//////////////////////////
// F_SQRT
F_SQRT::F_SQRT(Config &conf):
  config(&conf)
{
  if (!config->get<bool>("INITMB"))
    return;
  if (config->exist("SEMBFUNC")) {
    D_Base::get_grid(*config, "SEMBFUNC", sgrid);
  } else {
    sgrid.resize(1);
    sgrid[0] = 1.0;
  }
}
F_SQRT::F_SQRT():
  config(nullptr)
{
  sgrid.resize(1);
  sgrid[0] = 1.0;
}
double F_SQRT::calc_F(double rho,size_t p) {
  if (rho>0)
    return sgrid[p]*std::sqrt(rho);
  return 0;
}
double F_SQRT::calc_dF(double rho,size_t p) {
  if (rho>0)
    return 0.5*sgrid[p]/(std::sqrt(rho));
  return 0;
}
