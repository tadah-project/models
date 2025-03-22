#include <tadah/models/descriptors/dm/dm_func.h>

#include <tadah/core/context.h>
namespace tadah {
namespace models {


F_RLR::F_RLR(tadah::core::Context &conf):
  config(&conf)
{
  if (!config->get<bool>("INITMB")) return;
  D_Base::get_grid(*config,"SEMBFUNC",sgrid);
  D_Base::get_grid(*config,"CEMBFUNC",cgrid);

}
F_RLR::F_RLR():
  config(nullptr)
{
}
double F_RLR::calc_F(double rho,size_t c) {
  if (rho>0)
    return sgrid[c]*rho*log(cgrid[c]*rho);
  return 0;
}
double F_RLR::calc_dF(double rho,size_t c) {
  if (rho>0)
    return sgrid[c]*(log(cgrid[c]*rho)+1);
  return sgrid[c];
}

} // namespace models
} // namespace tadah

