#include <tadah/models/descriptors/dm/dm_func.h>

namespace tadah {
namespace models {


F_RLR::F_RLR(tadah::core::Context &c):
  context(&c)
{
  if (!context->get<bool>("INITMB")) return;
  D_Base::get_grid(*context,"SEMBFUNC",sgrid);
  D_Base::get_grid(*context,"CEMBFUNC",cgrid);

}
F_RLR::F_RLR():
  context(nullptr)
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

