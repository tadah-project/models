#include <tadah/models/functions/f_all.h>
#include <tadah/core/registry.h>
namespace tadah {
namespace models {


template<> CONFIG::Registry<tadah::core::Function_Base>::Map CONFIG::Registry<tadah::core::Function_Base>::registry{};
template<> CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Map CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::registry{};

CONFIG::Registry<tadah::core::Function_Base>::Register<BF_Linear> BF_Linear_1( "BF_Linear" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<BF_Linear> BF_Linear_2( "BF_Linear" );
CONFIG::Registry<tadah::core::Function_Base>::Register<BF_Polynomial2> BF_Polynomial2_1( "BF_Polynomial2" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<BF_Polynomial2> BF_Polynomial2_2( "BF_Polynomial2" );
CONFIG::Registry<tadah::core::Function_Base>::Register<Kern_Linear> Kern_Linear_1( "Kern_Linear" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<Kern_Linear> Kern_Linear_2( "Kern_Linear" );
CONFIG::Registry<tadah::core::Function_Base>::Register<Kern_LQ> Kern_LQ_1( "Kern_LQ" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<Kern_LQ> Kern_LQ_2( "Kern_LQ" );
CONFIG::Registry<tadah::core::Function_Base>::Register<Kern_Polynomial> Kern_Polynomial_1( "Kern_Polynomial" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<Kern_Polynomial> Kern_Polynomial_2( "Kern_Polynomial" );
CONFIG::Registry<tadah::core::Function_Base>::Register<Kern_Quadratic> Kern_Quadratic_1( "Kern_Quadratic" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<Kern_Quadratic> Kern_Quadratic_2( "Kern_Quadratic" );
CONFIG::Registry<tadah::core::Function_Base>::Register<Kern_RBF> Kern_RBF_1( "Kern_RBF" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<Kern_RBF> Kern_RBF_2( "Kern_RBF" );
CONFIG::Registry<tadah::core::Function_Base>::Register<Kern_Sigmoid> Kern_Sigmoid_1( "Kern_Sigmoid" );
CONFIG::Registry<tadah::core::Function_Base,tadah::core::Config&>::Register<Kern_Sigmoid> Kern_Sigmoid_2( "Kern_Sigmoid" );


} // namespace models
} // namespace tadah

