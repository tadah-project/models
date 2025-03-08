#include <tadah/models/descriptors/d3/d3_base.h>
#include <tadah/models/descriptors/d3/d3_dummy.h>
namespace tadah {
namespace models {


template<> CONFIG::Registry<D3_Base>::Map CONFIG::Registry<D3_Base>::registry{};
template<> CONFIG::Registry<D3_Base,tadah::core::Config&>::Map CONFIG::Registry<D3_Base,tadah::core::Config&>::registry{};

CONFIG::Registry<D3_Base>::Register<D3_Dummy> D3_Dummy_0( "D3_Dummy" );
CONFIG::Registry<D3_Base,tadah::core::Config&>::Register<D3_Dummy> D3_Dummy_1( "D3_Dummy" );

} // namespace models
} // namespace tadah

