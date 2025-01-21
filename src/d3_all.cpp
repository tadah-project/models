#include <tadah/models/descriptors/d3/d3_base.h>
#include <tadah/models/descriptors/d3/d3_dummy.h>

template<> CONFIG::Registry<D3_Base>::Map CONFIG::Registry<D3_Base>::registry{};
template<> CONFIG::Registry<D3_Base,Config&>::Map CONFIG::Registry<D3_Base,Config&>::registry{};

CONFIG::Registry<D3_Base>::Register<D3_Dummy> D3_Dummy_0( "D3_Dummy" );
CONFIG::Registry<D3_Base,Config&>::Register<D3_Dummy> D3_Dummy_1( "D3_Dummy" );
