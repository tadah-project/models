#include "d3_base.h"
#include "d3_dummy.h"

template<> CONFIG::Registry<D3_Base,Config&>::Map CONFIG::Registry<D3_Base,Config&>::registry{};

CONFIG::Registry<D3_Base,Config&>::Register<D3_Dummy> D3_Dummy_1( "D3_Dummy" );
