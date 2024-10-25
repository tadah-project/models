#include <tadah/models/descriptors/d2/d2_base.h>
#include <tadah/models/descriptors/d2/d2_dummy.h>
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/models/descriptors/d2/d2_bp.h>
#include <tadah/models/descriptors/d2/d2_eam.h>
#include <tadah/models/descriptors/d2/d2_lj.h>
#include <tadah/models/descriptors/d2/d2_mie.h>

template<> CONFIG::Registry<D2_Base,Config&>::Map CONFIG::Registry<D2_Base,Config&>::registry{};

CONFIG::Registry<D2_Base,Config&>::Register<D2_Blip> D2_Blip_1( "D2_Blip" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_BP> D2_BP_1( "D2_BP" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_Dummy> D2_Dummy_1( "D2_Dummy" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_EAM> D2_EAM_1( "D2_EAM" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_LJ> D2_LJ_1( "D2_LJ" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_MIE> D2_MIE_1( "D2_MIE" );

CONFIG::Registry<DM_Base, Config&>::Register<DM_mEAD<F_RLR>> DM_mEAM_1("DM_mEAD");