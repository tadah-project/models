#include "dm_base.h"
#include "dm_dummy.h"
#include "dm_eam.h"
#include "dm_ead.h"
#include "dm_blip.h"
#include "dm_mEAD.h"

template<> CONFIG::Registry<DM_Base,Config&>::Map CONFIG::Registry<DM_Base,Config&>::registry{};

CONFIG::Registry<DM_Base,Config&>::Register<DM_Blip> DM_Blip_1( "DM_Blip" );
CONFIG::Registry<DM_Base,Config&>::Register<DM_Dummy> DM_Dummy_1( "DM_Dummy" );
CONFIG::Registry<DM_Base,Config&>::Register<DM_EAD> DM_EAD_1( "DM_EAD" );
CONFIG::Registry<DM_Base,Config&>::Register<DM_EAM> DM_EAM_1( "DM_EAM" );

//template<> Registry<DM_Base,F_Base&,Config&>::Map Registry<DM_Base,F_Base&,Config&>::registry{};
//Registry<DM_Base,Config&>::Register<DM_mEAD<F_RLR>> DM_mEAD_1( "DM_mEAD" );
