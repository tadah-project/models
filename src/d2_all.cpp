#include <tadah/models/descriptors/d2/d2_base.h>
#include <tadah/models/descriptors/d2/d2_dummy.h>
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/models/descriptors/d2/d2_bp.h>
#include <tadah/models/descriptors/d2/d2_eam.h>
#include <tadah/models/descriptors/d2/d2_lj.h>
#include <tadah/models/descriptors/d2/d2_mie.h>
#include <tadah/models/descriptors/d2/d2_join.h>
#include <tadah/models/descriptors/d2/d2_join_n.h>

template<> CONFIG::Registry<D2_Base,Config&>::Map CONFIG::Registry<D2_Base,Config&>::registry{};

CONFIG::Registry<D2_Base,Config&>::Register<D2_Blip> D2_Blip_1( "D2_Blip" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_BP> D2_BP_1( "D2_BP" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_Dummy> D2_Dummy_1( "D2_Dummy" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_EAM> D2_EAM_1( "D2_EAM" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_LJ> D2_LJ_1( "D2_LJ" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_MIE> D2_MIE_1( "D2_MIE" );



// This is in development
CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_LJ,D2_BP>> D2_JOIN_1( "D2_LJ+BP" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_LJ,D2_Blip>> D2_JOIN_2( "D2_LJ+Blip" );

CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_MIE,D2_BP>> D2_JOIN_3( "D2_MIE+BP" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_MIE,D2_Blip>> D2_JOIN_4( "D2_MIE+Blip" );

CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_Blip,D2_Blip>> D2_JOIN_5( "D2_Blip_x2" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_BP,D2_BP>> D2_JOIN_6( "D2_BP_x2" );

CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_EAM,D2_BP>> D2_JOIN_7( "D2_EAM+BP" );
CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_EAM,D2_Blip>> D2_JOIN_8( "D2_EAM+Blip" );

CONFIG::Registry<D2_Base,Config&>::Register<D2_Join<D2_BP,D2_Blip>> D2_JOIN_9( "D2_BP+Blip" );

// CONFIG::Registry<D2_Base,Config&>::Register<D2_Join_N<D2_LJ,D2_LJ>> D2_JOIN_3( "D2_LJ_D2_LJ" );
//CONFIG::Registry<D2_Base,Config&>::Register<D2_Join_N<D2_LJ,D2_LJ,D2_LJ>> D2_JOIN_4( "D2_LJ_D2_LJ_D2_LJ" );
