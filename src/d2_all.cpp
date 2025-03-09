#include "tadah/models/descriptors/d2/d2_all.h"
namespace tadah {
namespace models {

template<> tadah::core::Registry<D_Base>::Map tadah::core::Registry<D_Base>::registry{};
template<> tadah::core::Registry<D_Base,tadah::core::Config&>::Map tadah::core::Registry<D_Base,tadah::core::Config&>::registry{};
template<> tadah::core::Registry<D2_Base,tadah::core::Config&>::Map tadah::core::Registry<D2_Base,tadah::core::Config&>::registry{};

tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_Blip> D2_Blip_00( "D2_Blip" );
tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_BP> D2_BP_00( "D2_BP" );
tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_Dummy> D2_Dummy_00( "D2_Dummy" );
tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_EAM> D2_EAM_00( "D2_EAM" );
tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_LJ> D2_LJ_00( "D2_LJ" );
tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_MIE> D2_MIE_00( "D2_MIE" );
tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_ZBL> D2_ZBL_00( "D2_ZBL" );
tadah::core::Registry<D_Base,tadah::core::Config&>::Register<D2_mJoin> D2_mJoin_00( "D2_mJoin" );

tadah::core::Registry<D_Base>::Register<D2_Blip> D2_Blip_0( "D2_Blip" );
tadah::core::Registry<D_Base>::Register<D2_BP> D2_BP_0( "D2_BP" );
tadah::core::Registry<D_Base>::Register<D2_Dummy> D2_Dummy_0( "D2_Dummy" );
tadah::core::Registry<D_Base>::Register<D2_EAM> D2_EAM_0( "D2_EAM" );
tadah::core::Registry<D_Base>::Register<D2_LJ> D2_LJ_0( "D2_LJ" );
tadah::core::Registry<D_Base>::Register<D2_MIE> D2_MIE_0( "D2_MIE" );
tadah::core::Registry<D_Base>::Register<D2_ZBL> D2_ZBL_0( "D2_ZBL" );
tadah::core::Registry<D_Base>::Register<D2_mJoin> D2_mJoin_0( "D2_mJoin" );

tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Blip> D2_Blip_1( "D2_Blip" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_BP> D2_BP_1( "D2_BP" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Dummy> D2_Dummy_1( "D2_Dummy" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_EAM> D2_EAM_1( "D2_EAM" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_LJ> D2_LJ_1( "D2_LJ" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_MIE> D2_MIE_1( "D2_MIE" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_ZBL> D2_ZBL_1( "D2_ZBL" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_mJoin> D2_mJoin_1( "D2_mJoin" );



// This is in development
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_LJ,D2_BP>> D2_JOIN_1( "D2_LJ+BP" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_LJ,D2_Blip>> D2_JOIN_2( "D2_LJ+Blip" );

tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_MIE,D2_BP>> D2_JOIN_3( "D2_MIE+BP" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_MIE,D2_Blip>> D2_JOIN_4( "D2_MIE+Blip" );

tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_Blip,D2_Blip>> D2_JOIN_5( "D2_Blip_x2" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_BP,D2_BP>> D2_JOIN_6( "D2_BP_x2" );

tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_EAM,D2_BP>> D2_JOIN_7( "D2_EAM+BP" );
tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_EAM,D2_Blip>> D2_JOIN_8( "D2_EAM+Blip" );

tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join<D2_BP,D2_Blip>> D2_JOIN_9( "D2_BP+Blip" );

// tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join_N<D2_LJ,D2_LJ>> D2_JOIN_3( "D2_LJ_D2_LJ" );
//tadah::core::Registry<D2_Base,tadah::core::Config&>::Register<D2_Join_N<D2_LJ,D2_LJ,D2_LJ>> D2_JOIN_4( "D2_LJ_D2_LJ_D2_LJ" );

} // namespace models
} // namespace tadah

