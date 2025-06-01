#include <tadah/models/descriptors/dm/dm_all.h>
namespace tadah {
namespace models {


template<> tadah::core::Registry<DM_Base>::Map tadah::core::Registry<DM_Base>::registry{};
template<> tadah::core::Registry<DM_Base,tadah::core::Context&>::Map tadah::core::Registry<DM_Base,tadah::core::Context&>::registry{};
template<> tadah::core::Registry<DM_Base,const tadah::core::Context&>::Map tadah::core::Registry<DM_Base,const tadah::core::Context&>::registry{};

tadah::core::Registry<D_Base>::Register<DM_Blip> DM_Blip_0( "DM_Blip" );
tadah::core::Registry<D_Base>::Register<DM_Dummy> DM_Dummy_0( "DM_Dummy" );
tadah::core::Registry<D_Base>::Register<DM_EAD> DM_EAD_0( "DM_EAD" );
tadah::core::Registry<D_Base>::Register<DM_EAM> DM_EAM_0( "DM_EAM" );
tadah::core::Registry<D_Base>::Register<DM_mEAD<F_RLR>> DM_mEAD_0("DM_mRLR");
tadah::core::Registry<D_Base>::Register<DM_mEAD<F_SQRT>> DM_mEAD_00("DM_mSQRT");
tadah::core::Registry<D_Base>::Register<DM_mJoin> DM_mJoin_0("DM_mJoin");

tadah::core::Registry<DM_Base>::Register<DM_Blip> DM_Blip_1( "DM_Blip" );
tadah::core::Registry<DM_Base>::Register<DM_Dummy> DM_Dummy_1( "DM_Dummy" );
tadah::core::Registry<DM_Base>::Register<DM_EAD> DM_EAD_1( "DM_EAD" );
tadah::core::Registry<DM_Base>::Register<DM_EAM> DM_EAM_1( "DM_EAM" );
tadah::core::Registry<DM_Base>::Register<DM_mEAD<F_RLR>> DM_mEAD_1("DM_mRLR");
tadah::core::Registry<DM_Base>::Register<DM_mEAD<F_SQRT>> DM_mEAD_11("DM_mSQRT");
tadah::core::Registry<DM_Base>::Register<DM_mJoin> DM_mJoin_1("DM_mJoin");

tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Blip> DM_Blip_2( "DM_Blip" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Dummy> DM_Dummy_2( "DM_Dummy" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_EAD> DM_EAD_2( "DM_EAD" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_EAM> DM_EAM_2( "DM_EAM" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_RLR>> DM_mEAD_2("DM_mRLR");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_SQRT>> DM_mEAD_22("DM_mSQRT");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mJoin> DM_mJoin_2( "DM_mJoin" );
//template<> Registry<DM_Base,F_Base&,tadah::core::Context&>::Map Registry<DM_Base,F_Base&,tadah::core::Context&>::registry{};
//Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_RLR>> DM_mEAD_1( "DM_mEAD" );

} // namespace models
} // namespace tadah
