#include <tadah/models/descriptors/dm/dm_all.h>

template<> CONFIG::Registry<DM_Base>::Map CONFIG::Registry<DM_Base>::registry{};
template<> CONFIG::Registry<DM_Base,Config&>::Map CONFIG::Registry<DM_Base,Config&>::registry{};

CONFIG::Registry<D_Base>::Register<DM_Blip> DM_Blip_0( "DM_Blip" );
CONFIG::Registry<D_Base>::Register<DM_Dummy> DM_Dummy_0( "DM_Dummy" );
CONFIG::Registry<D_Base>::Register<DM_EAD> DM_EAD_0( "DM_EAD" );
CONFIG::Registry<D_Base>::Register<DM_EAM> DM_EAM_0( "DM_EAM" );
CONFIG::Registry<D_Base>::Register<DM_mEAD<F_RLR>> DM_mEAD_0("DM_mRLR");
CONFIG::Registry<D_Base>::Register<DM_mEAD<F_SQRT>> DM_mEAD_00("DM_mSQRT");
CONFIG::Registry<D_Base>::Register<DM_mJoin> DM_mJoin_0("DM_mJoin");

CONFIG::Registry<DM_Base>::Register<DM_Blip> DM_Blip_1( "DM_Blip" );
CONFIG::Registry<DM_Base>::Register<DM_Dummy> DM_Dummy_1( "DM_Dummy" );
CONFIG::Registry<DM_Base>::Register<DM_EAD> DM_EAD_1( "DM_EAD" );
CONFIG::Registry<DM_Base>::Register<DM_EAM> DM_EAM_1( "DM_EAM" );
CONFIG::Registry<DM_Base>::Register<DM_mEAD<F_RLR>> DM_mEAD_1("DM_mRLR");
CONFIG::Registry<DM_Base>::Register<DM_mEAD<F_SQRT>> DM_mEAD_11("DM_mSQRT");
CONFIG::Registry<DM_Base>::Register<DM_mJoin> DM_mJoin_1("DM_mJoin");

CONFIG::Registry<DM_Base,Config&>::Register<DM_Blip> DM_Blip_2( "DM_Blip" );
CONFIG::Registry<DM_Base,Config&>::Register<DM_Dummy> DM_Dummy_2( "DM_Dummy" );
CONFIG::Registry<DM_Base,Config&>::Register<DM_EAD> DM_EAD_2( "DM_EAD" );
CONFIG::Registry<DM_Base,Config&>::Register<DM_EAM> DM_EAM_2( "DM_EAM" );
CONFIG::Registry<DM_Base,Config&>::Register<DM_mEAD<F_RLR>> DM_mEAD_2("DM_mRLR");
CONFIG::Registry<DM_Base,Config&>::Register<DM_mEAD<F_SQRT>> DM_mEAD_22("DM_mSQRT");
CONFIG::Registry<DM_Base,Config&>::Register<DM_mJoin> DM_mJoin_2( "DM_mJoin" );
//template<> Registry<DM_Base,F_Base&,Config&>::Map Registry<DM_Base,F_Base&,Config&>::registry{};
//Registry<DM_Base,Config&>::Register<DM_mEAD<F_RLR>> DM_mEAD_1( "DM_mEAD" );
