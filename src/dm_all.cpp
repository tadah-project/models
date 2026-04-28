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
tadah::core::Registry<D_Base>::Register<DM_mEAD<F_Blip>> DM_mEAD_000("DM_mBlip");
tadah::core::Registry<D_Base>::Register<DM_mEAD<F_mEnv<Cut_SinInv>>> DM_mEAD_0000("DM_mBlipEnv");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQ>>                    DM_PS_SQ_0   ("DM_PS_SQ");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQRT>>                  DM_PS_SQRT_0 ("DM_PS_SQRT");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_RLR>>                   DM_PS_RLR_0  ("DM_PS_RLR");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_Blip>>                  DM_PS_Blip_0 ("DM_PS_Blip");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQ,   GaussianRadial>>           DM_GPS_SQ_0  ("DM_GPS_SQ");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQRT, GaussianRadial>>           DM_GPS_SQRT_0("DM_GPS_SQRT");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_RLR,  GaussianRadial>>           DM_GPS_RLR_0 ("DM_GPS_RLR");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_Blip, GaussianRadial>>           DM_GPS_Blip_0("DM_GPS_Blip");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQ,   BlipRadial, OrbitalWise>>  DM_OW_SQ_0   ("DM_OW_SQ");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQRT, BlipRadial, OrbitalWise>>  DM_OW_SQRT_0 ("DM_OW_SQRT");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_RLR,  BlipRadial, OrbitalWise>>  DM_OW_RLR_0  ("DM_OW_RLR");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_Blip, BlipRadial, OrbitalWise>>  DM_OW_Blip_0 ("DM_OW_Blip");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQ,   GaussianRadial,OrbitalWise>> DM_GOW_SQ_0  ("DM_GOW_SQ");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_SQRT, GaussianRadial,OrbitalWise>> DM_GOW_SQRT_0("DM_GOW_SQRT");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_RLR,  GaussianRadial,OrbitalWise>> DM_GOW_RLR_0 ("DM_GOW_RLR");
tadah::core::Registry<D_Base>::Register<DM_Poly<F_Blip, GaussianRadial,OrbitalWise>> DM_GOW_Blip_0("DM_GOW_Blip");
tadah::core::Registry<D_Base>::Register<DM_mJoin> DM_mJoin_0("DM_mJoin");

tadah::core::Registry<DM_Base>::Register<DM_Blip> DM_Blip_1( "DM_Blip" );
tadah::core::Registry<DM_Base>::Register<DM_Dummy> DM_Dummy_1( "DM_Dummy" );
tadah::core::Registry<DM_Base>::Register<DM_EAD> DM_EAD_1( "DM_EAD" );
tadah::core::Registry<DM_Base>::Register<DM_EAM> DM_EAM_1( "DM_EAM" );
tadah::core::Registry<DM_Base>::Register<DM_mEAD<F_RLR>> DM_mEAD_1("DM_mRLR");
tadah::core::Registry<DM_Base>::Register<DM_mEAD<F_SQRT>> DM_mEAD_11("DM_mSQRT");
tadah::core::Registry<DM_Base>::Register<DM_mEAD<F_Blip>> DM_mEAD_111("DM_mBlip");
tadah::core::Registry<DM_Base>::Register<DM_mEAD<F_mEnv<Cut_SinInv>>> DM_mEAD_1111("DM_mBlipEnv");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQ>>                    DM_PS_SQ_1   ("DM_PS_SQ");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQRT>>                  DM_PS_SQRT_1 ("DM_PS_SQRT");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_RLR>>                   DM_PS_RLR_1  ("DM_PS_RLR");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_Blip>>                  DM_PS_Blip_1 ("DM_PS_Blip");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQ,   GaussianRadial>>           DM_GPS_SQ_1  ("DM_GPS_SQ");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQRT, GaussianRadial>>           DM_GPS_SQRT_1("DM_GPS_SQRT");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_RLR,  GaussianRadial>>           DM_GPS_RLR_1 ("DM_GPS_RLR");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_Blip, GaussianRadial>>           DM_GPS_Blip_1("DM_GPS_Blip");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQ,   BlipRadial, OrbitalWise>>  DM_OW_SQ_1   ("DM_OW_SQ");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQRT, BlipRadial, OrbitalWise>>  DM_OW_SQRT_1 ("DM_OW_SQRT");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_RLR,  BlipRadial, OrbitalWise>>  DM_OW_RLR_1  ("DM_OW_RLR");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_Blip, BlipRadial, OrbitalWise>>  DM_OW_Blip_1 ("DM_OW_Blip");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQ,   GaussianRadial,OrbitalWise>> DM_GOW_SQ_1  ("DM_GOW_SQ");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_SQRT, GaussianRadial,OrbitalWise>> DM_GOW_SQRT_1("DM_GOW_SQRT");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_RLR,  GaussianRadial,OrbitalWise>> DM_GOW_RLR_1 ("DM_GOW_RLR");
tadah::core::Registry<DM_Base>::Register<DM_Poly<F_Blip, GaussianRadial,OrbitalWise>> DM_GOW_Blip_1("DM_GOW_Blip");
tadah::core::Registry<DM_Base>::Register<DM_mJoin> DM_mJoin_1("DM_mJoin");

tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Blip> DM_Blip_2( "DM_Blip" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Dummy> DM_Dummy_2( "DM_Dummy" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_EAD> DM_EAD_2( "DM_EAD" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_EAM> DM_EAM_2( "DM_EAM" );
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_RLR>> DM_mEAD_2("DM_mRLR");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_SQRT>> DM_mEAD_22("DM_mSQRT");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_Blip>> DM_mEAD_222("DM_mBlip");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_mEnv<Cut_SinInv>>> DM_mEAD_2222("DM_mBlipEnv");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQ>>                    DM_PS_SQ_2   ("DM_PS_SQ");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQRT>>                  DM_PS_SQRT_2 ("DM_PS_SQRT");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_RLR>>                   DM_PS_RLR_2  ("DM_PS_RLR");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_Blip>>                  DM_PS_Blip_2 ("DM_PS_Blip");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQ,   GaussianRadial>>           DM_GPS_SQ_2  ("DM_GPS_SQ");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQRT, GaussianRadial>>           DM_GPS_SQRT_2("DM_GPS_SQRT");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_RLR,  GaussianRadial>>           DM_GPS_RLR_2 ("DM_GPS_RLR");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_Blip, GaussianRadial>>           DM_GPS_Blip_2("DM_GPS_Blip");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQ,   BlipRadial, OrbitalWise>>  DM_OW_SQ_2   ("DM_OW_SQ");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQRT, BlipRadial, OrbitalWise>>  DM_OW_SQRT_2 ("DM_OW_SQRT");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_RLR,  BlipRadial, OrbitalWise>>  DM_OW_RLR_2  ("DM_OW_RLR");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_Blip, BlipRadial, OrbitalWise>>  DM_OW_Blip_2 ("DM_OW_Blip");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQ,   GaussianRadial,OrbitalWise>> DM_GOW_SQ_2  ("DM_GOW_SQ");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_SQRT, GaussianRadial,OrbitalWise>> DM_GOW_SQRT_2("DM_GOW_SQRT");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_RLR,  GaussianRadial,OrbitalWise>> DM_GOW_RLR_2 ("DM_GOW_RLR");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_Poly<F_Blip, GaussianRadial,OrbitalWise>> DM_GOW_Blip_2("DM_GOW_Blip");
tadah::core::Registry<DM_Base,tadah::core::Context&>::Register<DM_mJoin> DM_mJoin_2( "DM_mJoin" );

tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Blip> cDM_Blip_2( "DM_Blip" );
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Dummy> cDM_Dummy_2( "DM_Dummy" );
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_EAD> cDM_EAD_2( "DM_EAD" );
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_EAM> cDM_EAM_2( "DM_EAM" );
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_mEAD<F_RLR>> cDM_mEAD_2("DM_mRLR");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_mEAD<F_SQRT>> cDM_mEAD_22("DM_mSQRT");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_mEAD<F_Blip>> cDM_mEAD_222("DM_mBlip");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_mEAD<F_mEnv<Cut_SinInv>>> cDM_mEAD_2222("DM_mBlipEnv");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQ>>                    cDM_PS_SQ_2   ("DM_PS_SQ");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQRT>>                  cDM_PS_SQRT_2 ("DM_PS_SQRT");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_RLR>>                   cDM_PS_RLR_2  ("DM_PS_RLR");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_Blip>>                  cDM_PS_Blip_2 ("DM_PS_Blip");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQ,   GaussianRadial>>           cDM_GPS_SQ_2  ("DM_GPS_SQ");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQRT, GaussianRadial>>           cDM_GPS_SQRT_2("DM_GPS_SQRT");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_RLR,  GaussianRadial>>           cDM_GPS_RLR_2 ("DM_GPS_RLR");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_Blip, GaussianRadial>>           cDM_GPS_Blip_2("DM_GPS_Blip");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQ,   BlipRadial, OrbitalWise>>  cDM_OW_SQ_2   ("DM_OW_SQ");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQRT, BlipRadial, OrbitalWise>>  cDM_OW_SQRT_2 ("DM_OW_SQRT");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_RLR,  BlipRadial, OrbitalWise>>  cDM_OW_RLR_2  ("DM_OW_RLR");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_Blip, BlipRadial, OrbitalWise>>  cDM_OW_Blip_2 ("DM_OW_Blip");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQ,   GaussianRadial,OrbitalWise>> cDM_GOW_SQ_2  ("DM_GOW_SQ");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_SQRT, GaussianRadial,OrbitalWise>> cDM_GOW_SQRT_2("DM_GOW_SQRT");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_RLR,  GaussianRadial,OrbitalWise>> cDM_GOW_RLR_2 ("DM_GOW_RLR");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_Poly<F_Blip, GaussianRadial,OrbitalWise>> cDM_GOW_Blip_2("DM_GOW_Blip");
tadah::core::Registry<DM_Base,const tadah::core::Context&>::Register<DM_mJoin> cDM_mJoin_2( "DM_mJoin" );
//template<> Registry<DM_Base,F_Base&,tadah::core::Context&>::Map Registry<DM_Base,F_Base&,tadah::core::Context&>::registry{};
//Registry<DM_Base,tadah::core::Context&>::Register<DM_mEAD<F_RLR>> DM_mEAD_1( "DM_mEAD" );

} // namespace models
} // namespace tadah
