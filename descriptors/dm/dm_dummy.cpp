#include "dm_dummy.h"

Registry<DM_Base,Config&>::Register<DM_Dummy> DM_Dummy_1( "DM_Dummy" );

DM_Dummy::DM_Dummy() {}
DM_Dummy::DM_Dummy(Config&) {}
void DM_Dummy::calc_aed(
        rho_rtype,
        aed_rtype ) {}
int DM_Dummy::calc_dXijdri_dXjidri(
        const double,
        const double,
        const Eigen::Vector3d &,
        const double,
        const double,
        rho_rtype,
        rho_rtype,
        fd_type &,
        const double,
        const double) { return 0; }
int DM_Dummy::calc_dXijdri(
        const double,
        const double,
        const Eigen::Vector3d &,
        const double,
        const double,
        rho_rtype,
        fd_type &) { return 0; }
size_t DM_Dummy::size() { return s;}
std::string DM_Dummy::label() {return lab;}
void DM_Dummy::init_rhoi(rho_rtype) {}
void DM_Dummy::calc_rho(
        const double,
        const double,
        const double,
        const Eigen::Vector3d &,
        rho_rtype) {}
size_t DM_Dummy::rhoi_size() { return 0; }
size_t DM_Dummy::rhoip_size() { return 0; }

