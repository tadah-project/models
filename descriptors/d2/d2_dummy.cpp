#include "d2_dummy.h"

Registry<D2_Base,Config&>::Register<D2_Dummy> D2_Dummy_1( "D2_Dummy" );

D2_Dummy::D2_Dummy() {}
D2_Dummy::D2_Dummy(Config &c):
    verbose(c.get<int>("VERBOSE"))
{}

void D2_Dummy::calc_aed(
        const double,
        const double,
        const double,
        aed_rtype ) {}
void D2_Dummy::calc_dXijdri(
        const double,
        const double,
        const double,
        const double,
        fd_type &) {}
void D2_Dummy::calc_all(
        const double,
        const double,
        const double,
        const double,
        aed_rtype ,
        fd_type &) {}
size_t D2_Dummy::size() { return s; }
std::string D2_Dummy::label() { return lab; }
