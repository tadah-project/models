#include "catch2/catch.hpp"
#include <limits>
#include <stdexcept>
#include "../cutoffs/cutoffs.h"
#include <string>

TEST_CASE( "Testing Cutoffs: Cut_Dummy", "[Cut_Dummy]" ) {

    double rcut2b = 6.2;
    double rcut2bsq = rcut2b*rcut2b;
    using Cut = Cut_Dummy;
    std::string cuttype = "Cut_Dummy";
    Cut_Base *c2b = new Cut( rcut2b );

    REQUIRE( c2b->label() == cuttype );

    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    // cutoff cannot be negative
    double temp = -0.1;
    REQUIRE_THROWS(Cut( temp ));
    REQUIRE_THROWS_AS(c2b->set_rcut(temp), std::runtime_error);
    REQUIRE_THROWS_AS(c2b->test_rcut(temp), std::runtime_error);

    // recheck after resetting cutoff
    rcut2b=3.4;
    rcut2bsq=rcut2b*rcut2b;
    c2b->set_rcut(100000);
    c2b->set_rcut(rcut2b);
    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    delete c2b;
}
TEST_CASE( "Testing Cutoffs: Cut_Cos", "[Cut_Cos]" ) {

    double rcut2b = 6.2;
    double rcut2bsq = rcut2b*rcut2b;
    using Cut = Cut_Cos;
    std::string cuttype = "Cut_Cos";
    Cut_Base *c2b = new Cut( rcut2b );

    REQUIRE( c2b->label() == cuttype );

    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    // cutoff cannot be negative
    double temp = -0.1;
    REQUIRE_THROWS(Cut( temp ));
    REQUIRE_THROWS_AS(c2b->set_rcut(temp), std::runtime_error);
    REQUIRE_THROWS_AS(c2b->test_rcut(temp), std::runtime_error);

    // recheck after resetting cutoff
    rcut2b=3.4;
    rcut2bsq=rcut2b*rcut2b;
    c2b->set_rcut(100000);
    c2b->set_rcut(rcut2b);
    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    delete c2b;
}
TEST_CASE( "Testing Cutoffs: Cut_Tanh", "[Cut_Tanh]" ) {

    double rcut2b = 6.2;
    double rcut2bsq = rcut2b*rcut2b;
    using Cut = Cut_Tanh;
    std::string cuttype = "Cut_Tanh";
    Cut_Base *c2b = new Cut( rcut2b );

    REQUIRE( c2b->label() == cuttype );

    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    // cutoff cannot be negative
    double temp = -0.1;
    REQUIRE_THROWS(Cut( temp ));
    REQUIRE_THROWS_AS(c2b->set_rcut(temp), std::runtime_error);
    REQUIRE_THROWS_AS(c2b->test_rcut(temp), std::runtime_error);

    // recheck after resetting cutoff
    rcut2b=3.4;
    rcut2bsq=rcut2b*rcut2b;
    c2b->set_rcut(100000);
    c2b->set_rcut(rcut2b);
    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    delete c2b;
}
TEST_CASE( "Testing Cutoffs: Cut_Poly2", "[Cut_Poly2]" ) {

    double rcut2b = 6.2;
    double rcut2bsq = rcut2b*rcut2b;
    using Cut = Cut_Poly2;
    std::string cuttype = "Cut_Poly2";
    Cut_Base *c2b = new Cut( rcut2b );

    REQUIRE( c2b->label() == cuttype );

    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    // cutoff cannot be negative
    double temp = -0.1;
    REQUIRE_THROWS(Cut( temp ));
    REQUIRE_THROWS_AS(c2b->set_rcut(temp), std::runtime_error);
    REQUIRE_THROWS_AS(c2b->test_rcut(temp), std::runtime_error);

    // recheck after resetting cutoff
    rcut2b=3.4;
    rcut2bsq=rcut2b*rcut2b;
    c2b->set_rcut(100000);
    c2b->set_rcut(rcut2b);
    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

    delete c2b;
}
