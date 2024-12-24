#include "catch2/catch.hpp"
#include <limits>
#include <stdexcept>
#include <tadah/models/cutoffs.h>
#include <string>

TEST_CASE( "Testing Cutoffs: Cut_Dummy", "[Cut_Dummy]" ) {

    REQUIRE_NOTHROW(Cut_Dummy());

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

    REQUIRE_NOTHROW(Cut_Cos());
    double rcut2b = 6.2;
    double rcut2bsq = rcut2b*rcut2b;
    using Cut = Cut_Cos;
    std::string cuttype = "Cut_Cos";
    Cut_Base *c2b = new Cut( rcut2b );

    REQUIRE( c2b->label() == cuttype );

    REQUIRE( c2b->calc(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( c2b->calc(0) == 1);
    REQUIRE( c2b->calc_prime(rcut2b) < std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut()-rcut2b)<std::numeric_limits<double>::min() );
    REQUIRE( std::abs(c2b->get_rcut_sq()-rcut2bsq)<std::numeric_limits<double>::min() );

  REQUIRE_THAT(c2b->calc(2.0),
               Catch::Matchers::WithinAbs(0.76448200516348, 1e-12));
  REQUIRE_THAT(c2b->calc_prime(2.0),
               Catch::Matchers::WithinAbs(-0.21500762619812, 1e-12));


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
    REQUIRE_NOTHROW(Cut_Tanh());

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
  REQUIRE_THAT(c2b->calc(2.0),
               Catch::Matchers::WithinAbs(0.20521106432925, 1e-12));
  REQUIRE_THAT(c2b->calc_prime(2.0),
               Catch::Matchers::WithinAbs(-0.1097751763050810, 1e-12));

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

    REQUIRE_NOTHROW(Cut_Poly2());
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

  REQUIRE_THAT(c2b->calc(2.0),
               Catch::Matchers::WithinAbs(1, 1e-12));
  REQUIRE_THAT(c2b->calc_prime(2.0),
               Catch::Matchers::WithinAbs(0, 1e-12));

  REQUIRE_THAT(c2b->calc(3.0),
               Catch::Matchers::WithinAbs(0.31744, 1e-12));
  REQUIRE_THAT(c2b->calc_prime(3.0),
               Catch::Matchers::WithinAbs(-1.728, 1e-12));
    delete c2b;
}
TEST_CASE( "Testing Cutoffs: Cut_Cos_s", "[Cut_Cos_S]" ) {

    REQUIRE_NOTHROW(Cut_Cos_S());
    double rcut2b = 6.2;
    double rcut2bsq = rcut2b*rcut2b;
    using Cut = Cut_Cos_S;
    std::string cuttype = "Cut_Cos_S";
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

  REQUIRE_THAT(c2b->calc(2.0),
               Catch::Matchers::WithinAbs(1, 1e-12));
  REQUIRE_THAT(c2b->calc_prime(2.0),
               Catch::Matchers::WithinAbs(0, 1e-12));

  REQUIRE_THAT(c2b->calc(3.0),
               Catch::Matchers::WithinAbs(0.3454915028, 1e-10));
  REQUIRE_THAT(c2b->calc_prime(3.0),
               Catch::Matchers::WithinAbs(-1.4939160824, 1e-10));
    delete c2b;
}
