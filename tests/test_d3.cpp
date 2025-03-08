#include "catch2/catch.hpp"
#include <tadah/models/descriptors/d3/d3_all.h>
#include <tadah/core/maths.h>

using namespace tadah::core;
using namespace tadah::models;
TEST_CASE( "Testing D3_Dummy", "[D3]" ) {
  Config c;
  c.remove("INIT3B");
  c.add("INIT3B", "true");

  REQUIRE_NOTHROW(D3_Dummy());
  REQUIRE_NOTHROW(D3_Dummy(c));

  D3_Dummy d(c);
  REQUIRE(d.size()==0);
  REQUIRE(d.label()=="D3_Dummy");
  aed_type aed;
  fd_type fd;
  REQUIRE_NOTHROW(d.calc_aed(0,0,0,0,0,aed));
  REQUIRE_NOTHROW(d.calc_all(0,0,0,0,0,0.0,0.0,aed,fd));
  REQUIRE_NOTHROW(d.calc_fd(0,0,0,0,0,0.0,0.0,fd));
}
