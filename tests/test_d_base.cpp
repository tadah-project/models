#include "catch2/catch.hpp"
#include <tadah/models/descriptors/d_base.h>
#include <tadah/core/config.h>
#include <tadah/core/core_types.h>

using namespace tadah::core;
using namespace tadah::models;
TEST_CASE( "Testing D_Base", "[D_Base]" ) {

  v_type sgrid(5);
  v_type cgrid(6);
  Config c;
  c.add("SGRID2B", "-1 5 0.1 12");
  c.add("SGRIDMB", "-1 5 0.1 12");
  c.add("CGRID2B", "-1 5 0.1 12");
  c.add("CGRIDMB", "-1 5 0.1 12");
  D_Base::get_grid(c, "SGRID2B", sgrid );
  std::vector<double> sg1 = {0.1,3.075,6.05,9.025,12};
  REQUIRE_THAT(sgrid, Catch::Matchers::Equals(sg1));
  D_Base::get_grid(c, "SGRIDMB", sgrid );


  D_Base::get_grid(c, "CGRID2B", cgrid );
  D_Base::get_grid(c, "CGRIDMB", cgrid );
  std::vector<double> cg1 = {0.1,3.075,6.05,9.025,12};
  REQUIRE_THAT(cgrid, Catch::Matchers::Equals(cg1));

  double margin = 1e-10;
  Config c2;
  c2.add("SGRID2B", "-2 5 0.1 12");
  c2.add("SGRIDMB", "-2 5 0.1 12");
  c2.add("CGRID2B", "-2 5 0.1 12");
  c2.add("CGRIDMB", "-2 5 0.1 12");
  D_Base::get_grid(c2, "SGRID2B", sgrid );
  D_Base::get_grid(c2, "SGRIDMB", sgrid );
  std::vector<double> sg2 = {0.1, 0.330975092, 1.095445115, 3.6256504768, 12.0 };
  REQUIRE_THAT(sgrid, Catch::Matchers::Approx(sg2).margin(margin));

  D_Base::get_grid(c2, "CGRID2B", cgrid );
  D_Base::get_grid(c2, "CGRIDMB", cgrid );
  std::vector<double> cg2 = {0.1, 0.330975092, 1.095445115, 3.6256504768, 12.0 };
  REQUIRE_THAT(cgrid, Catch::Matchers::Approx(cg2).margin(margin));

  Config c3;
  c3.add("SGRID2B", "-2.1 5 0.1 12");
  D_Base::get_grid(c3, "SGRID2B", sgrid );

  //throws
  Config c_err;
  c_err.add("SGRID2B", "-3 5 0.1 12");
  REQUIRE_THROWS(D_Base::get_grid(c_err,"SGRID2B",sgrid));
}
