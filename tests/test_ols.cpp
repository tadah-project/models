#include "catch2/catch.hpp"
#include <tadah/core/maths.h>
#include <tadah/models/ols.h>
#include <tadah/models/ridge_regression.h>
#include <iomanip>

TEST_CASE("Testing OLS") {
  Matrix Phi1(3, 3, {1, 4, 7, 2, 5, 8, 3, 6, 9});
  Matrix Phi2(3, 3, {1, 4, 7, 2, 5, 8, 3, 6, 9});
  aed_type2 b1(3);
  b1[0]=1;
  b1[1]=1;
  b1[2]=1;
  aed_type2 w(3);
  aed_type2 b2=b1;

  OLS::solve(Phi1,b1,w);
  aed_type2 p= Phi2*w;
  REQUIRE_THAT(p[0], Catch::Matchers::WithinRel(b2[0]));
  REQUIRE_THAT(p[1], Catch::Matchers::WithinRel(b2[1]));
  REQUIRE_THAT(p[2], Catch::Matchers::WithinRel(b2[2]));
}
