#include "catch2/catch.hpp"
#include <tadah/models/descriptors/dm/dm_func.h>

#include <cmath>
#include <cassert>

/**
 * @test Checks F_RLR behavior with a simple configuration.
 * Ensures no NaN is returned and verifies plausibility of results.
 */
TEST_CASE("F_RLR: Basic Function and Derivative", "[F_RLR]")
{
  // Reason: "INITMB=true" triggers parameter grid load in F_RLR constructor
  Config c;
  c.remove("INITMB");
  c.add("INITMB", "true");
  c.add("SEMBFUNC", 0.2);
  c.add("CEMBFUNC", 0.3);

  F_RLR rlr(c);

  double rho = 2.0;
  size_t p = 0; // Single parameter index
  double valF = rlr.calc_F(rho, p);
  double valdF = rlr.calc_dF(rho, p);

  // Reason: Checking for no NaN or inf
  REQUIRE_FALSE(std::isnan(valF));
  REQUIRE_FALSE(std::isinf(valF));
  REQUIRE_FALSE(std::isnan(valdF));
  REQUIRE_FALSE(std::isinf(valdF));

  // Basic range checks for sanity (not strict correctness)
  REQUIRE(valF != 0.0);
  REQUIRE(valdF != 0.0);
}

/**
 * @test Checks F_SQRT behavior with a simple configuration.
 * Ensures no NaN is returned and verifies plausibility of results.
 */
TEST_CASE("F_SQRT: Basic Function and Derivative", "[F_SQRT]")
{
  // Reason: "INITMB=true" triggers sgrid load in F_SQRT constructor
  Config c;
  c.remove("INITMB");
  c.add("INITMB", "true");
  c.add("SEMBFUNC", 0.2);

  F_SQRT sqrt_func(c);

  double rho = 4.0;
  size_t p = 0; // Single parameter index
  double valF = sqrt_func.calc_F(rho, p);
  double valdF = sqrt_func.calc_dF(rho, p);

  // Reason: Checking for no NaN or inf
  REQUIRE_FALSE(std::isnan(valF));
  REQUIRE_FALSE(std::isinf(valF));
  REQUIRE_FALSE(std::isnan(valdF));
  REQUIRE_FALSE(std::isinf(valdF));

  // Basic numeric checks
  // s = 0.2, sqrt(4.0) = 2.0, so expected = 0.4 (ideal, ignoring round-off)
  REQUIRE(valF == Approx(0.4).margin(1e-12));

  // Derivative: 0.5*s / sqrt(rho) => 0.5*0.2/2.0 = 0.05
  REQUIRE(valdF == Approx(0.05).margin(1e-12));
}
