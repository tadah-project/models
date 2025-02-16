#include "catch2/catch.hpp"
#include <tadah/core/lapack.h>
#include <tadah/core/maths.h>
#include <tadah/models/ols.h>
#include <iomanip>

TEST_CASE("Testing OLS") {
  using Matrix = Matrix_T<>; // Replace with your actual matrix type.
  using Vector = aed_type;   // Replace with your actual vector type.

  // Helper function to compare vectors
  auto vectors_are_close = [](const Vector& v1, const Vector& v2, double tolerance) {
    REQUIRE(v1.size() == v2.size());
    for (size_t i = 0; i < v1.size(); ++i) {
      REQUIRE_THAT(v1[i], Catch::Matchers::WithinAbs(v2[i], tolerance));
    }
  };

  // Prepare test data with corrected matrix initialization
  Matrix A(3, 3, {1, 4, 7,   // First column
    2, 5, 8,   // Second column
    3, 6, 10}); // Third column

  Vector expected_weights = {1, 2, 3};

  // Recompute B = A * expected_weights
  Vector B = A * expected_weights;

  SECTION("Solve using default algorithm (GELSD) without workspace") {
    Vector weights(3);
    OLS::solve(A, B, weights);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  // Other sections remain unchanged...

  SECTION("Test handling of rank-deficient matrix") {
    // Corrected initialization in column-major order
    Matrix A_singular(3, 3, {
      1, 4, 7,  // First column
      2, 5, 8,  // Second column
      3, 6, 9   // Third column
    });

    // Make second column a multiple of the first: second column = 2 * first column
    for (int i = 0; i < 3; ++i) {
      A_singular(i,1) = 2 * A_singular(i,0);
    }

    // Recompute B_singular to match the modified A_singular and expected weights
    Vector expected_weights_singular = {1, 0, 3}; // Since second column is dependent, set weight[1] = 0
    Vector B_singular = A_singular * expected_weights_singular;

    Matrix A_singular2 = A_singular;
    Vector B_singular2 = B_singular;

    Vector weights(3);

    // Solve using GELSD, which can handle rank deficiency
    OLS::solve(A_singular, B_singular, weights);

    // Since the system is rank-deficient, the solution is not unique
    // We can check that A_singular * weights is close to B_singular
    Vector B_computed = A_singular2 * weights;

    // Increase tolerance due to potential numerical inaccuracies
    vectors_are_close(B_computed, B_singular2, 1e-4);
  }

  SECTION("Test solve with overdetermined system") {
    // Corrected matrix initialization in column-major order
    int m = 5;
    int n = 3;
    Matrix A_over(m, n, {
      1, 4, 7, 1, 0, // First column
      2, 5, 8, 0, 1, // Second column
      3, 6, 10,0, 0  // Third column
    });
    Matrix A_over2 = A_over;

    Vector expected_weights_over = {1, 2, 3};

    // Recompute B_over
    Vector B_over = A_over * expected_weights_over;
    Vector B_over2 = B_over;;

    Vector weights(n);

    OLS::solve(A_over, B_over, weights);

    // Compute B_computed = A_over * weights
    Vector B_computed = A_over2 * weights;

    // Check that the residual is small
    for (int i = 0; i < m; ++i) {
      REQUIRE(std::abs(B_over2[i] - B_computed[i]) < 1e-3); // Adjusted tolerance
    }
  }

  SECTION("Test solve with underdetermined system") {
    // Corrected matrix initialization in column-major order
    int m = 3;
    int n = 5;
    Matrix A_under(m, n, {
      1, 6, 11, // First column
      2, 7, 12, // Second column
      3, 8, 13, // Third column
      4, 9, 14, // Fourth column
      5,10,15  // Fifth column
    });

    Vector B_under = {1, 2, 3};
    // Enlarge B_under to size max(m, n) by appending zeros
    B_under.resize(std::max(m, n));

    Vector weights(n);

    OLS::solve(A_under, B_under, weights);

    // Compute B_computed = A_under * weights
    Vector B_computed = A_under * weights;

    // Since B_under has been resized, compare only the first m elements
    Vector B_under_trimmed = Vector(B_under.ptr(),m);
    Vector B_computed_trimmed = Vector(B_computed.ptr(),m);

    vectors_are_close(B_computed_trimmed, B_under_trimmed, 1e-4);
  }

  // Rest of the tests...
}

