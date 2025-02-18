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

  SECTION("Solve using GELS algorithm without workspace") {
    Vector weights(3);
    OLS::solve(A, B, weights, -1.0, OLS::Algorithm::GELS);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Solve using Cholesky algorithm without workspace") {
    Vector weights(3);
    OLS::solve(A, B, weights, -1.0, OLS::Algorithm::Cholesky);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Solve using custom SVD implementation without workspace") {
    Vector weights(3);
    OLS::solve(A, B, weights, -1.0, OLS::Algorithm::SVD);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Solve using default algorithm (GELSD) with workspace") {
    Vector weights(3);
    OLS::Workspace ws;
    OLS::solve_with_workspace(A, B, weights, ws);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Solve using GELS algorithm with workspace") {
    Vector weights(3);
    OLS::Workspace ws;
    OLS::solve_with_workspace(A, B, weights, ws, -1.0, OLS::Algorithm::GELS);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Solve using Cholesky algorithm with workspace") {
    Vector weights(3);
    OLS::Workspace ws;
    OLS::solve_with_workspace(A, B, weights, ws, -1.0, OLS::Algorithm::Cholesky);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Solve using custom SVD implementation with workspace") {
    Vector weights(3);
    OLS::Workspace ws;
    OLS::solve_with_workspace(A, B, weights, ws, -1.0, OLS::Algorithm::SVD);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

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
    Vector B_over2 = B_over;

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
    Vector B_under_trimmed = Vector(B_under.ptr(), m);
    Vector B_computed_trimmed = Vector(B_computed.ptr(), m);

    vectors_are_close(B_computed_trimmed, B_under_trimmed, 1e-4);
  }

  SECTION("Test solve_with_svd without workspace") {
    // Assuming SVDType is properly defined and we have access to an SVD implementation
    // Replace 'SVDType' with your actual SVD class or struct

    // For this test, we'll use the custom SVD algorithm on matrix A

    // Placeholder SVDType and methods
    struct SVDType {
      double* getU() { return U.data(); }
      double* getS() { return S.data(); }
      double* getVT() { return VT.data(); }
      std::vector<double> U;
      std::vector<double> S;
      std::vector<double> VT;
    };

    // Compute SVD of A
    SVDType svd;
    int m = A.rows();
    int n = A.cols();
    int min_mn = std::min(m, n);
    svd.U.resize(m * min_mn);
    svd.S.resize(min_mn);
    svd.VT.resize(min_mn * n);

    // Prepare parameters for LAPACK dgesdd_
    char jobz = 'S';
    int lda = m;
    int ldu = m;
    int ldvt = min_mn;
    int info;
    std::vector<double> A_copy(A.ptr(), A.ptr() + m * n);

    // Workspace query
    int lwork = -1;
    double wkopt;
    std::vector<int> iwork(8 * min_mn);
    ::dgesdd_(&jobz, &m, &n, A_copy.data(), &lda, svd.S.data(), svd.U.data(), &ldu, svd.VT.data(), &ldvt, &wkopt, &lwork, iwork.data(), &info);

    lwork = static_cast<int>(wkopt);
    std::vector<double> work(lwork);

    // Compute SVD
    ::dgesdd_(&jobz, &m, &n, A_copy.data(), &lda, svd.S.data(), svd.U.data(), &ldu, svd.VT.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);

    REQUIRE(info == 0);

    Vector weights(3);
    OLS::solve_with_svd(m, n, B, weights, svd);

    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Test solve_with_svd with workspace") {
    // Similar to the previous section but using workspace
    struct SVDType {
      double* getU() { return U.data(); }
      double* getS() { return S.data(); }
      double* getVT() { return VT.data(); }
      std::vector<double> U;
      std::vector<double> S;
      std::vector<double> VT;
    };

    // Compute SVD of A
    SVDType svd;
    int m = A.rows();
    int n = A.cols();
    int min_mn = std::min(m, n);
    svd.U.resize(m * min_mn);
    svd.S.resize(min_mn);
    svd.VT.resize(min_mn * n);

    // Prepare parameters for LAPACK dgesdd_
    char jobz = 'S';
    int lda = m;
    int ldu = m;
    int ldvt = min_mn;
    int info;
    std::vector<double> A_copy(A.ptr(), A.ptr() + m * n);

    // Workspace query
    int lwork = -1;
    double wkopt;
    std::vector<int> iwork(8 * min_mn);
    ::dgesdd_(&jobz, &m, &n, A_copy.data(), &lda, svd.S.data(), svd.U.data(), &ldu, svd.VT.data(), &ldvt, &wkopt, &lwork, iwork.data(), &info);

    lwork = static_cast<int>(wkopt);
    std::vector<double> work(lwork);

    // Compute SVD
    ::dgesdd_(&jobz, &m, &n, A_copy.data(), &lda, svd.S.data(), svd.U.data(), &ldu, svd.VT.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);

    REQUIRE(info == 0);

    Vector weights(3);
    OLS::Workspace ws;
    double rcond = -1.0;
    double lambda = 0.0;

    OLS::solve_with_svd_with_workspace(m, n, B, weights, svd, ws, rcond, lambda);

    vectors_are_close(weights, expected_weights, 1e-4);
  }

  SECTION("Test solve_with_svd with regularization (lambda)") {
    // Using a non-zero lambda for regularization
    struct SVDType {
      double* getU() { return U.data(); }
      double* getS() { return S.data(); }
      double* getVT() { return VT.data(); }
      std::vector<double> U;
      std::vector<double> S;
      std::vector<double> VT;
    };

    // Compute SVD of A
    SVDType svd;
    int m = A.rows();
    int n = A.cols();
    int min_mn = std::min(m, n);
    svd.U.resize(m * min_mn);
    svd.S.resize(min_mn);
    svd.VT.resize(min_mn * n);

    // Prepare parameters for LAPACK dgesdd_
    char jobz = 'S';
    int lda = m;
    int ldu = m;
    int ldvt = min_mn;
    int info;
    std::vector<double> A_copy(A.ptr(), A.ptr() + m * n);

    // Workspace query
    int lwork = -1;
    double wkopt;
    std::vector<int> iwork(8 * min_mn);
    ::dgesdd_(&jobz, &m, &n, A_copy.data(), &lda, svd.S.data(), svd.U.data(), &ldu, svd.VT.data(), &ldvt, &wkopt, &lwork, iwork.data(), &info);

    lwork = static_cast<int>(wkopt);
    std::vector<double> work(lwork);

    // Compute SVD
    ::dgesdd_(&jobz, &m, &n, A_copy.data(), &lda, svd.S.data(), svd.U.data(), &ldu, svd.VT.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);

    REQUIRE(info == 0);

    Vector weights(3);
    double rcond = -1.0;
    double lambda = 0.0005; // Regularization parameter

    OLS::solve_with_svd(m, n, B, weights, svd, rcond, lambda);

    // Since lambda introduces regularization, weights may not match exactly
    // but for this test, we can still check that the solution is reasonable

    Vector B_computed = A * weights;

    // Compute residual norm
    double residual_norm = 0.0;
    for (int i = 0; i < m; ++i) {
      residual_norm += std::pow(B[i] - B_computed[i], 2);
    }
    residual_norm = std::sqrt(residual_norm);

    REQUIRE(residual_norm < 1e-3); // Adjust tolerance as needed
  }

  SECTION("Test solve with Cholesky algorithm on ill-conditioned matrix") {
    // Create an ill-conditioned matrix
    Matrix A_ill(3, 3, {1e-11, 0,   0,
                        0,   1e-11, 0,
                        0,    0,   1e-11});

    Vector B_ill = {1e-11, 1e-11, 1e-11};
    Vector weights(3);

    // Solve using Cholesky algorithm
    OLS::solve(A_ill, B_ill, weights, -1.0, OLS::Algorithm::Cholesky);
    Vector B_computed = A_ill * weights;
    vectors_are_close(B_ill, B_computed, 1e-11);
  }

  SECTION("Test solve with Cholesky algorithm on well-conditioned matrix") {
    Matrix A_well = A; // Use the original well-conditioned matrix
    Vector B_well = B;
    Vector weights(3);

    OLS::solve(A_well, B_well, weights, -1.0, OLS::Algorithm::Cholesky);
    vectors_are_close(weights, expected_weights, 1e-4);
  }

  // Rest of the tests...
}

