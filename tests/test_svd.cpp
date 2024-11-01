#include "catch2/catch.hpp"
#include <tadah/core/maths.h>
#include <tadah/models/svd.h>
#include <iomanip>

void multiplyMatrices(const double* A, const double* B, double* C, int m, int n, int p) {
  // C = A * B
  // A is m x n, B is n x p, C is m x p
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < p; ++j) {
      C[i + j * m] = 0;
      for (int k = 0; k < n; ++k) {
        C[i + j * m] += A[i + k * m] * B[k + j * n]; // Column-major order
      }
    }
  }
}


TEST_CASE("SVD Constructor and Destruction") {
  Matrix Phi(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
  SVD svd(Phi);

  REQUIRE(svd.get_info() == 0);
  REQUIRE(svd.sizeU() == 9);
  REQUIRE(svd.sizeVT() == 9);

  // Ensure that singular values are correctly allocated
  double* s = svd.getS();
  for (size_t i = 0; i < svd.sizeS(); ++i) {
    REQUIRE(s[i] >= 0.0);
  }
}

TEST_CASE("SVD Copy Constructor") {
  Matrix Phi(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
  SVD svd(Phi);    
  SVD svd_copy(svd);

  REQUIRE(svd_copy.get_info() == svd.get_info());
  REQUIRE(svd_copy.sizeU() == svd.sizeU());
  REQUIRE(svd_copy.sizeVT() == svd.sizeVT());

  // Verify deep copy
  for (size_t i = 0; i < svd.sizeS(); ++i) {
    REQUIRE(svd_copy.getS()[i] == svd.getS()[i]);
  }
}

TEST_CASE("SVD Numeric Accuracy 2x2") {
  Matrix Phi(2, 2, {1, 2, 3, 4});
  SVD svd(Phi);
  REQUIRE(svd.get_info() == 0);

  double* U = svd.getU();
  double* S = svd.getS();
  double* VT = svd.getVT();

  // S are sorted S[0] > S[1]
  REQUIRE_THAT(5.464985704219043, Catch::Matchers::WithinRel(S[0]));
  REQUIRE_THAT(0.365966190626258, Catch::Matchers::WithinRel(S[1]));
}
TEST_CASE("SVD Numeric Accuracy - Column-Major Order") {
  Matrix Phi(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
  Matrix Phi_cpy(Phi); // as svd destroys original matrix
  SVD svd(Phi);


  REQUIRE(svd.get_info() == 0);

  double* U = svd.getU();
  double* S = svd.getS();
  double* VT = svd.getVT();

  auto [ldu, ucol] = svd.shapeU();
  auto [ldvt, n] = svd.shapeVT();

  Matrix mU (U, ldu,ucol);
  Matrix mVT (VT,ldvt,n);

  aed_type2 vec(S, n);
  REQUIRE_THAT(16.848103352614209, Catch::Matchers::WithinRel(S[0]));
  REQUIRE_THAT(1.068369514554709, Catch::Matchers::WithinRel(S[1]));
  REQUIRE_THAT(0, Catch::Matchers::WithinAbs(S[2],1e-15));

  double USV[9];
  double SV[9];

  // S * VT: We can multiply S into VT directly since S is diagonal
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      SV[i + j * 3] = S[i] * VT[i + j * 3]; // Diagonal multiplication
    }
  }

  multiplyMatrices(U, SV, USV, 3, 3, 3);

  for (int i = 0; i < 9; ++i)
    REQUIRE_THAT(USV[i], Catch::Matchers::WithinRel(Phi_cpy[i]));

}
