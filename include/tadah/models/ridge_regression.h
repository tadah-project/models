#ifndef RIDGE_REGRESSION_H
#define RIDGE_REGRESSION_H

#include <tadah/core/core_types.h>
#include <tadah/core/maths.h>
#include <tadah/core/lapack.h>
#include <tadah/models/svd.h>

/**
 * @class RidgeRegression
 * @brief Implements ridge regression using Singular Value Decomposition (SVD).
 *
 * Provides methods to perform ridge regression, a type of linear regression that includes
 * regularization to prevent overfitting.
 */
class RidgeRegression {
  public:
    /**
     * @brief Inverts and multiplies sigma values by lambda for regularization.
     *
     * Performs an element-wise operation on the sigma values to handle regularization.
     *
     * @tparam V Vector type for sigma and result.
     * @param sigma Input vector of sigma values.
     * @param result Output vector for inverted and multiplied results.
     * @param n Size of the vectors.
     * @param lambda Regularization parameter.
     */
    template <typename V>
    static void invertMultiplySigmaLambda(const V& sigma, V& result, int n, double lambda) {
      for (int i = 0; i < n; ++i) {
        result[i] = sigma[i] != 0.0 ? sigma[i] / (sigma[i] * sigma[i] + lambda) : 0.0;
      }
    }

    /**
     * @brief Solves the ridge regression problem using SVD components.
     *
     * Computes the weights that minimize the regularized least squares error.
     *
     * @tparam V Vector type for b.
     * @tparam W Vector type for weights.
     * @param svd SVD object containing decomposed matrix components.
     * @param b Vector of target values.
     * @param weights Output vector for computed weights.
     * @param lambda Regularization parameter.
     */
    template <typename V, typename W>
    static void solve(const SVD &svd, V b, W &weights, const double lambda) {
      double *U = svd.getU(); // Matrix U from SVD (m x m)
      double *s = svd.getS(); // Singular values (as a vector)
      double *VT = svd.getVT(); // Matrix V^T from SVD (n x n)

      // Dynamic memory allocation for intermediate calculations
      double *d = new double[svd.sizeS()];

      int m = svd.shapeA().first;
      int n = svd.shapeA().second;

      double alpha_ = 1.0;  // Scalar used in calculations
      double beta_ = 0.0;   // Scalar used in calculations

      // Step 1: Compute U^T * b
      char trans = 'T';
      double *UTb = new double[n];
      int incx = 1;
      dgemv_(&trans, &m, &n, &alpha_, U, &m, b.ptr(), &incx, &beta_, UTb, &incx);

      // Step 2: Element-wise multiply and invert sigma in result
      invertMultiplySigmaLambda(s, d, n, lambda);
      for (int i = 0; i < n; ++i) {
        d[i] *= UTb[i];
      }

      // Step 3: Compute V * (D * (U^T * b))
      trans = 'T';
      weights.resize(n);
      dgemv_(&trans, &n, &n, &alpha_, VT, &n, d, &incx, &beta_, weights.ptr(), &incx);

      // Cleanup dynamic memory
      delete[] d;
      delete[] UTb;
    }
};

#endif // RIDGE_REGRESSION_H
