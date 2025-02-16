#ifndef RIDGE_REGRESSION_H
#define RIDGE_REGRESSION_H

#include <tadah/core/core_types.h>
#include <tadah/core/maths.h>
#include <tadah/core/lapack.h>
#include <tadah/models/svd.h>
#include <algorithm> // For std::max_element

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
     * @brief Inverts and multiplies sigma values by lambda for regularization,
     *        considering the rcond threshold.
     *
     * Performs an element-wise operation on the sigma values to handle regularization.
     * Singular values below the threshold are discarded (set to zero).
     *
     * @tparam V Vector type for sigma and result.
     * @param sigma Input vector of sigma values.
     * @param result Output vector for inverted and multiplied results.
     * @param n Size of the vectors.
     * @param lambda Regularization parameter.
     * @param rcond Reciprocal of the condition number; threshold for singular values.
     */
    template <typename V>
    static void invertMultiplySigmaLambdaRcond(const V& sigma, V& result, int n, double lambda, double rcond) {
      // Find the maximum singular value
      double smax = 0.0;
      for (int i = 0; i < n; ++i) {
        if (sigma[i] > smax) {
          smax = sigma[i];
        }
      }

      // Threshold for singular values
      double thresh = rcond > 0.0 ? rcond * smax : 0.0;

      for (int i = 0; i < n; ++i) {
        if (sigma[i] > thresh) {
          result[i] = sigma[i] / (sigma[i] * sigma[i] + lambda);
        } else {
          // Singular value is too small; set result to zero
          result[i] = 0.0;
        }
      }
    }

    /**
     * @brief Solves the ridge regression problem using SVD components,
     *        incorporating rcond for singular value thresholding.
     *
     * Computes the weights that minimize the regularized least squares error.
     *
     * @tparam V Vector type for b.
     * @tparam W Vector type for weights.
     * @param svd SVD object containing decomposed matrix components.
     * @param b Vector of target values.
     * @param weights Output vector for computed weights.
     * @param lambda Regularization parameter.
     * @param rcond Reciprocal of the condition number; threshold for singular values.
     */
    template <typename V, typename W>
    static void solve(const SVD &svd, V b, W &weights, double lambda, double rcond) {
      double *U = svd.getU();    // Matrix U from SVD (m x n)
      double *s = svd.getS();    // Singular values (as a vector)
      double *VT = svd.getVT();  // Matrix V^T from SVD (n x n)

      int m = svd.shapeA().first;
      int n = svd.shapeA().second;

      // Allocate memory for intermediate calculations
      double *d = new double[n];   // For inverted singular values
      double *UTb = new double[n]; // For U^T * b

      // Step 1: Compute U^T * b
      char trans = 'T';
      double alpha_ = 1.0;  // Scalar for multiplication
      double beta_ = 0.0;   // Scalar for addition
      int incx = 1;
      int incy = 1;
      dgemv_(&trans, &m, &n, &alpha_, U, &m, b.ptr(), &incx, &beta_, UTb, &incy);

      // Step 2: Invert and multiply sigma values, applying rcond threshold
      invertMultiplySigmaLambdaRcond(s, d, n, lambda, rcond);

      // Element-wise multiplication: d[i] = d[i] * (U^T * b)[i]
      for (int i = 0; i < n; ++i) {
        d[i] *= UTb[i];
      }

      // Step 3: Compute weights = V * d
      trans = 'T'; // Since VT is V^T, transposing VT gives V
      weights.resize(n);
      dgemv_(&trans, &n, &n, &alpha_, VT, &n, d, &incx, &beta_, weights.ptr(), &incy);

      // Cleanup dynamic memory
      delete[] d;
      delete[] UTb;
    }
};

#endif // RIDGE_REGRESSION_H
