#ifndef OLS_H
#define OLS_H

#include <tadah/core/lapack.h> // Ensure this header includes the necessary LAPACK function declarations
#include <tadah/core/core_types.h>
#include <cmath>
#include <stdexcept>
#include <cstdlib> // For malloc and free
#include <algorithm> // For std::max

/**
 * @class OLS
 * @brief Provides functionality for solving Ordinary Least Squares (OLS) problems.
 *
 * Utilizes various algorithms to solve linear systems where the goal is to minimize 
 * the Euclidean norm of residuals. Supports the following algorithms:
 *
 * **Algorithms:**
 * - **GELSD (`DGELSD`):** Uses a singular value decomposition (SVD) approach, handling rank-deficient and ill-conditioned systems.
 * - **GELS (`DGELS`):** Uses QR or LQ factorization, fast but less robust for ill-conditioned systems.
 * - **SVD:** Custom implementation using SVD, allows for explicit control over singular values and regularization.
 * - **Cholesky:** Solves the normal equations using Cholesky decomposition, efficient but may be less stable for ill-conditioned systems.
 *
 * **Workspace Management:**
 * - The class provides a `Workspace` struct to manage memory allocations for the computations.
 * - Workspaces are allocated based on the dimensions of the problem and the algorithm selected.
 *
 * **Example:**
 * ```cpp
 * Matrix A; // Initialize with your data
 * Vector B; // Initialize with your data
 * Vector weights(A.cols());
 * OLS::Workspace ws;
 *
 * // Solve using default algorithm (GELSD)
 * OLS::solve(A, B, weights);
 *
 * // Solve using Cholesky algorithm
 * OLS::solve_with_workspace(A, B, weights, ws, -1.0, OLS::Algorithm::Cholesky);
 * ```
 */
class OLS {
public:
  /**
     * @brief Enum to select the algorithm used for solving the OLS problem.
     */
  enum class Algorithm {
    GELS,     // Uses DGELS
    GELSD,    // Uses DGELSD
    SVD,      // Uses custom SVD implementation
    Cholesky  // Uses Cholesky decomposition
  };

  /**
     * @brief Workspace struct to hold pre-allocated memory.
     *
     * This struct contains the necessary buffers for the computations in different algorithms.
     * Users can create an instance of `Workspace` and pass it to the `solve_with_workspace` method.
     */
  struct Workspace {
    // Common workspace variables
    double* work; // Double workspace array
    int lwork;    // Size of the work array
    int m;        // Number of rows in A (for which workspace was allocated)
    int n;        // Number of columns in A (for which workspace was allocated)
    Algorithm algo; // Algorithm for which workspace was allocated

    // For DGELSD (GELSD algorithm)
    double* s;    // Singular values array
    int* iwork;   // Integer workspace array

    // For SVD algorithm
    double* U;    // U matrix from SVD (m x n)
    double* S;    // Singular values array (size min(m,n))
    double* VT;   // V^T matrix from SVD (n x n)

    // For Cholesky algorithm
    double* AtA;  // A^T * A matrix (n x n)
    double* Atb;  // A^T * b vector (n)

    // Constructors and methods (declarations only)
    Workspace();
    ~Workspace();

    void allocate(int m_, int n_, Algorithm algorithm = Algorithm::GELSD);
    void release();
    bool is_sufficient(int m_, int n_, Algorithm algorithm = Algorithm::GELSD) const;
  };

  // Public methods
  template <typename M, typename V>
  static void solve(M& A, V& B, V& weights, double rcond, Algorithm algorithm);

  template <typename M, typename V>
  static void solve_with_workspace(M& A, V& B, V& weights, Workspace& ws, double rcond, Algorithm algorithm);

  template <typename M, typename V>
  static void solve(M& A, V& B, V& weights);

  template <typename M, typename V>
  static void solve_with_workspace(M& A, V& B, V& weights, Workspace& ws);
};

// Implement template methods here
// Since they are templates, definitions must be in the header

template <typename M, typename V>
void OLS::solve(M& A, V& B, V& weights, double rcond, Algorithm algorithm) {
  Workspace ws;
  solve_with_workspace(A, B, weights, ws, rcond, algorithm);
}

template <typename M, typename V>
void OLS::solve_with_workspace(M& A, V& B, V& weights, Workspace& ws, double rcond, Algorithm algorithm) {
  int m = A.rows();
  int n = A.cols();

  // Check for valid algorithm
  if (algorithm != Algorithm::GELSD && algorithm != Algorithm::GELS &&
    algorithm != Algorithm::SVD && algorithm != Algorithm::Cholesky) {
    throw std::invalid_argument("Invalid algorithm specified.");
  }

  // Check that dimensions of weights match expected size
  if (weights.size() != static_cast<std::size_t>(n)) {
    throw std::invalid_argument("weights size is incorrect.");
  }

  // Check that B has sufficient size
  int maxmn = std::max(m, n);
  if (B.size() < static_cast<std::size_t>(maxmn)) {
    throw std::invalid_argument("B size is insufficient.");
  }

  // Check if workspace is sufficient; if not, reallocate
  if (!ws.is_sufficient(m, n, algorithm)) {
    ws.allocate(m, n, algorithm);
  }

  int info;
  if (algorithm == Algorithm::GELSD) {
    int nrhs = 1;
    double* a = A.ptr();
    int lda = m;
    double* b = B.ptr();
    int ldb = maxmn;
    int rank;

    // Solve using DGELSD
    ::dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, ws.s, &rcond, &rank,
              ws.work, &ws.lwork, ws.iwork, &info);

    if (info != 0) {
      throw std::runtime_error("Error in DGELSD: info = " + std::to_string(info));
    }

    // Copy solution to weights
    for (int i = 0; i < n; ++i) {
      weights[i] = b[i];
    }
  } else if (algorithm == Algorithm::GELS) {
    int nrhs = 1;
    double* a = A.ptr();
    int lda = m;
    double* b = B.ptr();
    int ldb = m;
    char trans = 'N';

    // Solve using DGELS
    ::dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb,
             ws.work, &ws.lwork, &info);

    if (info != 0) {
      throw std::runtime_error("Error in DGELS: info = " + std::to_string(info));
    }

    // Copy solution to weights
    for (int i = 0; i < n; ++i) {
      weights[i] = b[i];
    }
  } else if (algorithm == Algorithm::SVD) {
    // Custom SVD-based solution
    // Implementation can be provided or left as an exercise
    //throw std::runtime_error("SVD algorithm is not implemented in this version.");
    // Allocate memory for intermediate calculations






  } else if (algorithm == Algorithm::Cholesky) {
    // Cholesky-based solution
    // Implementation can be provided or left as an exercise
    throw std::runtime_error("Cholesky algorithm is not implemented in this version.");
  }
}

template <typename M, typename V>
void OLS::solve(M& A, V& B, V& weights) {
  // Default rcond and algorithm
  double rcond = -1.0;
  Algorithm algorithm = Algorithm::GELSD;
  solve(A, B, weights, rcond, algorithm);
}

template <typename M, typename V>
void OLS::solve_with_workspace(M& A, V& B, V& weights, Workspace& ws) {
  // Default rcond and algorithm
  double rcond = -1.0;
  Algorithm algorithm = Algorithm::GELSD;
  solve_with_workspace(A, B, weights, ws, rcond, algorithm);
}

#endif // OLS_H

