#ifndef OLS_H
#define OLS_H

#include <tadah/core/lapack.h> // Ensure this header includes the necessary LAPACK function declarations
#include <tadah/core/core_types.h>
#include <cmath>
#include <stdexcept>
#include <cstdlib>    // For malloc and free
#include <algorithm>  // For std::max
#include <iostream>   // For std::ostream

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
 *
 * // Solve using precomputed SVD and regularization
 * SVD svd(A);
 * double lambda = 0.1;
 * OLS::solve_with_svd(A, B, weights, svd, lambda);
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

    // Overload the << operator for the Algorithm enum class
    friend std::ostream& operator<<(std::ostream& os, Algorithm algo) {
        switch (algo) {
            case Algorithm::GELS:
                os << "GELS";
                break;
            case Algorithm::GELSD:
                os << "GELSD";
                break;
            case Algorithm::SVD:
                os << "SVD";
                break;
            case Algorithm::Cholesky:
                os << "Cholesky";
                break;
            default:
                os << "Unknown Algorithm";
                break;
        }
        return os;
    }

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
        double lambda;
        double* U;    // U matrix from SVD (m x min(m,n))
        double* S;    // Singular values array (size min(m,n))
        double* VT;   // V^T matrix from SVD (min(m,n) x n)
        int* iwork_svd; // Integer workspace array for SVD
        double* U_t_b;  // Store result of U^T * b

        // Ownership flags for U, S, VT
        bool owns_U;
        bool owns_S;
        bool owns_VT;

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

    // Original solve methods
    template <typename M, typename V>
    static void solve(M& A, V& B, V& weights);

    template <typename M, typename V>
    static void solve(M& A, V& B, V& weights, double rcond, Algorithm algorithm);

    template <typename M, typename V>
    static void solve_with_workspace(M& A, V& B, V& weights, Workspace& ws);

    template <typename M, typename V>
    static void solve_with_workspace(M& A, V& B, V& weights, Workspace& ws, double rcond, Algorithm algorithm);

    // New methods: solve_with_svd

    // Solve using precomputed SVD (without workspace)
    template <typename V1 ,typename V2, typename SVDType>
      static void solve_with_svd(int m_, int n_, V1& B, V2& weights, SVDType& svd);

    template <typename V1 ,typename V2, typename SVDType>
      static void solve_with_svd(int m_, int n_, V1& B, V2& weights, SVDType& svd, double rcond, double lambda);

    // Solve using precomputed SVD with workspace
    template <typename V1 ,typename V2, typename SVDType>
      static void solve_with_svd_with_workspace(int m_, int n_, V1& B, V2& weights, SVDType& svd, Workspace& ws);

    template <typename V1 ,typename V2, typename SVDType>
      static void solve_with_svd_with_workspace(int m_, int n_, V1& B, V2& weights, SVDType& svd, Workspace& ws, double rcond, double lambda);
};

// Implement template methods here (definitions must be in the header)

// Original solve methods

template <typename M, typename V>
void OLS::solve(M& A, V& B, V& weights) {
    double rcond = -1.0;
    Algorithm algorithm = Algorithm::GELSD;
    solve(A, B, weights, rcond, algorithm);
}

template <typename M, typename V>
void OLS::solve(M& A, V& B, V& weights, double rcond, Algorithm algorithm) {
    Workspace ws;
    solve_with_workspace(A, B, weights, ws, rcond, algorithm);
}

template <typename M, typename V>
void OLS::solve_with_workspace(M& A, V& B, V& weights, Workspace& ws) {
    double rcond = -1.0;
    Algorithm algorithm = Algorithm::GELSD;
    solve_with_workspace(A, B, weights, ws, rcond, algorithm);
}

// The main 'solve_with_workspace' template method
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
        int min_mn = std::min(m, n);

        if (!ws.U || !ws.S || !ws.VT) {
            // Perform SVD

            // Ensure workspace is allocated
            if (!ws.is_sufficient(m, n, algorithm)) {
                ws.allocate(m, n, algorithm);
            }

            char jobz = 'S'; // Compute the first min(m,n) singular values and vectors
            int lda = m;
            int ldu = m;
            int ldvt = min_mn; // VT is min(m,n) x n

            // Now perform SVD
            ::dgesdd_(&jobz, &m, &n, A.ptr(), &lda, ws.S, ws.U, &ldu, ws.VT, &ldvt, ws.work, &ws.lwork, ws.iwork_svd, &info);

            if (info != 0) {
                throw std::runtime_error("Error in DGESDD: info = " + std::to_string(info));
            }

            ws.owns_U = true;
            ws.owns_S = true;
            ws.owns_VT = true;
        }

        if (ws.U_t_b == nullptr) {
            ws.U_t_b = new double[min_mn];
        }

        // Compute Uᵗ * b
        char trans = 'T'; // Transpose U
        double alpha = 1.0;
        double beta = 0.0;
        int incx = 1;
        int incy = 1;
        int lda_u = m;
        ::dgemv_(&trans, &m, &min_mn, &alpha, ws.U, &lda_u, B.ptr(), &incx, &beta, ws.U_t_b, &incy);

        // Find the maximum singular value
        double smax = 0.0;
        for (int i = 0; i < min_mn; ++i) {
            if (ws.S[i] > smax) {
                smax = ws.S[i];
            }
        }

        // Threshold for singular values
        double thresh = rcond > 0.0 ? rcond * smax : 0.0;

        // Apply regularization using smax and lambda from workspace
        for (int i = 0; i < min_mn; ++i) {
            if (ws.S[i] > thresh) {
                ws.U_t_b[i] = ws.U_t_b[i] * ws.S[i] / (ws.S[i] * ws.S[i] + ws.lambda);
            } else {
                ws.U_t_b[i] = 0.0;
            }
        }

        // Compute x = V * (regularized solution)
        trans = 'T'; // Transpose VT to get V
        int lda_vt = min_mn;
        ::dgemv_(&trans, &min_mn, &n, &alpha, ws.VT, &lda_vt, ws.U_t_b, &incx, &beta, weights.ptr(), &incy);

    } else if (algorithm == Algorithm::Cholesky) {
        // Cholesky-based solution
        throw std::runtime_error("Cholesky algorithm is not implemented in this version.");
    }
}

// New methods: solve_with_svd

// Overload without lambda, defaults to lambda = 0.0
template <typename V1 ,typename V2, typename SVDType>
void OLS::solve_with_svd(int m_, int n_, V1& B, V2& weights, SVDType& svd) {
    double rcond = -1.0;
    double lambda = 0.0;
    Workspace ws;
    solve_with_svd_with_workspace(m_, n_, B, weights, svd, ws, rcond, lambda);
}

// Overload with rcond and lambda
template <typename V1 ,typename V2, typename SVDType>
void OLS::solve_with_svd(int m_, int n_, V1& B, V2& weights, SVDType& svd, double rcond, double lambda) {
    Workspace ws;
    solve_with_svd_with_workspace(m_, n_, B, weights, svd, ws, rcond, lambda);
}

// Overload with workspace and default lambda = 0.0
template <typename V1 ,typename V2, typename SVDType>
void OLS::solve_with_svd_with_workspace(int m_, int n_, V1& B, V2& weights, SVDType& svd, Workspace& ws) {
    double rcond = -1.0;
    double lambda = 0.0;
    solve_with_svd_with_workspace(m_, n_, B, weights, svd, ws, rcond, lambda);
}

// Main method with workspace, rcond, and lambda
template <typename V1 ,typename V2, typename SVDType>
void OLS::solve_with_svd_with_workspace(int m_, int n_, V1& B, V2& weights, SVDType& svd, Workspace& ws, double rcond, double lambda) {
    int m = m_;
    int n = n_;

    int min_mn = std::min(m, n);

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
    if (!ws.is_sufficient(m, n, Algorithm::SVD)) {
        ws.allocate(m, n, Algorithm::SVD);
    }

    // Copy precomputed SVD components into the workspace
    ws.U = svd.getU();   // Assume getU() returns a pointer to U (m x min(m,n))
    ws.S = svd.getS();   // Assume getS() returns a pointer to singular values (min(m,n))
    ws.VT = svd.getVT(); // Assume getVT() returns a pointer to VT (min(m,n) x n)
    ws.lambda = lambda;

    // Indicate that we don't own the SVD component pointers
    ws.owns_U = false;
    ws.owns_S = false;
    ws.owns_VT = false;

    // Proceed with the computations
    if (ws.U_t_b == nullptr) {
        ws.U_t_b = new double[min_mn];
    }

    // Compute Uᵗ * B
    char trans = 'T'; // Transpose U
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1;
    int incy = 1;
    int lda = m;
    ::dgemv_(&trans, &m, &min_mn, &alpha, ws.U, &lda, B.ptr(), &incx, &beta, ws.U_t_b, &incy);

    // Find the maximum singular value
    double smax = 0.0;
    for (int i = 0; i < min_mn; ++i) {
        if (ws.S[i] > smax) {
            smax = ws.S[i];
        }
    }

    // Threshold for singular values
    double thresh = rcond > 0.0 ? rcond * smax : 0.0;

    // Apply regularization using 'smax' approach and lambda
    for (int i = 0; i < min_mn; ++i) {
        if (ws.S[i] > thresh) {
            ws.U_t_b[i] = ws.U_t_b[i] * ws.S[i] / (ws.S[i] * ws.S[i] + ws.lambda);
        } else {
            ws.U_t_b[i] = 0.0;
        }
    }

    // Compute x = V * (regularized solution)
    trans = 'T'; // Transpose VT to get V
    int lda_vt = min_mn;
    ::dgemv_(&trans, &min_mn, &n, &alpha, ws.VT, &lda_vt, ws.U_t_b, &incx, &beta, weights.ptr(), &incy);
}

#endif // OLS_H
