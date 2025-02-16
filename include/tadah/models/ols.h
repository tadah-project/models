#ifndef OLS_H
#define OLS_H

#include <tadah/core/lapack.h> // Ensure this header includes the necessary LAPACK function declarations
#include <tadah/core/core_types.h>
#include <cmath>
#include <stdexcept>
#include <cstdlib> // For malloc and free

/**
 * @class OLS
 * @brief Provides functionality for solving Ordinary Least Squares (OLS) problems.
 *
 * Utilizes LAPACK's DGELS and DGELSD functions to solve linear systems where the goal is to minimize 
 * the Euclidean norm of residuals. By default, it uses DGELSD, but users can select DGELS if desired.
 *
 * **Algorithm Selection:**
 * - **GELSD (`DGELSD`):** Uses a singular value decomposition (SVD) approach, which is more robust and
 *   can handle rank-deficient or ill-conditioned systems. It can provide a minimum-norm solution and
 *   handle cases where the matrix `A` does not have full rank. `DGELSD` is recommended when precision
 *   is important, or when the system may be rank-deficient or ill-conditioned.
 * - **GELS (`DGELS`):** Uses QR or LQ factorization, which is generally faster but less robust for
 *   ill-conditioned problems. It assumes that `A` has full rank. `DGELS` is suitable for well-conditioned
 *   systems where performance is critical and the matrix `A` is expected to have full rank.
 *
 * **Usage Recommendation:**
 * - Use **GELSD** (`Algorithm::GELSD`) when dealing with potentially rank-deficient or ill-conditioned systems,
 *   or when you require the most accurate solution at the expense of computational time.
 * - Use **GELS** (`Algorithm::GELS`) when working with well-conditioned systems where speed is more important
 *   than handling rank deficiency.
 *
 * **Workspace Management:**
 * - The class provides methods to solve the OLS problem with and without providing a preallocated workspace.
 *   Preallocating the workspace can improve performance when solving multiple problems of the same size.
 * - If the workspace is insufficient, the `solve_with_workspace` method will automatically reallocate it to
 *   the required size.
 *
 * **Example:**
 * ```cpp
 * Matrix A; // Initialize with your data
 * Vector B; // Initialize with your data
 * Vector weights(A.cols());
 *
 * // Solve using default algorithm (GELSD)
 * OLS::solve(A, B, weights);
 *
 * // Solve using GELS algorithm
 * OLS::solve(A, B, weights, -1.0, OLS::Algorithm::GELS);
 * ```
 */
class OLS {
public:
    /**
     * @brief Enum to select the algorithm used for solving the OLS problem.
     */
    enum class Algorithm {
        GELS,   // Uses DGELS
        GELSD   // Uses DGELSD
    };

    /**
     * @brief Workspace struct to hold pre-allocated memory.
     *
     * This struct contains the necessary buffers for the LAPACK functions.
     * Users can create an instance of `Workspace` and pass it to the `solve_with_workspace` method.
     */
    struct Workspace {
        double* s;    // Singular values array (used by DGELSD)
        int* iwork;   // Integer workspace array (used by DGELSD)
        double* work; // Double workspace array
        int lwork;    // Size of the work array
        int m;        // Number of rows in A (for which workspace was allocated)
        int n;        // Number of columns in A (for which workspace was allocated)
        Algorithm algo; // Algorithm for which workspace was allocated

        Workspace();
        ~Workspace();

        /**
         * @brief Allocates memory for the workspace based on problem dimensions and algorithm.
         * 
         * @param m_ Number of rows in matrix A.
         * @param n_ Number of columns in matrix A.
         * @param algorithm The algorithm for which to allocate workspace.
         */
        void allocate(int m_, int n_, Algorithm algorithm = Algorithm::GELSD);

        /**
         * @brief Frees the allocated memory.
         */
        void release();

        /**
         * @brief Checks if the workspace is sufficient for the given dimensions and algorithm.
         * 
         * @param m_ Number of rows required.
         * @param n_ Number of columns required.
         * @param algorithm The algorithm required.
         * @return True if workspace is sufficient, false otherwise.
         */
        bool is_sufficient(int m_, int n_, Algorithm algorithm = Algorithm::GELSD) const;
    };

    /**
     * @brief Solves the OLS problem for the given matrix and vector without requiring external workspace.
     *
     * @tparam M Type for the matrix A.
     * @tparam V Type for the vectors B and weights.
     * @param A Input matrix (may be modified during computation).
     * @param B Input target vector; contains the solution upon return.
     * @param weights Output vector containing the computed weights.
     * @param rcond The reciprocal of the condition number threshold (used for DGELSD). Default is -1.0.
     * @param algorithm Algorithm to use (GELSD by default).
     */
    template <typename M, typename V>
    static void solve(M& A, V& B, V& weights, double rcond, Algorithm algorithm);

    /**
     * @brief Solves the OLS problem for the given matrix and vector using provided workspace.
     *
     * @tparam M Type for the matrix A.
     * @tparam V Type for the vectors B and weights.
     * @param A Input matrix (may be modified during computation).
     * @param B Input target vector; contains the solution upon return.
     * @param weights Output vector containing the computed weights.
     * @param ws Reference to the Workspace instance.
     * @param rcond The reciprocal of the condition number threshold (used for DGELSD). Default is -1.0.
     * @param algorithm Algorithm to use (GELSD by default).
     */
    template <typename M, typename V>
    static void solve_with_workspace(M& A, V& B, V& weights, Workspace& ws, double rcond, Algorithm algorithm);

    /**
     * @brief Solves the OLS problem using default settings (GELSD algorithm and default rcond).
     *
     * @tparam M Type for the matrix A.
     * @tparam V Type for the vectors B and weights.
     * @param A Input matrix (may be modified during computation).
     * @param B Input target vector; contains the solution upon return.
     * @param weights Output vector containing the computed weights.
     */
    template <typename M, typename V>
    static void solve(M& A, V& B, V& weights);

    /**
     * @brief Solves the OLS problem using default settings and provided workspace.
     *
     * @tparam M Type for the matrix A.
     * @tparam V Type for the vectors B and weights.
     * @param A Input matrix (may be modified during computation).
     * @param B Input target vector; contains the solution upon return.
     * @param weights Output vector containing the computed weights.
     * @param ws Reference to the Workspace instance.
     */
    template <typename M, typename V>
    static void solve_with_workspace(M& A, V& B, V& weights, Workspace& ws);
};

// Implement templated methods here

template <typename M, typename V>
void OLS::solve(M& A, V& B, V& weights, double rcond, Algorithm algorithm)
{
    Workspace ws;
    solve_with_workspace(A, B, weights, ws, rcond, algorithm);
}

template <typename M, typename V>
void OLS::solve_with_workspace(M& A, V& B, V& weights, Workspace& ws, double rcond, Algorithm algorithm)
{
    int m = A.rows();
    int n = A.cols();

    // Check for valid algorithm
    if (algorithm != Algorithm::GELSD && algorithm != Algorithm::GELS) {
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

    int nrhs = 1;
    double* a = A.ptr();
    int lda = m;
    double* b = B.ptr();
    int ldb = maxmn;
    int info;

    if (algorithm == Algorithm::GELSD) {
        int rank;

        // Solve using DGELSD
        ::dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, ws.s, &rcond, &rank,
                ws.work, &ws.lwork, ws.iwork, &info);
    } else { // Algorithm::GELS
        // Solve using DGELS
        char trans = 'N';
        ::dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb,
               ws.work, &ws.lwork, &info);
    }

    if (info != 0) {
        throw std::runtime_error("Error in LAPACK function: info = " + std::to_string(info));
    }

    // Copy solution to weights
    for (int i = 0; i < n; ++i) {
        weights[i] = b[i];
    }
}

template <typename M, typename V>
void OLS::solve(M& A, V& B, V& weights)
{
    // Default rcond and algorithm
    double rcond = -1.0;
    Algorithm algorithm = Algorithm::GELSD;
    solve(A, B, weights, rcond, algorithm);
}

template <typename M, typename V>
void OLS::solve_with_workspace(M& A, V& B, V& weights, Workspace& ws)
{
    // Default rcond and algorithm
    double rcond = -1.0;
    Algorithm algorithm = Algorithm::GELSD;
    solve_with_workspace(A, B, weights, ws, rcond, algorithm);
}

#endif // OLS_H
