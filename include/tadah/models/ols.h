#ifndef OLS_H
#define OLS_H

#include <tadah/core/core_types.h>
#include <tadah/core/lapack.h>

/**
 * @class OLS
 * @brief Provides functionality for solving Ordinary Least Squares (OLS) problems.
 *
 * Utilizes LAPACK's DGELSD function to solve linear systems where the goal is to minimize 
 * the Euclidean norm of residuals.
 */
class OLS {
  public:
    /**
     * @brief Solves the OLS problem for the given matrix and vector.
     *
     * This function computes the weights that minimize the difference between 
     * the target vector and the matrix product.
     *
     * @tparam M Type for the matrix A.
     * @tparam V Type for the vector B.
     * @param A Input matrix, which will be destroyed upon exit.
     * @param B Input target vector, resized if smaller than matrix columns; contains the solution upon return.
     * @param weights Output vector containing the computed weights.
     */
    template <typename M, typename V>
    static void solve(M &A, V &B, V &weights) {
        // Resize B if necessary to match A's column count.
        if (B.size() < A.cols())
            B.resize(A.cols());

        // Variables for LAPACK computation
        int m = A.rows();
        int n = A.cols();
        int nrhs = 1; // Number of target columns
        double *a = A.ptr();
        int lda = m;
        double *b = B.ptr();
        int ldb = std::max(m, n);
        double *s = new double[std::min(m, n)]; // Singular values
        double rcond = 1e-8; // Condition for singularity
        int rank;
        double *work;
        int lwork = -1; // Workspace query
        int smlsiz = 25;
        int nlvl = std::max(0, int(std::log2(std::min(m, n) / (smlsiz + 1))) + 1);
        int *iwork = new int[3 * std::min(m, n) * nlvl + 11 * std::min(m, n)];
        int info;
        double wkopt;

        // Query optimal workspace size
        dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
                &wkopt, &lwork, iwork, &info);
        
        // Allocate workspace memory
        lwork = (int)wkopt;
        work = new double[lwork];

        // Solve the OLS problem
        dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
                work, &lwork, iwork, &info);

        // Copy solution to weights
        weights = B.copy(A.cols());

        // Cleanup memory
        delete[] s;
        delete[] work;
        delete[] iwork;
    }
};

#endif // OLS_H
