#ifndef TADAH_MODELS_OLS_H
#define TADAH_MODELS_OLS_H

#include <tadah/models/Algorithm.h> // Include the Algorithm enum
#include <tadah/models/memory/IModelsWorkspaceManager.h>
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>
#include <iostream>

namespace tadah {
namespace models {

/**
 * @class OLS
 * @brief Provides functionality for solving Ordinary Least Squares (OLS) problems.
 *
 * The OLS class implements various algorithms to solve linear systems where
 * the goal is to minimize the Euclidean norm of residuals.
 *
 * **Algorithms:**
 * - **GELS (`DGELS`):** Uses QR or LQ factorization, fast but less robust for ill-conditioned systems.
 * - **GELSD (`DGELSD`):** Uses Singular Value Decomposition (SVD) approach, handling rank-deficient and ill-conditioned systems; suitable for both overdetermined and underdetermined systems.
 * - **Pseudoinverse:** Custom implementation using SVD to compute the pseudoinverse; allows for explicit control over singular values and regularization; suitable for underdetermined systems.
 * - **NormalEquations:** Solves the normal equations \f$ A^\top A x = A^\top b \f$ using Cholesky decomposition.
 */
class OLS {
public:
    /**
     * @brief Overload the << operator for the Algorithm enum class.
     *
     * @param os The output stream.
     * @param algo The algorithm to output.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, Algorithm algo);

    /**
     * @brief Default constructor.
     *
     * Creates an instance of OLS with its own workspace manager.
     */
    OLS();

    /**
     * @brief Constructor with external workspace manager.
     *
     * @param workspaceManager Reference to an external IModelsWorkspaceManager.
     */
    explicit OLS(memory::IModelsWorkspaceManager& workspaceManager);

    /**
     * @brief Destructor.
     *
     * Deletes the workspace manager if owned.
     */
    ~OLS();

    /**
     * @brief Solve the OLS problem using the default algorithm (GELSD).
     *
     * @tparam M Type of the matrix A.
     * @tparam V Type of the vectors B and weights.
     * @param A The matrix A.
     * @param B The vector B.
     * @param weights Vector to store the solution weights.
     */
    template <typename M, typename V>
    void solve(M& A, V& B, V& weights);

    /**
     * @brief Solve the OLS problem using the specified algorithm.
     *
     * @tparam M Type of the matrix A.
     * @tparam V Type of the vectors B and weights.
     * @param A The matrix A.
     * @param B The vector B.
     * @param weights Vector to store the solution weights.
     * @param rcond Reciprocal condition number; used to determine effective rank.
     * @param algorithm The algorithm to use for solving.
     */
    template <typename M, typename V>
    void solve(M& A, V& B, V& weights, double rcond, Algorithm algorithm);

    /**
     * @brief Solve the OLS problem using a precomputed SVD (Pseudoinverse method).
     *
     * @tparam V1 Type of the vector B.
     * @tparam V2 Type of the vector weights.
     * @param m Number of rows in the matrix A.
     * @param n Number of columns in the matrix A.
     * @param B The vector B.
     * @param weights Vector to store the solution weights.
     * @param svdWorkspace Reference to an SVDWorkspace containing precomputed SVD.
     * @param rcond Reciprocal condition number; used to determine effective rank (default is -1.0).
     */
template <typename V1, typename V2>
void solveWithPseudoinverse(int m, int n, V1& B, V2& weights, memory::SVDWorkspace& svdWorkspace, double rcond, double lambda);
template <typename V1, typename V2>
void solveWithPseudoinverse(int m, int n, V1& B, V2& weights, memory::SVDWorkspace& svdWorkspace);


private:
    memory::IModelsWorkspaceManager* workspaceManager_; ///< Pointer to the workspace manager.
    bool ownsWorkspaceManager_; ///< Indicates if this instance owns the workspace manager.

    // Internal implementation methods

    /**
     * @brief Internal implementation of the solve method.
     *
     * @tparam M Type of the matrix A.
     * @tparam V Type of the vectors B and weights.
     * @param A The matrix A.
     * @param B The vector B.
     * @param weights Vector to store the solution weights.
     * @param rcond Reciprocal condition number; used to determine effective rank.
     * @param algorithm The algorithm to use for solving.
     */
    template <typename M, typename V>
    void solveImpl(M& A, V& B, V& weights, double rcond, Algorithm algorithm);

};

} // namespace models
} // namespace tadah

#include <tadah/models/ols_impl.h> // Include implementation

} // namespace models
} // namespace tadah

namespace tadah {
namespace models {

#endif // TADAH_MODELS_OLS_H

