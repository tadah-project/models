#ifndef SVD_MULTIPLIER_H
#define SVD_MULTIPLIER_H
#include <tadah/models/svd.h>
#include <tadah/core/lapack.h>
#include <vector>
#include <algorithm>

/**
 * @class SVDMultiplier
 * @brief Class to compute A * x using the SVD decomposition of A without allocating memory during computation.
 */
class SVDMultiplier {
public:
    /**
     * @brief Constructor that initializes the multiplier with the SVD decomposition.
     *
     * @param svd The SVD object containing the decomposition of A.
     */
    explicit SVDMultiplier(const tadah::models::SVD& svd)
        : m_(svd.rows()), n_(svd.cols()), min_mn_(std::min(m_, n_)),
          U_(svd.getU()), S_(svd.getS()), VT_(svd.getVT()),
          y_(min_mn_), z_(min_mn_) {}

    /**
     * @brief Computes the matrix-vector product A * x using the SVD decomposition.
     *
     * @param x Pointer to the input vector x (of size n).
     * @param result Pointer to the output vector where A * x will be stored (of size m).
     */
    void compute(const double* x, double* result) {
        // Step 1: Compute y = V^T * x
        {
            char trans = 'N'; // V^T is already transposed
            int m = min_mn_;
            int n = n_;
            double alpha = 1.0;
            int lda = m;
            int incx = 1;
            double beta = 0.0;
            int incy = 1;

      tadah::core::dgemv_(&trans, &m, &n, &alpha, VT_, &lda, x, &incx, &beta, y_.data(), &incy);
        }

        // Step 2: Compute z = S * y (element-wise multiplication)
        for (int i = 0; i < min_mn_; ++i) {
            z_[i] = S_[i] * y_[i];
        }

        // Step 3: Compute result = U * z
        {
            char trans = 'N';
            int m = m_;
            int n = min_mn_;
            double alpha = 1.0;
            int lda = m;
            int incx = 1;
            double beta = 0.0;
            int incy = 1;

      tadah::core::dgemv_(&trans, &m, &n, &alpha, U_, &lda, z_.data(), &incx, &beta, result, &incy);
        }
    }

private:
    int m_;                  // Number of rows in A
    int n_;                  // Number of columns in A
    int min_mn_;             // Minimum of m and n

    const double* U_;        // Pointer to U matrix
    const double* S_;        // Pointer to singular values
    const double* VT_;       // Pointer to V^T matrix

    std::vector<double> y_;  // Intermediate vector y (size min_mn_)
    std::vector<double> z_;  // Intermediate vector z (size min_mn_)
};
#endif

#ifndef M_EA_H
#define M_EA_H

#include <tadah/core/context.h>
#include <tadah/models/svd.h>
#include <tadah/models/ols.h>
#include <tadah/core/core_types.h>
#include <tadah/core/lapack.h>

namespace tadah {
namespace models {

class EA {
  private:
    const tadah::core::Context &config;
    tadah::models::SVD &svd;
    tadah::core::t_type &T;
    int verbose;
  public:
    EA(const tadah::core::Context &c, tadah::models::SVD &svd, tadah::core::t_type &T):
      config(c),
      svd(svd),
      T(T),
      verbose(c.get<int>("VERBOSE"))
  {}

    int run(double &alpha, double &beta) {
      // Phi = USV^T
      size_t N=std::max(svd.rows(), svd.cols());
      size_t sizeS = std::min(svd.rows(), svd.cols());
      size_t maxsteps=40;
      const double EPS2 = 5e-12;
      int convergence_status = 1;
      double test1 = 2.0*EPS2;
      double test2 = 2.0*EPS2;
      double alpha_old = alpha;
      double beta_old = beta;
      size_t counter=0;
      double lambda = alpha/beta;

      double *s = svd.getS();

      // square of the S matrix contains all eigenvalues
      // of Phi^T*Phi = (VSU^T)(USV^T) = VS^2V^T = S^2
      double *eval = new double[sizeS];
      for (size_t i=0; i<sizeS; ++i)
        eval[i] = s[i]*s[i];

      if (verbose) printf("%9s  %9s  %9s  %9s  %9s  %9s\n","del alpha", "del beta", "alpha", "beta", "lambda", "gamma");
      // D is a filter factors matrix
      double *d = new double[sizeS];

      int outprec=config.get<int>("OUTPREC");

      tadah::core::t_type m_n(sizeS);
      int m = svd.rows();
      int n = svd.cols();
      tadah::core::t_type Tpred(m);
    double rcond = config.size("OALGO")==2 ? config.get<double>("OALGO",1) : 1e-8;
    //
        // Create a workspace manager
        tadah::models::memory::ModelsWorkspaceManager workspaceManager;

        // Create OLS instance with the workspace manager
        tadah::models::OLS ols(workspaceManager);

        SVDMultiplier svdMultiplier(svd);
      while (test1>EPS2 || test2>EPS2) {

        // regularised least square problem (Tikhonov Regularization):
        /*RidgeRegression::solve(svd,T,m_n,lambda,rcond);*/
        /*OLS::solve_with_svd_with_workspace(m, n, T, m_n, svd, ws, rcond, lambda);*/
        ols.solveWithPseudoinverse(m, n, T, m_n, svd.getSVDWorkspace(), rcond, lambda);

        double gamma = 0.0;
        for (size_t j=0; j<sizeS; j++) {
          if (beta*eval[j] > std::numeric_limits<double>::min()) {
            gamma += (beta*eval[j])/(beta*eval[j]+alpha);
          }
        }
        alpha = gamma / (m_n*m_n);

        /*tadah::core::t_type Tpred = svd.predict(m_n,t,x);*/
        svdMultiplier.compute(m_n.data(), Tpred.data());
        double temp = (T-Tpred)*(T-Tpred);

        beta = (static_cast<double>(N)-gamma)/temp;

        test1 = std::fabs(alpha_old-alpha)/alpha;
        test2 = std::fabs(beta_old-beta)/beta;
        if (verbose) printf("%.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n",test1, test2,alpha, beta, alpha/beta, gamma);

        lambda = alpha/beta;
        alpha_old = alpha;
        beta_old = beta;

        if (++counter>maxsteps) {
          convergence_status = 0;
          break;
        }
      }
      if (convergence_status && verbose) std::cout << "Number of steps for convergence: " + tadah::core::to_string(counter,outprec) << std::endl;
      if (!convergence_status && verbose) std::cout << "Max number of steps reached: " + tadah::core::to_string(counter,outprec) << std::endl;

      delete [] eval;
      delete [] d;


      return convergence_status;
    };
};
} // namespace models
} // namespace tadah
#endif

