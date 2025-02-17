#ifndef M_EA_H
#define M_EA_H

#include <tadah/core/config.h>
#include <tadah/models/svd.h>
#include <tadah/models/ridge_regression.h>
#include <tadah/models/ols.h>
#include <tadah/core/core_types.h>
#include <tadah/core/lapack.h>

class EA {
  private:
    const Config &config;
    const SVD &svd;
    t_type &T;
    int verbose;
  public:
    EA(const Config &c, const SVD &svd, t_type &T):
      config(c),
      svd(svd),
      T(T),
      verbose(c.get<int>("VERBOSE"))
  {}

    int run(double &alpha, double &beta) {
      // Phi = USV^T
      size_t N=svd.shapeA().first;
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
      double *eval = new double[svd.sizeS()];
      for (size_t i=0; i<svd.sizeS(); ++i)
        eval[i] = s[i]*s[i];

      if (verbose) printf("%9s  %9s  %9s  %9s  %9s  %9s\n","del alpha", "del beta", "alpha", "beta", "lambda", "gamma");
      // D is a filter factors matrix
      double *d = new double[svd.sizeS()];

      int outprec=config.get<int>("OUTPREC");

      double *t = new double[svd.shapeA().first];
      double *x = new double[svd.shapeA().second];
      t_type m_n((size_t)svd.sizeS());
      OLS::Workspace ws;
      auto shapeA = svd.shapeA();
      int m = shapeA.first;
      int n = shapeA.second;
    double rcond = config.size("OALGO")==2 ? config.get<double>("OALGO",1) : 1e-8;
      while (test1>EPS2 || test2>EPS2) {

        // regularised least square problem (Tikhonov Regularization):
        /*RidgeRegression::solve(svd,T,m_n,lambda,rcond);*/
        OLS::solve_with_svd_with_workspace(m, n, T, m_n, svd, ws, rcond, lambda);

        double gamma = 0.0;
        for (size_t j=0; j<svd.sizeS(); j++) {
          if (beta*eval[j] > std::numeric_limits<double>::min()) {
            gamma += (beta*eval[j])/(beta*eval[j]+alpha);
          }
        }
        alpha = gamma / (m_n*m_n);

        t_type Tpred = svd.predict(m_n,t,x);
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
      if (convergence_status && verbose) std::cout << "Number of steps for convergence: " + to_string(counter,outprec) << std::endl;
      if (!convergence_status && verbose) std::cout << "Max number of steps reached: " + to_string(counter,outprec) << std::endl;

      delete [] eval;
      delete [] d;

      delete[] t;
      delete[] x;

      return convergence_status;
    };
};
#endif

