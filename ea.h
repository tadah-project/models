#ifndef M_EA_H
#define M_EA_H

#include "../LIBS/Eigen/Dense"
#include "../CORE/config/config.h"
#include "../CORE/typedefs.h"

class EA {
    private:
        Config &config;
        int verbose;
    public:
        EA(Config &c):
            config(c),
            verbose(c.get<int>("VERBOSE"))
    {}

        int run  (const mat &Phi,const vec &T,
                const mat &S, const mat &U,const mat &V,
                double &alpha, double &beta) {

            // Phi = USV^T
            size_t N=Phi.rows();
            size_t maxsteps=40;
            const double EPS2 = 5e-12;
            int convergence_status = 1;
            double test1 = 2.0*EPS2;
            double test2 = 2.0*EPS2;
            double alpha_old = alpha;
            double beta_old = beta;
            size_t counter=0;
            double lambda = alpha/beta;

            // square of the S matrix contains all eigenvalues
            // of Phi^T*Phi = (VSU^T)(USV^T) = VS^2V^T = S^2
            const vec eval = S.array().pow(2);


            if (verbose) printf("%9s  %9s  %9s  %9s  %9s  %9s\n","del alpha", "del beta", "alpha", "beta", "lambda", "gamma");
            // D is a filter factors matrix
            mat D(S.rows(),S.rows());
            D.setZero();
            int outprec=config.get<int>("OUTPREC");

            while (test1>EPS2 || test2>EPS2) {

                for (long int i=0; i<S.rows(); ++i) {
                    double s = S(i,0);
                    if (s>std::numeric_limits<double>::min())
                        D(i,i) = (s)/(s*s+lambda);
                    else
                        D(i,i) = 0.0;
                }
                // regularised least square problem (Tikhonov Regularization):
                vec m_n = V*D*U.transpose()*T;

                double gamma = 0.0;
                for (long int j=0; j<S.rows(); j++) {
                    // exlude eigenvalues which are smaller than ~1e-10
                    // for numerical stability b/c matrix is not regularized
                    if (beta*eval(j) > std::numeric_limits<double>::min()) {
                        gamma += (beta*eval(j))/(beta*eval(j)+alpha);
                    }
                }
                alpha = gamma / m_n.dot(m_n);

                double temp = (T-Phi*m_n).dot(T-Phi*m_n);
                beta = (static_cast<double>(N)-gamma)/temp;

                test1 = fabs(alpha_old-alpha)/alpha;
                test2 = fabs(beta_old-beta)/beta;
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
            return convergence_status;

        }
};
#endif

