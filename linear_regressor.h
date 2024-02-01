#ifndef LINEAR_REGRESSOR_H
#define LINEAR_REGRESSOR_H

#include "../CORE/config/config.h"
#include "../CORE/eigen_types.h"
#include "../CORE/core_types.h"
#include "../LIBS/Eigen/Dense"
#include "ea.h"
#include <vector>
class LinearRegressor {
            // do training...
            // In general there are 3 methods to solve least squares problem:
            // 1) Normal equation -direct inverse - fast but least accurate:
            //Eigen::MatrixXd I = Eigen::MatrixXd::Identity(Phi.cols(),Phi.cols());
            //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> PhiPhi = Phi.transpose()*Phi;
            //weights = (PhiPhi+lambda*I).inverse() * Phi.transpose() * T;
            // 2) QR - in between 1) and 3) in terms of speed and accuracy
            //weights = Phi.colPivHouseholderQr().solve(T);
            // 3) SVD - most accurate but the slowest
            // w = V*D*S(*)*U^T*T
            // where S(*) = S.asDiagonal().inverse()
            // and D is a filter matrix with diagonal entires: S(i)^2/(S(i)^2+lambda)
    public:
        static void train(Config &config, const phi_type &Phi, const t_type &T,
                t_type &weights, mat &Sigma) {

            // Convert t_type to Eigen type
           vec ww=Eigen::Map<vec>(weights.data.data(),weights.data.size());
           vec TT=Eigen::Map<const vec>(T.data.data(),T.data.size());

            int verbose = config.get<int>("VERBOSE");

            double lambda=config.get<double>("LAMBDA");

            if (lambda == 0)
                ww = Phi.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(TT);
            else {
                double alpha = config.get<double>("ALPHA");
                double beta = config.get<double>("BETA");
                //Eigen::JacobiSVD<mat> svd(Phi, Eigen::ComputeThinU | Eigen::ComputeThinV);
                Eigen::BDCSVD<mat> svd(Phi, Eigen::ComputeThinU | Eigen::ComputeThinV);
                mat S =svd.singularValues();
                const mat &U =svd.matrixU();
                const mat &V =svd.matrixV();

                // Filter out v.small singularvalues from S and its inverse
                for (long int i=0; i<S.rows(); ++i) {
                    if (S(i,0)<std::numeric_limits<double>::min()) {
                        S(i,0) = 0.0;
                    }
                }

                if (lambda<0) {
                    EA ea(config);
                    ea.run(Phi,TT,S,U,V,alpha,beta);
                    lambda=alpha/beta;

                    mat P_inv(S.rows(),S.rows());
                    P_inv.setZero();

                    for (long int i=0; i<S.rows(); ++i) {
                        double s = S(i,0);
                        P_inv(i,i) = 1.0/(beta*s*s+alpha);
                    }

                    Sigma = V*P_inv*V.transpose();
                    config.remove("ALPHA");
                    config.remove("BETA");
                    config.add("ALPHA",alpha);
                    config.add("BETA",beta);
                }
                mat D(S.rows(),S.rows());
                D.setZero();
                for (long int i=0; i<S.rows(); ++i) {
                    double s = S(i,0);
                    if (s>std::numeric_limits<double>::min())
                        D(i,i) = (s)/(s*s+lambda);
                    else
                        D(i,i) = 0.0;
                }
                if (verbose) std::cout << std::endl << "REG LAMBDA: " << lambda << std::endl;
                ww = V*D*U.transpose()*TT;
                //weights = beta*Sigma*V*S.asDiagonal()*U.transpose()*TT;
                // add Sigma to the config file
                dump_sigma(config,Sigma);
            }
            weights = t_type(ww);
        }
        static void dump_sigma(Config &config, mat &Sigma) {
            vec Sigma_as_vector(Eigen::Map<vec>(Sigma.data(), Sigma.cols()*Sigma.rows()));
            config.add("SIGMA",Sigma.rows());
            config.add<vec>("SIGMA",Sigma_as_vector);

        }
        static void read_sigma(Config &config, mat &Sigma) {
            using vec=std::vector<double>;
            // get N = number of rows or cols, Sigma is NxN
            size_t N;
            try {
                N = config.get<size_t>("SIGMA");
            }
            catch (const std::runtime_error& error) {
                throw std::runtime_error("Sigma matrix is not computed. \n\
                        Hint: It is only computed for LAMBDA -1.\n");
            }
            // +1 as first number is N - a number of rows and cols
            vec Sigma_as_vector(N*N+1);
            // read sigma including N, it is sotred col by col
            config.get<vec>("SIGMA",Sigma_as_vector);
            Sigma_as_vector.erase(Sigma_as_vector.begin());
            Sigma.resize(N,N);
            for (size_t c=0; c<N; ++c)
                for (size_t r=0; r<N; ++r)
                    Sigma(r,c) = Sigma_as_vector.at(N*c+r);

        }
};
#endif
