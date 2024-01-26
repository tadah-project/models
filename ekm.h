#ifndef EKM_H
#define EKM_H

#include "../LIBS/Eigen/Dense"
#include "../CORE/typedefs.h"

template <typename K>
class EKM {
    private:
        K kernel;
    public:
        Eigen::MatrixXd KK;
        template <typename T>
            EKM(T &t):
                kernel(t)
    {}
        void configure(aeds_type &basis) {
            // KK is temp matrix
            KK = Eigen::MatrixXd(basis.size(), basis.size());

            for (long int i=0; i<KK.rows(); i++) {
                for (long int j=0; j<KK.cols(); j++) {
                    double val = kernel(basis.col(i),basis.col(j));
                    KK(i,j) = val;
                }
            }

            // KK is a kernel matrix
            // we compute cholesky for KK^-1
            // then KK = L^-1
            KK = KK.inverse();
            Eigen::LLT<Eigen::MatrixXd> llt(KK);
            KK = llt.matrixL();
            KK.transposeInPlace();
        }

        void project(phi_type &Phi) {
            for (long int i=0; i<Phi.rows(); ++i) {
                Phi.row(i)= KK*Phi.row(i).transpose();
            }
        }
        aed_type project(aed_rtype aed, aeds_type &basis) {
            aed_type temp(basis.size());
            for (long int b=0; b<basis.cols(); b++) {
                temp(b) = kernel(aed,basis[b]);
            }
            return  KK*temp;
        }
};
#endif
