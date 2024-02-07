#ifndef EKM_H
#define EKM_H

#include "../CORE/core_types.h"

template <typename K>
class EKM {
    private:
        K kernel;
    public:
        Matrix KK;
        template <typename T>
            EKM(T &t):
                kernel(t)
    {}
        void configure(Matrix &basis) {
            KK = Matrix(basis.cols(), basis.cols());

            // Compute lower diagonal part only
            for (size_t i=0; i<KK.cols(); i++) {
                for (size_t j=i+1; j<KK.cols(); j++) {
                    double val = kernel(basis.col(i),basis.col(j));
                    KK(i,j) = val;
                }
            }

            // KK is a kernel matrix
            // we compute cholesky for KK^-1
            inverse(KK);
            cholesky(KK,'L');

        }
        void inverse(Matrix &A) {
            double *a = &A.data()[0];
            int n = A.rows();   // A is square
            int *ipiv = new int[n];
            int lwork = n*n;
            double *work = new double[lwork];
            int info;

            dgetrf_(&n,&n,a,&n,ipiv,&info);
            dgetri_(&n,a,&n,ipiv,work,&lwork,&info);

            delete[] ipiv;
            delete[] work;
        }
        void cholesky(Matrix &A, char uplo) {
            // matrix must be either U or L
            int n = A.cols(); // A is square
            double *a = &A.data()[0];
            int info;
            dpotrf_(&uplo, &n, a, &n, &info);
        }

        void project(phi_type &Phi) {
            char transa = 'N';
            char transb = 'N';
            int m = Phi.rows();
            int n = Phi.cols();
            double *a = &Phi.data()[0];
            double *b = &KK.data()[0];
            double *c = new double[m*n];
            double alpha = 1.0;
            double beta = 0.0;
            dgemm_(&transa, &transb, &m, &n, &n, &alpha, a, &m, b, &n,
                    &beta, c, &m);

            Phi = Matrix(c,m,n);

            delete[] c;

        }
        aed_type2 project(aed_type2 &aed, Matrix &basis) {
            aed_type2 temp(basis.cols());
            for (size_t b=0; b<basis.cols(); b++) {
                temp(b) = kernel(aed,basis[b]);
            }
            return  KK*temp;
        }
};
#endif
