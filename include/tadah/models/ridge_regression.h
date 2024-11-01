#ifndef RIDGE_REGRESSION_H
#define RIDGE_REGRESSION_H

#include <tadah/core/core_types.h>
#include <tadah/core/maths.h>
#include <tadah/core/lapack.h>
#include <tadah/models/svd.h>

class RidgeRegression {
  public:
    /*template <typename V, typename W>*/
    /*  static void solve(const SVD &svd, V T, W &weights,*/
    /*      const double lambda) {*/
    /**/
    /*    double *u = svd.getU();*/
    /*    double *s = svd.getS();*/
    /*    double *vt = svd.getVT();*/
    /**/
    /*    double *d = new double[svd.sizeS()];*/
    /**/
    /*    for (size_t i=0; i<svd.sizeS(); ++i) {*/
    /*      double sv = s[i];*/
    /*      if (sv>std::numeric_limits<double>::min())*/
    /*        d[i] = (sv)/(sv*sv+lambda);*/
    /*      else*/
    /*        d[i] = 0.0;*/
    /*    }*/
    /**/
    /*    int m = svd.shapeA().first;*/
    /*    int n = svd.shapeA().second;*/
    /*    int ldu = svd.shapeU().first;*/
    /*    int ucol = svd.shapeU().second;*/
    /*    int lda = m;*/
    /*    int ldvt = svd.shapeVT().first;*/
    /**/
    /*    // SOLVE: V D U^T T*/
    /**/
    /*    //STEP1: U^T T*/
    /*    char trans = 'T';*/
    /*    double alpha_ = 1.0;*/
    /*    double *x = T.ptr();*/
    /*    double *y = new double[lda];*/
    /*    int incx = 1;*/
    /*    double beta_ = 0.0;*/
    /*    dgemv_(&trans, &ldu, &ucol, &alpha_, u, &lda, x, &incx,*/
    /*        &beta_, y, &incx);*/
    /**/
    /*    // STEP2 : D (U^T T)*/
    /*    for (int r=0; r<ldvt; ++r) {*/
    /*      y[r]*=d[r];*/
    /*    }*/
    /**/
    /*    // STEP3: V (D U^T T)*/
    /*    trans='T';*/
    /*    dgemv_(&trans, &ldvt, &n, &alpha_, vt, &ldvt, y, &incx,*/
    /*        &beta_, x, &incx);*/
    /**/
    /*    weights = t_type(x, (size_t)ldvt, true);*/
    /**/
    /*    delete[] d;*/
    /*    delete[] y;*/
    /*  }*/

    //template <typename V>
    //static void invertMultiplySigma(const V& sigma, V& result, int n) {
    //  // TODO test s*s!=0
    //  for (int i = 0; i < n; ++i) {
    //    if (sigma[i] != 0.0) {
    //      result[i] = (sigma[i] * sigma[i] ? 1.0 / (sigma[i] * sigma[i]) : 0.0);
    //    } else {
    //      result[i] = 0.0;
    //    }
    //  }
    //}
      template <typename V>
      static void invertMultiplySigmaLambda(const V& sigma, V& result, int n, double lambda) {
        for (int i = 0; i < n; ++i) {
          result[i] = sigma[i] != 0.0 ? sigma[i] / (sigma[i] * sigma[i] + lambda) : 0.0;
        }
      }
    template <typename V, typename W>
      static void solve(const SVD &svd, V b, W &weights,
          const double lambda) {

        double *U = svd.getU(); // m x m
        double *s = svd.getS(); // SIGMA in m x n, here as a vector of n values
        double *VT = svd.getVT(); // n x n

        double *d = new double[svd.sizeS()];

        int m = svd.shapeA().first;
        int n = svd.shapeA().second;

        double alpha_ = 1.0;  // for dgemv_
        double beta_ = 0.0;   // for dgemv_

        // TODO current implemention is without regularization
        /*for (size_t i=0; i<n; ++i) {*/
        /*  double sv = s[i];*/
        /*  if (sv>std::numeric_limits<double>::min())*/
        /*    d[i] = (sv)/(sv*sv+lambda);*/
        /*  else*/
        /*    d[i] = 0.0;*/
        /*}*/


        // SOLVE: V D U^T b

        //STEP1: U^T b
        char trans = 'T';
        double *UTb = new double[n];
        int incx = 1;
        double *bb = b.data();
        dgemv_(&trans, &m, &n, &alpha_, U, &m, b.ptr(), &incx,
            &beta_, UTb, &incx);

        // Step 2: Element-wise multiply and invert sigma in result
        invertMultiplySigmaLambda(s, d, n,lambda);
        for (int i = 0; i < n; ++i) {
          d[i] *= UTb[i];
        }

        trans='T';
        weights.resize(n);
        dgemv_(&trans, &n, &n, &alpha_, VT, &n, d, &incx, &beta_, weights.ptr(), &incx);

        //// STEP2 : D (U^T b)
        //for (int r=0; r<ldvt; ++r) {
        //  UTb[r]*=d[r];
        //}

        //// STEP3: V (D U^T b)
        //trans='T';
        //dgemv_(&trans, &ldvt, &n, &alpha_, vt, &ldvt, y, &incx,
        //    &beta_, x, &incx);

        //weights = t_type(x, (size_t)ldvt, true);
      
        delete[] d;
        delete[] UTb;
      }
};
#endif
