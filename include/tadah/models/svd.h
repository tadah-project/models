#ifndef SVD_H
#define SVD_H

#include <tadah/core/core_types.h>
#include <tadah/core/maths.h>
#include <tadah/core/lapack.h>

/**
 * @class SVD
 * @brief Computes the Singular Value Decomposition (SVD) of a matrix.
 *
 * Decomposes a matrix into U, S, and V^T such that M = USV^T. Provides methods for accessing
 * components and predicting the effect of this decomposition.
 */
class SVD {
  private:
    char jobz = 'S'; ///< Compute the first min(M,N) columns of U
    int m, n, lda, ldu, ldvt, ucol, lwork, info;

    // Arrays to store SVD components
    double *a = nullptr;
    double *s = nullptr;
    double *u = nullptr;
    double *vt = nullptr;
    double *work = nullptr;
    int *iwork = nullptr;

    bool work_alloc = false;

  public:
    /**
     * @brief Destructor to free allocated resources.
     */
    ~SVD() {
      delete[] s;
      delete[] u;
      delete[] vt;

      if (work_alloc)
        delete[] iwork;

      if (work_alloc && lwork > 0)
        delete[] work;
    }

    /**
     * @brief Copy constructor for SVD class.
     * 
     * Performs a deep copy of another SVD instance.
     * 
     * @param other The SVD instance to copy.
     */
    SVD(const SVD &other) {
      jobz = other.jobz;
      m = other.m;
      n = other.n;
      lda = other.lda;
      ldu = other.ldu;
      ldvt = other.ldvt;
      ucol = other.ucol;
      lwork = other.lwork;
      info = other.info;

      // Allocate memory
      s = new double[sizeS()];
      u = new double[sizeU()];
      vt = new double[sizeVT()];

      // No need for these arrays in a copy
      a = nullptr;
      work = nullptr;
      iwork = nullptr;

      work_alloc = false;

      // Hard copy the data
      std::memcpy(s, other.s, sizeS() * sizeof(double));
      std::memcpy(u, other.u, sizeU() * sizeof(double));
      std::memcpy(vt, other.vt, sizeVT() * sizeof(double));
    }

    /**
     * @brief Constructs the SVD of the given matrix.
     *
     * @param Phi Input matrix to decompose.
     */
    SVD(Matrix &Phi)
      : m(Phi.rows()),
        n(Phi.cols()),
        lda(Phi.rows()),
        ldu(Phi.rows()),
        ldvt(std::min(m, n)),
        ucol(std::min(m, n)),
        a(&Phi.data()[0]),
        s(new double[sizeS()]),
        u(new double[sizeU()]),
        vt(new double[sizeVT()]),
        iwork(new int[8 * std::min(m, n)]),
        work_alloc(true)
    {
      lwork = -1; // Query optimal work size
      double wkopt;
      dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          &wkopt, &lwork, iwork, &info);
      lwork = (int)wkopt;
      work = new double[lwork];
      dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, iwork, &info);

      // Filter out very small singular values
      for (size_t i = 0; i < sizeS(); ++i) {
        if (s[i] < std::numeric_limits<double>::min()) {
          s[i] = 0.0;
        }
      }
    }

    /**
     * @brief Returns the computational info from the LAPACK routine.
     * @return Info status.
     */
    int get_info() const {
      return info;
    }

    /**
     * @brief Accessor for the original matrix A.
     * @return Pointer to matrix A.
     */
    double *getA() const {
      return a;
    }

    /**
     * @brief Accessor for the U matrix.
     * @return Pointer to matrix U.
     */
    double *getU() const {
      return u;
    }

    /**
     * @brief Accessor for the V^T matrix.
     * @return Pointer to matrix V^T.
     */
    double *getVT() const {
      return vt;
    }

    /**
     * @brief Accessor for the singular values matrix S.
     * @return Pointer to singular values.
     */
    double *getS() const {
      return s;
    }

    /**
     * @brief Returns the size of matrix A.
     * @return Size of matrix A.
     */
    size_t sizeA() const {
      return m * n;
    }

    /**
     * @brief Returns the size of matrix U.
     * @return Size of matrix U.
     */
    size_t sizeU() const {
      return ldu * ucol;
    }

    /**
     * @brief Returns the size of matrix V^T.
     * @return Size of matrix V^T.
     */
    size_t sizeVT() const {
      return ldvt * n;
    }

    /**
     * @brief Returns the number of singular values.
     * @return Number of singular values.
     */
    size_t sizeS() const {
      return std::min(m, n);
    }

    /**
     * @brief Returns the dimensions of matrix A.
     * @return Pair of dimensions (rows, cols).
     */
    std::pair<size_t, size_t> shapeA() const {
      return std::pair<size_t, size_t>(m, n);
    }

    /**
     * @brief Returns the dimensions of matrix U.
     * @return Pair of dimensions (rows, cols).
     */
    std::pair<size_t, size_t> shapeU() const {
      return std::pair<size_t, size_t>(ldu, ucol);
    }

    /**
     * @brief Returns the dimensions of matrix V^T.
     * @return Pair of dimensions (rows, cols).
     */
    std::pair<size_t, size_t> shapeVT() const {
      return std::pair<size_t, size_t>(ldvt, n);
    }

    /**
     * @brief Returns the dimensions of the singular values.
     * @return Pair of dimensions (count, 1).
     */
    std::pair<size_t, size_t> shapeS() const {
      return std::pair<size_t, size_t>(std::min(m, n), 1);
    }

    /**
     * @brief Predicts the result of multiplying U, S, V^T with a vector.
     *
     * @tparam T Type of vector w.
     * @param w Vector to be transformed.
     * @return Transformed vector.
     */
    template <typename T>
    T predict(const T &w) const {
      double *t = new double[shapeA().first];
      double *x = new double[shapeA().second];

      T res = predict(w, t, x);

      delete[] t;
      delete[] x;
      return res;
    }

    /**
     * @brief Predicts the result with intermediate storage.
     *
     * @tparam T Type of vector w.
     * @param w Vector to be transformed.
     * @param t Temporary storage for transformation.
     * @param x Temporary storage for transformation.
     * @return Transformed vector.
     */
    template <typename T>
    T predict(const T &w, double *t, double *x) const {
      const double *w_ptr = &w.data()[0];
      char trans = 'N';
      int m = static_cast<int>(shapeA().first);
      int n = static_cast<int>(shapeA().second);
      double alpha_ = 1.0;
      double beta_ = 0.0;
      int incx = 1;

      // 1. V^T w
      dgemv_(&trans, &n, &n, &alpha_, vt, &n, w_ptr, &incx, &beta_, x, &incx);

      // 2. S (V^T w)
      for (int i = 0; i < n; ++i) {
        x[i] *= s[i];
      }

      // 3. U (S V^T w)
      dgemv_(&trans, &m, &n, &alpha_, u, &m, x, &incx, &beta_, t, &incx);
      return T(t, m);
    }
};

#endif // SVD_H
