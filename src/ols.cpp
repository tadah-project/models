#include <tadah/models/ols.h>
#include <cmath>
#include <stdexcept>
#include <cstdlib>    // For malloc and free
#include <algorithm>  // For std::max

// Implementation of non-template methods

// Constructor
OLS::Workspace::Workspace()
    : work(nullptr), lwork(-1), m(0), n(0), algo(Algorithm::GELSD),
      s(nullptr), iwork(nullptr),
      lambda(0.0), U(nullptr), S(nullptr), VT(nullptr), iwork_svd(nullptr), U_t_b(nullptr),
      owns_U(false), owns_S(false), owns_VT(false),
      AtA(nullptr), Atb(nullptr) {
}

// Destructor
OLS::Workspace::~Workspace() {
    release();
}

// allocate method
void OLS::Workspace::allocate(int m_, int n_, Algorithm algorithm) {
    release(); // Release any existing memory

    m = m_;
    n = n_;
    algo = algorithm;
    lambda = 0.0; // Reset lambda

    int info;

    if (algorithm == Algorithm::GELSD) {
        // Allocate workspace for DGELSD
        int minmn = std::min(m, n);
        int maxmn = std::max(m, n);
        int nlvl = std::max(0, static_cast<int>(std::log2(minmn / 25.0)) + 1);

        // Allocate s (singular values)
        s = new double[minmn];

        // Allocate iwork
        int liwork = 3 * minmn * nlvl + 11 * minmn;
        iwork = new int[liwork];

        // Workspace query
        double work_query;
        int lwork_query = -1;
        int nrhs = 1;
        double* a_dummy = nullptr;
        double* b_dummy = nullptr;
        ::dgelsd_(&m, &n, &nrhs, a_dummy, &m, b_dummy, &maxmn, s, nullptr,
                  nullptr, &work_query, &lwork_query, iwork, &info);
        lwork = static_cast<int>(work_query);
        work = new double[lwork];
    } else if (algorithm == Algorithm::GELS) {
        // Allocate workspace for DGELS
        double work_query;
        int lwork_query = -1;
        int nrhs = 1;
        double* a_dummy = nullptr;
        double* b_dummy = nullptr;
        char trans = 'N';
        ::dgels_(&trans, &m, &n, &nrhs, a_dummy, &m, b_dummy, &m, &work_query, &lwork_query, &info);
        lwork = static_cast<int>(work_query);
        work = new double[lwork];
    } else if (algorithm == Algorithm::SVD) {

        // Allocate U, S, VT, and iwork_svd
        int min_mn = std::min(m, n);

        // Only allocate if we own the data
        if (owns_U) {
            U = new double[m * min_mn];      // U is m x min(m,n)
        }
        if (owns_S) {
            S = new double[min_mn];          // Singular values
        }
        if (owns_VT) {
            VT = new double[min_mn * n];     // VT is min(m,n) x n
        }
        iwork_svd = new int[8 * min_mn];     // Workspace for dgesdd_

        // Query for optimal workspace size
        char jobz = 'S';
        int lda = m;
        int ldu = m;
        int ldvt = min_mn;
        lwork = -1;
        double wkopt;
        int info_query;
        double* a_dummy = nullptr;
        ::dgesdd_(&jobz, &m, &n, a_dummy, &lda, S, U, &ldu, VT, &ldvt, &wkopt, &lwork, iwork_svd, &info_query);

        if (info_query != 0) {
            throw std::runtime_error("Error in DGESDD (workspace query): info = " + std::to_string(info_query));
        }

        lwork = static_cast<int>(wkopt);
        work = new double[lwork];

    } else if (algorithm == Algorithm::Cholesky) {
        // Allocate space for AtA (n x n) and Atb (n)
        AtA = new double[n * n];
        Atb = new double[n];

        // No additional workspace required for Cholesky decomposition
        work = nullptr;
        lwork = 0;
    } else {
        throw std::invalid_argument("Invalid algorithm specified in Workspace::allocate().");
    }
}

// release method
void OLS::Workspace::release() {
    if (owns_U && U != nullptr) {
        delete[] U;
        U = nullptr;
    }
    if (owns_S && S != nullptr) {
        delete[] S;
        S = nullptr;
    }
    if (owns_VT && VT != nullptr) {
        delete[] VT;
        VT = nullptr;
    }
    delete[] s;
    delete[] iwork;
    delete[] U_t_b;
    delete[] iwork_svd;
    delete[] work;
    delete[] AtA;
    delete[] Atb;

    s = nullptr;
    iwork = nullptr;
    U_t_b = nullptr;
    iwork_svd = nullptr;
    work = nullptr;
    AtA = nullptr;
    Atb = nullptr;

    lwork = -1;
    m = 0;
    n = 0;
    algo = Algorithm::GELSD;
    lambda = 0.0;

    owns_U = false;
    owns_S = false;
    owns_VT = false;
}

// is_sufficient method
bool OLS::Workspace::is_sufficient(int m_, int n_, Algorithm algorithm) const {
    if (m != m_ || n != n_ || algo != algorithm) {
        return false;
    }
    return true;
}
