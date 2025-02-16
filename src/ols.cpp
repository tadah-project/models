#include <tadah/models/ols.h>
#include <cmath>
#include <stdexcept>
#include <cstdlib> // For malloc and free

// Ensure that LAPACK functions are declared in the global namespace
// If not, include the necessary headers or declarations

OLS::Workspace::Workspace()
    : s(nullptr), iwork(nullptr), work(nullptr), lwork(0), m(0), n(0), algo(Algorithm::GELSD)
{
}

OLS::Workspace::~Workspace()
{
    release();
}

void OLS::Workspace::allocate(int m_, int n_, Algorithm algorithm)
{
    // Check for valid algorithm
    if (algorithm != Algorithm::GELSD && algorithm != Algorithm::GELS) {
        throw std::invalid_argument("Invalid algorithm specified in Workspace::allocate().");
    }

    // Release existing memory if allocated
    release();

    m = m_;
    n = n_;
    algo = algorithm;

    if (algorithm == Algorithm::GELSD) {
        int minmn = std::min(m, n);

        // Allocate s for singular values
        s = reinterpret_cast<double*>(malloc(minmn * sizeof(double)));
        if (!s) throw std::bad_alloc();

        // Compute nlvl and allocate iwork
        int smlsiz = 25;
        int nlvl = std::max(0, int(std::log2(static_cast<double>(minmn) / (smlsiz + 1))) + 1);
        int iwork_size = 3 * minmn * nlvl + 11 * minmn;
        iwork = reinterpret_cast<int*>(malloc(iwork_size * sizeof(int)));
        if (!iwork) {
            release();
            throw std::bad_alloc();
        }
    } else {
        // For DGELS
        s = nullptr;
        iwork = nullptr;
    }

    // Query for optimal work array size
    double wkopt;
    int lwork_query = -1;
    int nrhs = 1;
    int lda = m;
    int ldb = std::max(m, n);

    // Allocate minimal dummy arrays for A and B
    double* adummy = reinterpret_cast<double*>(malloc(m * n * sizeof(double)));
    double* bdummy = reinterpret_cast<double*>(malloc(ldb * nrhs * sizeof(double)));
    if (!adummy || !bdummy) {
        if (adummy) ::free(adummy);
        if (bdummy) ::free(bdummy);
        release();
        throw std::bad_alloc();
    }

    int info;

    if (algorithm == Algorithm::GELSD) {
        int rank;
        double rcond_query = -1;  // Use default rcond

        // Set up dummy LAPACK call
        ::dgelsd_(&m, &n, &nrhs, adummy, &lda, bdummy, &ldb, s, &rcond_query, &rank,
                &wkopt, &lwork_query, iwork, &info);
    } else { // Algorithm::GELS
        char trans = 'N';
        ::dgels_(&trans, &m, &n, &nrhs, adummy, &lda, bdummy, &ldb, &wkopt, &lwork_query, &info);
    }

    // Free dummy arrays
    ::free(adummy);
    ::free(bdummy);

    if (info != 0) {
        release();
        throw std::runtime_error("Error in workspace query: info = " + std::to_string(info));
    }

    lwork = static_cast<int>(wkopt);
    work = reinterpret_cast<double*>(malloc(lwork * sizeof(double)));
    if (!work) {
        release();
        throw std::bad_alloc();
    }
}

void OLS::Workspace::release()
{
    if (s) ::free(s);
    if (iwork) ::free(iwork);
    if (work) ::free(work);
    s = nullptr;
    iwork = nullptr;
    work = nullptr;
    lwork = 0;
    m = 0;
    n = 0;
    algo = Algorithm::GELSD;
}

bool OLS::Workspace::is_sufficient(int m_, int n_, Algorithm algorithm) const
{
    // Check for valid algorithm
    if (algorithm != Algorithm::GELSD && algorithm != Algorithm::GELS) {
        throw std::invalid_argument("Invalid algorithm specified in Workspace::is_sufficient().");
    }

    // Check if current workspace dimensions and algorithm are sufficient
    if (algo != algorithm || m < m_ || n < n_) {
        return false;
    }

    int maxmn = std::max(m_, n_);

    if (algorithm == Algorithm::GELSD) {
        // Check if s is allocated
        if (!s) return false;

        // For work size, check if lwork is sufficient
        // We can perform a workspace query to find the required lwork
        double wkopt;
        int lwork_query = -1;
        int nrhs = 1;
        int lda = m_;
        int ldb = maxmn;

        // Allocate minimal dummy arrays for A and B
        double* adummy = reinterpret_cast<double*>(malloc(m_ * n_ * sizeof(double)));
        double* bdummy = reinterpret_cast<double*>(malloc(ldb * nrhs * sizeof(double)));
        if (!adummy || !bdummy) {
            if (adummy) ::free(adummy);
            if (bdummy) ::free(bdummy);
            throw std::bad_alloc();
        }
        int rank;
        int info;

        // Query optimal lwork
        double rcond_query = -1;  // Use default rcond
        ::dgelsd_(&m_, &n_, &nrhs, adummy, &lda, bdummy, &ldb, s, &rcond_query, &rank,
                &wkopt, &lwork_query, iwork, &info);

        // Free dummy arrays
        ::free(adummy);
        ::free(bdummy);

        if (info != 0) {
            throw std::runtime_error("Error in workspace sufficiency check: info = " + std::to_string(info));
        }

        int required_lwork = static_cast<int>(wkopt);

        if (lwork < required_lwork) {
            return false;
        }

    } else { // Algorithm::GELS
        // For DGELS
        // Query optimal lwork
        double wkopt;
        int lwork_query = -1;
        int nrhs = 1;
        int lda = m_;
        int ldb = maxmn;
        int info;
        char trans = 'N';

        // Allocate minimal dummy arrays for A and B
        double* adummy = reinterpret_cast<double*>(malloc(m_ * n_ * sizeof(double)));
        double* bdummy = reinterpret_cast<double*>(malloc(ldb * nrhs * sizeof(double)));
        if (!adummy || !bdummy) {
            if (adummy) ::free(adummy);
            if (bdummy) ::free(bdummy);
            throw std::bad_alloc();
        }

        ::dgels_(&trans, &m_, &n_, &nrhs, adummy, &lda, bdummy, &ldb, &wkopt, &lwork_query, &info);

        // Free dummy arrays
        ::free(adummy);
        ::free(bdummy);

        if (info != 0) {
            throw std::runtime_error("Error in workspace sufficiency check: info = " + std::to_string(info));
        }

        int required_lwork = static_cast<int>(wkopt);
        if (lwork < required_lwork) {
            return false;
        }
    }

    // All checks passed
    return true;
}
