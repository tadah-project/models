#include <tadah/models/ols.h>
#include <cmath>
#include <stdexcept>
#include <cstdlib> // For malloc and free
#include <algorithm> // For std::max

// Implementation of non-template methods

// Constructor
OLS::Workspace::Workspace()
  : work(nullptr), lwork(-1), m(0), n(0), algo(Algorithm::GELSD),
  s(nullptr), iwork(nullptr), U(nullptr), S(nullptr), VT(nullptr),
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
    // Allocate workspace for custom SVD implementation using DGESDD
    int minmn = std::min(m, n);
    int maxmn = std::max(m, n);

    // Set jobz parameter to compute singular vectors
    char jobz = 'S'; // Compute the first min(m,n) singular vectors

    // Leading dimensions
    int lda = m;      // Leading dimension of A
    int ldu = m;      // Leading dimension of U
    int ldvt = minmn; // Leading dimension of VT

    // Allocate arrays for U, S, and VT
    U = new double[ldu * minmn];       // U is m x min(m,n)
    S = new double[minmn];             // Singular values vector
    VT = new double[ldvt * n];         // VT is min(m,n) x n

    // Allocate integer workspace for DGESDD
    int* iwork = new int[8 * minmn];

    // Workspace query for DGESDD
    double work_query;
    int lwork = -1;
    int info;

    // Create a dummy A matrix for the workspace query
    double* a_dummy = nullptr;

    // Perform the workspace query
    ::dgesdd_(&jobz, &m, &n, a_dummy, &lda, S, U, &ldu, VT, &ldvt,
              &work_query, &lwork, iwork, &info);

    // Check for errors
    if (info != 0) {
      delete[] iwork;
      throw std::runtime_error("Error in DGESDD workspace query: INFO = " + std::to_string(info));
    }

    // Allocate the optimal workspace
    lwork = static_cast<int>(work_query);
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
  delete[] s;
  delete[] iwork;
  delete[] work;
  delete[] U;
  delete[] S;
  delete[] VT;
  delete[] AtA;
  delete[] Atb;

  s = nullptr;
  iwork = nullptr;
  work = nullptr;
  U = nullptr;
  S = nullptr;
  VT = nullptr;
  AtA = nullptr;
  Atb = nullptr;

  lwork = -1;
  m = 0;
  n = 0;
  algo = Algorithm::GELSD;
}

// is_sufficient method
bool OLS::Workspace::is_sufficient(int m_, int n_, Algorithm algorithm) const {
  if (m != m_ || n != n_ || algo != algorithm) {
    return false;
  }
  // Additional checks can be added here if necessary
  return true;
}

