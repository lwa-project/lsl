#ifndef CORRELATOR_BLAS_H_INCLUDE_GUARD_
#define CORRELATOR_BLAS_H_INCLUDE_GUARD_

#include <complex.h>
extern "C" {
    #include <cblas.h>
}

/*
  Easy BLAS for C++
*/

// cblas_[sd]scal wrappers
inline void blas_scal(const int N, 
                      const float alpha, 
                      float *X, 
                      const int incX) {
    cblas_sscal(N, alpha, X, incX);
}
inline void blas_scal(const int N, 
                      const double alpha, 
                      double *X, 
                      const int incX) {
    cblas_dscal(N, alpha, X, incX);
}

// cbals_[cz]dotc_sub wrappers
inline void blas_dotc_sub(const int N, 
                          const float complex* X, const int incX, 
                          const float complex* Y, const int incY,
                          float complex* dotc) {
    cblas_cdotc_sub(N, X, incX, Y, incY, dotc);
}
inline void blas_dotc_sub(const int N, 
                          const double complex* X, const int incX, 
                          const double complex* Y, const int incY,
                          double complex* dotc) {
    cblas_zdotc_sub(N, X, incX, Y, incY, dotc);
}


#endif // CORRELATOR_BLAS_H_INCLUDE_GUARD_