#ifndef CORRELATOR_BLAS_H_INCLUDE_GUARD_
#define CORRELATOR_BLAS_H_INCLUDE_GUARD_

#include <complex>

#include "common.h"

/*
  Easy BLAS for C++
*/

// cblas_[sd]scal replacement
template<typename IType>
inline void blas_scal(const int N, 
                      const IType alpha, 
                      IType *X, 
                      const int incX) {
    for(int i=0; i<N; i++) {
        *X *= alpha;
        X += incX;
    }
    
    //cblas_sscal(N, alpha, X, incX);
}

// cbals_[cz]dotc_sub replacement
template<typename IType>
inline void blas_dotc_sub(const int N, 
                          const IType* X, const int incX, 
                          const IType* Y, const int incY,
                          IType* dotc) {
    *dotc = 0.0;
    for(int i=0; i<N; i++) {
        *dotc += conj(*X) * *Y;
        X += incX;
        Y += incY;
    }
    
    //cblas_cdotc_sub(N, X, incX, Y, incY, dotc);
}

#endif // CORRELATOR_BLAS_H_INCLUDE_GUARD_