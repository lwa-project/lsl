#ifndef COMPLEX_COMPLEX_INT32_H_INCLUDE_GUARD_
#define COMPLEX_COMPLEX_INT32_H_INCLUDE_GUARD_

#ifdef __cplusplus
extern "C" {
#endif

#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
    
typedef struct {
    short int real;
    short int imag;
} complex_int32;

// Access methods
static NPY_INLINE void lsl_unpack_ci32(complex_int32 packed, short int* real, short int* imag) {
    *real = packed.real;
    *imag = packed.imag;
}

static NPY_INLINE void lsl_pack_ci32(short int real, short int imag, complex_int32* packed) {
    packed->real = real;
    packed->imag = imag;
}

// Unary bool operators
static NPY_INLINE int complex_int32_nonzero(complex_int32 c) {
    return c.real != 0 || c.imag != 0;
}

static NPY_INLINE int complex_int32_isnan(complex_int32 c) {
    return 0;
}

static NPY_INLINE int complex_int32_isinf(complex_int32 c) {
    return 0;
}

static NPY_INLINE int complex_int32_isfinite(complex_int32 c) {
    return 1;
}

// Binary bool operators
static NPY_INLINE int complex_int32_equal(complex_int32 c1, complex_int32 c2) {
    return 
        !complex_int32_isnan(c1) &&
        !complex_int32_isnan(c2) &&
        c1.real == c2.real &&
        c1.imag == c2.imag;
}

static NPY_INLINE int complex_int32_not_equal(complex_int32 c1, complex_int32 c2) {
    return !complex_int32_equal(c1, c2);
}

static NPY_INLINE int complex_int32_less(complex_int32 c1, complex_int32 c2) {
    return
        (!complex_int32_isnan(c1) &&
         !complex_int32_isnan(c2)) && (
            c1.real != c2.real ? c1.real < c2.real :
            c1.imag != c2.imag ? c1.imag < c2.imag : 0);
}

static NPY_INLINE int complex_int32_greater(complex_int32 c1, complex_int32 c2) {
    return
        (!complex_int32_isnan(c1) &&
         !complex_int32_isnan(c2)) && (
            c1.real != c2.real ? c1.real > c2.real :
            c1.imag != c2.imag ? c1.imag > c2.imag : 0);
}

static NPY_INLINE int complex_int32_less_equal(complex_int32 c1, complex_int32 c2) {
    return
        (!complex_int32_isnan(c1) &&
         !complex_int32_isnan(c2)) && (
            c1.real != c2.real ? c1.real < c2.real :
            c1.imag != c2.imag ? c1.imag < c2.imag : 1);
}

static NPY_INLINE int complex_int32_greater_equal(complex_int32 c1, complex_int32 c2) {
    return
        (!complex_int32_isnan(c1) &&
         !complex_int32_isnan(c2)) && (
            c1.real != c2.real ? c1.real > c2.real :
            c1.imag != c2.imag ? c1.imag > c2.imag : 1);
}

// Unary int returners
static NPY_INLINE long complex_int32_real(complex_int32 c) {
    return (long) c.real;
}

static NPY_INLINE long complex_int32_imag(complex_int32 c) {
    return (long) c.imag;
}

// Unary float returners
static NPY_INLINE double complex_int32_absolute(complex_int32 c) {
    return sqrt(((int) c.real)*c.real + ((int) c.imag)*c.imag);
}

static NPY_INLINE double complex_int32_kludgy_arctan2(long i, complex_int32 c) {
    return atan2((int) c.imag, (int) c.real);
}

// Unary complex_int32 returners
static NPY_INLINE complex_int32 complex_int32_negative(complex_int32 c) {
    short int real = -c.real;
    short int imag = -c.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE complex_int32 complex_int32_conjugate(complex_int32 c) {
    short int real =  c.real;
    short int imag = -c.imag;
    return (complex_int32) {real, imag};
}

// complex_int32-complex_int32/complex_int32-scalar/scalar-complex_int32 binary complex_int32 returners
static NPY_INLINE complex_int32 complex_int32_add(complex_int32 c1, complex_int32 c2) {
    short int real = c1.real + c2.real;
    short int imag = c1.imag + c2.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_add(complex_int32* c1, complex_int32 c2) {
    short int real = c1->real + c2.real;
    short int imag = c1->imag + c2.imag;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_scalar_add(long s, complex_int32 c) {
    short int real = s + c.real;
    short int imag = 0 + c.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_scalar_add(long s, complex_int32* c) {
    short int real = s + c->real;
    short int imag = s + c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_add_scalar(complex_int32 c, long s) {
    short int real = s + c.real;
    short int imag = 0 + c.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_add_scalar(complex_int32* c, long s) {
    short int real = s + c->real;
    short int imag = 0 + c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_subtract(complex_int32 c1, complex_int32 c2) {
    short int real = c1.real - c2.real;
    short int imag = c1.imag - c2.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_subtract(complex_int32* c1, complex_int32 c2) {
    short int real = c1->real - c2.real;
    short int imag = c1->imag - c2.imag;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_scalar_subtract(long s, complex_int32 c) {
    short int real = s - c.real;
    short int imag = 0 - c.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_scalar_subtract(long s, complex_int32* c) {
    short int real = s - c->real;
    short int imag = s - c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_subtract_scalar(complex_int32 c, long s) {
    short int real = -s + c.real;
    short int imag = -0 + c.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_subtract_scalar(complex_int32* c, long s) {
    short int real = -s + c->real;
    short int imag = -s + c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_multiply(complex_int32 c1, complex_int32 c2) {
    short int real = c1.real*c2.real - c1.imag*c2.imag;
    short int imag = c1.imag*c2.real + c1.real*c2.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_multiply(complex_int32* c1, complex_int32 c2) {
    short int real = c1->real*c2.real - c1->imag*c2.imag;
    short int imag = c1->imag*c2.real + c1->real*c2.imag;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_scalar_multiply(long s, complex_int32 c) {
    short int real = s*c.real - 0*c.imag;
    short int imag = 0*c.real + s*c.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_scalar_multiply(long s, complex_int32* c) {
    short int real = s*c->real - 0*c->imag;
    short int imag = 0*c->real + s*c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_multiply_scalar(complex_int32 c, long s) {
    short int real = s*c.real - 0*c.imag;
    short int imag = 0*c.real + s*c.imag;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_multiply_scalar(complex_int32* c, long s) {
    short int real = s*c->real - 0*c->imag;
    short int imag = 0*c->real + s*c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_divide(complex_int32 c1, complex_int32 c2) {
    long mag2 = ((int) c2.real)*c2.real + ((int) c2.imag)*c2.imag;
    short int real = (c1.real*c2.real + c1.imag*c2.imag) / mag2;
    short int imag = (c1.imag*c2.real - c1.real*c2.imag) / mag2;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_divide(complex_int32* c1, complex_int32 c2) {
    long mag2 = ((int) c2.real)*c2.real + ((int) c2.imag)*c2.imag;
    short int real = (c1->real*c2.real + c1->imag*c2.imag) / mag2;
    short int imag = (c1->imag*c2.real - c1->real*c2.imag) / mag2;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_scalar_divide(long s, complex_int32 c) {
    long mag2 = ((int) c.real)*c.real + ((int) c.imag)*c.imag;
    short int real = (s*c.real + 0*c.imag) / mag2;
    short int imag = (s*c.real - 0*c.imag) / mag2;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_scalar_divide(long s, complex_int32* c) {
    long mag2 = ((int) c->real)*c->real + ((int) c->imag)*c->imag;
    short int real = (s*c->real + 0*c->imag) / mag2;
    short int imag = (s*c->real - 0*c->imag) / mag2;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int32 complex_int32_divide_scalar(complex_int32 c, long s) {
    long mag2 = s*s + 0*0;
    short int real = (c.real*s - c.imag*0) / mag2;
    short int imag = (c.imag*s + c.real*0) / mag2;
    return (complex_int32) {real, imag};
}

static NPY_INLINE void complex_int32_inplace_divide_scalar(complex_int32* c, long s) {
    long mag2 = s*s + 0*0;
    short int real = (c->real*s - c->imag*0) / mag2;
    short int imag = (c->imag*s + c->real*0) / mag2;
    c->real = real;
    c->imag = imag;
}

#ifdef __cplusplus
}
#endif

#endif // COMPLEX_COMPLEX_INT32_H_INCLUDE_GUARD_
