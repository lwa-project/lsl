#ifndef COMPLEX_COMPLEX_INT16_H_INCLUDE_GUARD_
#define COMPLEX_COMPLEX_INT16_H_INCLUDE_GUARD_

#ifdef __cplusplus
extern "C" {
#endif

#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

typedef struct {
    signed char real;
    signed char imag;
} complex_int16;

// Access methods
static NPY_INLINE void lsl_unpack_ci16(complex_int16 packed, signed char* real, signed char* imag) {
    *real = packed.real;
    *imag = packed.imag;
}

static NPY_INLINE void lsl_pack_ci16(signed char real, signed char imag, complex_int16* packed) {
    packed->real = real;
    packed->imag = imag;
}

// Unary bool operators
static NPY_INLINE int complex_int16_nonzero(complex_int16 c) {
    return c.real != 0 || c.imag != 0;
}

static NPY_INLINE int complex_int16_isnan(complex_int16 c) {
    return 0;
}

static NPY_INLINE int complex_int16_isinf(complex_int16 c) {
    return 0;
}

static NPY_INLINE int complex_int16_isfinite(complex_int16 c) {
    return 1;
}

// Binary bool operators
static NPY_INLINE int complex_int16_equal(complex_int16 c1, complex_int16 c2) {
    return 
        !complex_int16_isnan(c1) &&
        !complex_int16_isnan(c2) &&
        c1.real == c2.real &&
        c1.imag == c2.imag; 
}

static NPY_INLINE int complex_int16_not_equal(complex_int16 c1, complex_int16 c2) {
    return !complex_int16_equal(c1, c2);
}

static NPY_INLINE int complex_int16_less(complex_int16 c1, complex_int16 c2) {
    return
        (!complex_int16_isnan(c1) &&
         !complex_int16_isnan(c2)) && (
            c1.real != c2.real ? c1.real < c2.real :
            c1.imag != c2.imag ? c1.imag < c2.imag : 0);
}

static NPY_INLINE int complex_int16_greater(complex_int16 c1, complex_int16 c2) {
    return
        (!complex_int16_isnan(c1) &&
         !complex_int16_isnan(c2)) && (
            c1.real != c2.real ? c1.real > c2.real :
            c1.imag != c2.imag ? c1.imag > c2.imag : 0);
}

static NPY_INLINE int complex_int16_less_equal(complex_int16 c1, complex_int16 c2) {
    return
        (!complex_int16_isnan(c1) &&
         !complex_int16_isnan(c2)) && (
            c1.real != c2.real ? c1.real < c2.real :
            c1.imag != c2.imag ? c1.imag < c2.imag : 1);
}

static NPY_INLINE int complex_int16_greater_equal(complex_int16 c1, complex_int16 c2) {
    return
        (!complex_int16_isnan(c1) &&
         !complex_int16_isnan(c2)) && (
            c1.real != c2.real ? c1.real > c2.real :
            c1.imag != c2.imag ? c1.imag > c2.imag : 1);
}

// Unary int returners
static NPY_INLINE long complex_int16_real(complex_int16 c) {
    return (long) c.real;
}

static NPY_INLINE long complex_int16_imag(complex_int16 c) {
    return (long) c.imag;
}

// Unary float returners
static NPY_INLINE double complex_int16_absolute(complex_int16 c) {
    return sqrt(((int) c.real)*c.real + ((int) c.imag)*c.imag);
}

static NPY_INLINE double complex_int16_kludgy_arctan2(long i, complex_int16 c) {
    return atan2((int) c.imag, (int) c.real);
}

// Unary complex_int16 returners
static NPY_INLINE complex_int16 complex_int16_negative(complex_int16 c) {
    signed char real = -c.real;
    signed char imag = -c.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE complex_int16 complex_int16_conjugate(complex_int16 c) {
    signed char real =  c.real;
    signed char imag = -c.imag;
    return (complex_int16) {real, imag};
}

// complex_int16-complex_int16/complex_int16-scalar/scalar-complex_int16 binary complex_int16 returners
static NPY_INLINE complex_int16 complex_int16_add(complex_int16 c1, complex_int16 c2) {
    signed char real = c1.real + c2.real;
    signed char imag = c1.imag + c2.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_add(complex_int16* c1, complex_int16 c2) {
    signed char real = c1->real + c2.real;
    signed char imag = c1->imag + c2.imag;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_scalar_add(long s, complex_int16 c) {
    signed char real = s + c.real;
    signed char imag = 0 + c.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_scalar_add(long s, complex_int16* c) {
    signed char real = s + c->real;
    signed char imag = 0 + c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_add_scalar(complex_int16 c, long s) {
    signed char real = s + c.real;
    signed char imag = 0 + c.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_add_scalar(complex_int16* c, long s) {
    signed char real = s + c->real;
    signed char imag = 0 + c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_subtract(complex_int16 c1, complex_int16 c2) {
    signed char real = c1.real - c2.real;
    signed char imag = c1.imag - c2.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_subtract(complex_int16* c1, complex_int16 c2) {
    signed char real = c1->real - c2.real;
    signed char imag = c1->imag - c2.imag;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_scalar_subtract(long s, complex_int16 c) {
    signed char real = s - c.real;
    signed char imag = 0 - c.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_scalar_subtract(long s, complex_int16* c) {
    signed char real = s - c->real;
    signed char imag = 0 - c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_subtract_scalar(complex_int16 c, long s) {
    signed char real = -s + c.real;
    signed char imag = -0 + c.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_subtract_scalar(complex_int16* c, long s) {
    signed char real = -s + c->real;
    signed char imag = -0 + c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_multiply(complex_int16 c1, complex_int16 c2) {
    signed char real = c1.real*c2.real - c1.imag*c2.imag;
    signed char imag = c1.imag*c2.real + c1.real*c2.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_multiply(complex_int16* c1, complex_int16 c2) {
    signed char real = c1->real*c2.real - c1->imag*c2.imag;
    signed char imag = c1->imag*c2.real + c1->real*c2.imag;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_scalar_multiply(long s, complex_int16 c) {
    signed char real = s*c.real;
    signed char imag = s*c.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_scalar_multiply(long s, complex_int16* c) {
    signed char real = s*c->real;
    signed char imag = s*c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_multiply_scalar(complex_int16 c, long s) {
    signed char real = s*c.real;
    signed char imag = s*c.imag;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_multiply_scalar(complex_int16* c, long s) {
    signed char real = s*c->real;
    signed char imag = s*c->imag;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_divide(complex_int16 c1, complex_int16 c2) {
    long mag2 = ((int) c2.real)*c2.real + ((int) c2.imag)*c2.imag;
    signed char real = (c1.real*c2.real + c1.imag*c2.imag) / mag2;
    signed char imag = (c1.imag*c2.real - c1.real*c2.imag) / mag2;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_divide(complex_int16* c1, complex_int16 c2) {
    long mag2 = ((int) c2.real)*c2.real + ((int) c2.imag)*c2.imag;
    signed char real = (c1->real*c2.real + c1->imag*c2.imag) / mag2;
    signed char imag = (c1->imag*c2.real - c1->real*c2.imag) / mag2;
    c1->real = real;
    c1->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_scalar_divide(long s, complex_int16 c) {
    long mag2 = ((int) c.real)*c.real + ((int) c.imag)*c.imag;
    signed char real = (s*c.real + 0*c.imag) / mag2;
    signed char imag = (0*c.real - s*c.imag) / mag2;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_scalar_divide(long s, complex_int16* c) {
    long mag2 = ((int) c->real)*c->real + ((int) c->imag)*c->imag;
    signed char real = (s*c->real + 0*c->imag) / mag2;
    signed char imag = (0*c->real - s*c->imag) / mag2;
    c->real = real;
    c->imag = imag;
}

static NPY_INLINE complex_int16 complex_int16_divide_scalar(complex_int16 c, long s) {
    long mag2 = s*s + 0*0;
    signed char real = (c.real*s - c.imag*0) / mag2;
    signed char imag = (c.imag*s + c.real*0) / mag2;
    return (complex_int16) {real, imag};
}

static NPY_INLINE void complex_int16_inplace_divide_scalar(complex_int16* c, long s) {
    long mag2 = s*s + 0*0;
    signed char real = (c->real*s - c->imag*0) / mag2;
    signed char imag = (c->imag*s + c->real*0) / mag2;
    c->real = real;
    c->imag = imag;
}

#ifdef __cplusplus
}
#endif

#endif // COMPLEX_COMPLEX_INT16_H_INCLUDE_GUARD_
