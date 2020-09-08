#ifndef COMPLEX_COMPLEX_INT8_H_INCLUDE_GUARD_
#define COMPLEX_COMPLEX_INT8_H_INCLUDE_GUARD_

#ifdef __cplusplus
extern "C" {
#endif

#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

typedef struct {
	unsigned char real_imag;
} complex_int8;

// Unpacking/repacking functions
extern signed char fourBitLUT[256][2];
void complex_int8_fillLUT();

static NPY_INLINE complex_int8 pack_ci8(signed char real, signed char imag) {
    complex_int8 c;
    c.real_imag  = (real * 16) & 0xF0;
    c.real_imag |= ((imag * 16) >> 4) & 0x0F;
    return c;
}

static NPY_INLINE void inplace_pack_ci8(signed char real, signed char imag, complex_int8 *c) {
    c->real_imag  = (real * 16) & 0xF0;
    c->real_imag |= ((imag * 16) >> 4) & 0x0F;
}

// Unary bool operators
static NPY_INLINE int complex_int8_nonzero(complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    return sc[0] != 0 || sc[1] != 0;
}

static NPY_INLINE int complex_int8_isnan(complex_int8 c) {
    return 0;
}

static NPY_INLINE int complex_int8_isinf(complex_int8 c) {
    return 0;
}

static NPY_INLINE int complex_int8_isfinite(complex_int8 c) {
    return 1;
}

// Binary bool operators
static NPY_INLINE int complex_int8_equal(complex_int8 c1, complex_int8 c2) {
    return 
        !complex_int8_isnan(c1) &&
        !complex_int8_isnan(c2) &&
        c1.real_imag == c2.real_imag; 
}

static NPY_INLINE int complex_int8_not_equal(complex_int8 c1, complex_int8 c2) {
    return !complex_int8_equal(c1, c2);
}

static NPY_INLINE int complex_int8_less(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    return
        (!complex_int8_isnan(c1) &&
         !complex_int8_isnan(c2)) && (
            sc1[0] != sc2[0] ? sc1[0] < sc2[0] :
            sc1[1] != sc2[1] ? sc1[1] < sc2[1] : 0);
}

static NPY_INLINE int complex_int8_greater(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    return
        (!complex_int8_isnan(c1) &&
         !complex_int8_isnan(c2)) && (
            sc1[0] != sc2[0] ? sc1[0] > sc2[0] :
            sc1[1] != sc2[1] ? sc1[1] > sc2[1] : 0);
}

static NPY_INLINE int complex_int8_less_equal(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    return
        (!complex_int8_isnan(c1) &&
         !complex_int8_isnan(c2)) && (
            sc1[0] != sc2[0] ? sc1[0] < sc2[0] :
            sc1[1] != sc2[1] ? sc1[1] < sc2[1] : 1);
}

static NPY_INLINE int complex_int8_greater_equal(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    return
        (!complex_int8_isnan(c1) &&
         !complex_int8_isnan(c2)) && (
            sc1[0] != sc2[0] ? sc1[0] > sc2[0] :
            sc1[1] != sc2[1] ? sc1[1] > sc2[1] : 1);
}

// Unary int returners
static NPY_INLINE long complex_int8_real(complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    return (long) sc[0];
}

static NPY_INLINE long complex_int8_imag(complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    return (long) sc[1];
}

// Unary float returners
static NPY_INLINE double complex_int8_absolute(complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    return sqrt((sc[0]*1.0)*sc[0] + (sc[1]*1.0)*sc[1]);
}

static NPY_INLINE double complex_int8_kludgy_arctan2(long i, complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    return atan2((double) sc[1], (double) sc[0]);
}

// Unary complex_int8 returners
static NPY_INLINE complex_int8 complex_int8_negative(complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real = -sc[0];
    signed char imag = -sc[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE complex_int8 complex_int8_conjugate(complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real =  sc[0];
    signed char imag = -sc[1];
    return pack_ci8(real, imag);
}

// complex_int8-complex_int8/complex_int8-scalar/scalar-complex_int8 binary complex_int8 returners
static NPY_INLINE complex_int8 complex_int8_add(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    signed char real = sc1[0] + sc2[0];
    signed char imag = sc1[1] + sc2[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_add(complex_int8* c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1->real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    signed char real = sc1[0] + sc2[0];
    signed char imag = sc1[1] + sc2[1];
    inplace_pack_ci8(real, imag, c1);
}

static NPY_INLINE complex_int8 complex_int8_scalar_add(long s, complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real = s + sc[0];
    signed char imag = 0 + sc[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_scalar_add(long s, complex_int8* c) {
    const signed char* sc = fourBitLUT[c->real_imag];
    signed char real = s + sc[0];
    signed char imag = 0 + sc[1];
    inplace_pack_ci8(real, imag, c);
}

static NPY_INLINE complex_int8 complex_int8_add_scalar(complex_int8 c, long s) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real = s + sc[0];
    signed char imag = 0 + sc[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_add_scalar(complex_int8* c, long s) {
    const signed char* sc = fourBitLUT[c->real_imag];
    signed char real = s + sc[0];
    signed char imag = 0 + sc[1];
    inplace_pack_ci8(real, imag, c);
}

static NPY_INLINE complex_int8 complex_int8_subtract(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    signed char real = sc1[0] - sc2[0];
    signed char imag = sc1[1] - sc2[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_subtract(complex_int8* c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1->real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    signed char real = sc1[0] - sc2[0];
    signed char imag = sc1[1] - sc2[1];
    inplace_pack_ci8(real, imag, c1);
}

static NPY_INLINE complex_int8 complex_int8_scalar_subtract(long s, complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real = s - sc[0];
    signed char imag = 0 - sc[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_scalar_subtract(long s, complex_int8* c) {
    const signed char* sc = fourBitLUT[c->real_imag];
    signed char real = s - sc[0];
    signed char imag = 0 - sc[1];
    inplace_pack_ci8(real, imag, c);
}

static NPY_INLINE complex_int8 complex_int8_subtract_scalar(complex_int8 c, long s) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real = -s + sc[0];
    signed char imag = -0 + sc[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_subtract_scalar(complex_int8* c, long s) {
    const signed char* sc = fourBitLUT[c->real_imag];
    signed char real = -s + sc[0];
    signed char imag = -0 + sc[1];
    inplace_pack_ci8(real, imag, c);
}

static NPY_INLINE complex_int8 complex_int8_multiply(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    signed char real = sc1[0]*sc2[0] - sc1[1]*sc2[1];
    signed char imag = sc1[1]*sc2[0] + sc1[0]*sc2[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_multiply(complex_int8* c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1->real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    signed char real = sc1[0]*sc2[0] - sc1[1]*sc2[1];
    signed char imag = sc1[1]*sc2[0] + sc1[0]*sc2[1];
    inplace_pack_ci8(real, imag, c1);
}

static NPY_INLINE complex_int8 complex_int8_scalar_multiply(long s, complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real = s*sc[0];
    signed char imag = s*sc[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_scalar_multiply(long s, complex_int8* c) {
    const signed char* sc = fourBitLUT[c->real_imag];
    signed char real = s*sc[0];
    signed char imag = s*sc[1];
    inplace_pack_ci8(real, imag, c);
}

static NPY_INLINE complex_int8 complex_int8_multiply_scalar(complex_int8 c, long s) {
    const signed char* sc = fourBitLUT[c.real_imag];
    signed char real = s*sc[0];
    signed char imag = s*sc[1];
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_multiply_scalar(complex_int8* c, long s) {
    const signed char* sc = fourBitLUT[c->real_imag];
    signed char real = s*sc[0];
    signed char imag = s*sc[1];
    inplace_pack_ci8(real, imag, c);
}

static NPY_INLINE complex_int8 complex_int8_divide(complex_int8 c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1.real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    long mag2 = ((int) sc2[0])*sc2[0] + ((int) sc2[1])*sc2[1];
    signed char real = (sc1[0]*sc2[0] + sc1[1]*sc2[1]) / mag2;
    signed char imag = (sc1[1]*sc2[0] - sc1[0]*sc2[1]) / mag2;
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_divide(complex_int8* c1, complex_int8 c2) {
    const signed char* sc1 = fourBitLUT[c1->real_imag];
    const signed char* sc2 = fourBitLUT[c2.real_imag];
    long mag2 = ((int) sc2[0])*sc2[0] + ((int) sc2[1])*sc2[1];
    signed char real = (sc1[0]*sc2[0] + sc1[1]*sc2[1]) / mag2;
    signed char imag = (sc1[1]*sc2[0] - sc1[0]*sc2[1]) / mag2;
    inplace_pack_ci8(real, imag, c1);
}

static NPY_INLINE complex_int8 complex_int8_scalar_divide(long s, complex_int8 c) {
    const signed char* sc = fourBitLUT[c.real_imag];
    long mag2 = ((int) sc[0])*sc[0] + ((int) sc[1])*sc[1];
    signed char real = (s*sc[0] + 0*sc[1]) / mag2;
    signed char imag = (0*sc[0] - s*sc[1]) / mag2;
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_scalar_divide(long s, complex_int8* c) {
    const signed char* sc = fourBitLUT[c->real_imag];
    long mag2 = ((int) sc[0])*sc[0] + ((int) sc[1])*sc[1];
    signed char real = (s*sc[0] + 0*sc[1]) / mag2;
    signed char imag = (0*sc[0] - s*sc[1]) / mag2;
    inplace_pack_ci8(real, imag, c);
}

static NPY_INLINE complex_int8 complex_int8_divide_scalar(complex_int8 c, long s) {
    const signed char* sc = fourBitLUT[c.real_imag];
    long mag2 = s*s + 0*0;
    signed char real = (sc[0]*s - sc[1]*0) / mag2;
    signed char imag = (sc[1]*s + sc[0]*0) / mag2;
    return pack_ci8(real, imag);
}

static NPY_INLINE void complex_int8_inplace_divide_scalar(complex_int8* c, long s) {
    const signed char* sc = fourBitLUT[c->real_imag];
    long mag2 = s*s + 0*0;
    signed char real = (sc[0]*s - sc[1]*0) / mag2;
    signed char imag = (sc[1]*s + sc[0]*0) / mag2;
    inplace_pack_ci8(real, imag, c);
}

#ifdef __cplusplus
}
#endif

#endif // COMPLEX_COMPLEX_INT8_H_INCLUDE_GUARD_
