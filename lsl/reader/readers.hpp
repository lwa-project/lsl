#pragma once

#include <cstdint>

#include "Python.h"

#if defined(__linux__)
/* Linux */
#include <byteswap.h>
#elif defined(__APPLE__) && defined(__MACH__)
/* OSX */
#include <libkern/OSByteOrder.h>
#define __bswap_16 OSSwapInt16
#define __bswap_32 OSSwapInt32
#define __bswap_64 OSSwapInt64
#endif


/*
  Minimum function for two values
*/

#define min(a,b) (((a) < (b)) ? (a) : (b))


/*
  Exceptions for the Go Fast! Readers
*/

// gofast.cpp
extern PyObject *SyncError;
extern PyObject *EOFError;


/*
  LWA1/LWA-SV Look up tables
*/

// gofast.cpp
extern int16_t tbw4LUT[256][2];
extern int8_t  tbnLUT[256];
extern int8_t  drx8LUT[256];
extern int8_t  drxLUT[256][2];
extern int8_t  tbfLUT[256][2];


/* 
  Support Functions
*/

// Validate a collection of Mark 5C sync words.  Return true if all are valid.
inline bool validSync5C(uint32_t syncWord) {
    if( syncWord == 0x5CDEC0DE ) {
        return true;
    }
    return false;
}

  

// vdif.cpp
extern void initVDIFLUTs(void);

/*
 ADP COR mode channel information
*/

#define COR_NCHAN 72

/*
  Reader Functions and Documentation
*/

// tbw.cpp
extern PyObject *read_tbw(PyObject*, PyObject*);
extern char read_tbw_doc[];
// tbn.cpp
extern PyObject *read_tbn_cf32(PyObject*, PyObject*);
extern char read_tbn_cf32_doc[];
extern PyObject *read_tbn_ci8(PyObject*, PyObject*);
extern char read_tbn_ci8_doc[];
// drx.cpp
extern PyObject *read_drx_cf32(PyObject*, PyObject*);
extern char read_drx_cf32_doc[];
extern PyObject *read_drx_ci8(PyObject*, PyObject*);
extern char read_drx_ci8_doc[];
// drx8.cpp
extern PyObject *read_drx8_cf32(PyObject*, PyObject*);
extern char read_drx8_cf32_doc[];
extern PyObject *read_drx8_ci8(PyObject*, PyObject*);
extern char read_drx8_ci8_doc[];
// drspec.cpp
extern PyObject *read_drspec(PyObject*, PyObject*);
extern char read_drspec_doc[];

// vdif.cpp
extern PyObject *read_vdif_f32(PyObject*, PyObject*, PyObject*);
extern char read_vdif_f32_doc[];
extern PyObject *read_vdif_i8(PyObject*, PyObject*, PyObject*);
extern char read_vdif_i8_doc[];

// tbf.cpp
extern PyObject *read_tbf_cf32(PyObject*, PyObject*);
extern char read_tbf_cf32_doc[];
extern PyObject *read_tbf_ci8(PyObject*, PyObject*);
extern char read_tbf_ci8_doc[];
// cor.cpp
extern PyObject *read_cor(PyObject*, PyObject*);
extern char read_cor_doc[];
