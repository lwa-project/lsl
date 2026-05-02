#include "Python.h"
#include <cmath>
#include <cstdint>
#include <complex>
#include <algorithm>
#include <fftw3.h>

#ifdef _OPENMP
    #include <omp.h>
    
    // OpenMP scheduling method
    #ifndef OMP_SCHEDULER
    #define OMP_SCHEDULER static
    #endif
#endif

#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "../correlator/common.hpp"
#include "../correlator/blas.hpp"

/* 
 Function to multiple two complex number together
*/

template<typename RealType>
inline std::complex<RealType> Complexm(std::complex<RealType> x, std::complex<RealType> y) {
    RealType re, im;
    re = (x.real() * y.real());
    im = (x.real() * y.imag());

    re = -(x.imag() * y.imag()) + re;
    im =  (x.imag() * y.real()) + im;

    return std::complex<RealType>(re, im);
}


/*
 Time block size control for the inner loops of the functions below
*/

#define TIME_BLOCK_DIM 32


/*
 Beamformer for complex coefficients and real data
*/

template<typename InType, typename DataType>
void phase_and_sum_real(long ninput,
                        long nchan,
                        long ntime,
                        InType const* coeffs,
                        DataType const* data,
                        unsigned char const* valid0,
                        unsigned char const* valid1,
                        Complex32* beam) {
  using RealType = typename InType::value_type;

  long norm0 = 0, norm1 = 0;
  for(long i=0; i<ninput; i++) {
      norm0 += (valid0[i] != 0);
      norm1 += (valid1[i] != 0);
  }
  RealType inv0 = norm0 > 0 ? (RealType) 1 / (RealType) norm0 : (RealType) 0;
  RealType inv1 = norm1 > 0 ? (RealType) 1 / (RealType) norm1 : (RealType) 0;

  Py_BEGIN_ALLOW_THREADS

  for(long c=0; c<nchan; c++) {
      Complex32* out0 = beam + 0*nchan*ntime + c*ntime;
      Complex32* out1 = beam + 1*nchan*ntime + c*ntime;

      InType* pol0 = (InType*) aligned64_malloc(TIME_BLOCK_DIM * sizeof(InType));
      InType* pol1 = (InType*) aligned64_malloc(TIME_BLOCK_DIM * sizeof(InType));

      for(long t=0; t<ntime; t+=TIME_BLOCK_DIM) {
          long bntime = std::min((long) TIME_BLOCK_DIM, ntime-t);

          ::memset((void*) pol0, 0, TIME_BLOCK_DIM * sizeof(InType));
          ::memset((void*) pol1, 0, TIME_BLOCK_DIM * sizeof(InType));

          for(long i=0; i<ninput; i++) {
              InType bc = coeffs[i*nchan + c];
              RealType v0 = (RealType) (valid0[i] != 0);
              RealType v1 = (RealType) (valid1[i] != 0);

              DataType const* row = data + i*nchan*ntime + c*ntime + t;
              for(long tb=0; tb<bntime; tb++) {
                  InType dp = InType((RealType) row[tb], (RealType) 0);
                  InType in = Complexm(bc, dp);
                  pol0[tb] += in*v0;
                  pol1[tb] += in*v1;
              }
          }

          for(long tb=0; tb<bntime; tb++) {
              out0[t+tb] = (Complex32) (pol0[tb] * inv0);
              out1[t+tb] = (Complex32) (pol1[tb] * inv1);
          }
      }

      aligned64_free(pol0);
      aligned64_free(pol1);
  }

  Py_END_ALLOW_THREADS
}


/*
 Beamformer for complex coefficients and complex data
*/

template<typename InType, typename DataType>
void phase_and_sum_complex(long ninput,
                           long nchan,
                           long ntime,
                           InType const* coeffs,
                           DataType const* data,
                           unsigned char const* valid0,
                           unsigned char const* valid1,
                           Complex32* beam) {
  using RealType = typename InType::value_type;

  long norm0 = 0, norm1 = 0;
  for(long i=0; i<ninput; i++) {
      norm0 += (valid0[i] != 0);
      norm1 += (valid1[i] != 0);
  }
  RealType inv0 = norm0 > 0 ? (RealType) 1 / (RealType) norm0 : (RealType) 0;
  RealType inv1 = norm1 > 0 ? (RealType) 1 / (RealType) norm1 : (RealType) 0;

  Py_BEGIN_ALLOW_THREADS

  for(long c=0; c<nchan; c++) {
      Complex32* out0 = beam + 0*nchan*ntime + c*ntime;
      Complex32* out1 = beam + 1*nchan*ntime + c*ntime;

      InType* pol0 = (InType*) aligned64_malloc(TIME_BLOCK_DIM * sizeof(InType));
      InType* pol1 = (InType*) aligned64_malloc(TIME_BLOCK_DIM * sizeof(InType));

      for(long t=0; t<ntime; t+=TIME_BLOCK_DIM) {
          long bntime = std::min((long) TIME_BLOCK_DIM, ntime-t);

          ::memset((void*) pol0, 0, TIME_BLOCK_DIM * sizeof(InType));
          ::memset((void*) pol1, 0, TIME_BLOCK_DIM * sizeof(InType));

          for(long i=0; i<ninput; i++) {
              InType bc = coeffs[i*nchan + c];
              RealType v0 = (RealType) (valid0[i] != 0);
              RealType v1 = (RealType) (valid1[i] != 0);

              DataType const* row = data + i*nchan*ntime*2 + c*ntime*2 + t*2;
              for(long tb=0; tb<bntime; tb++) {
                  InType dp = InType((RealType) row[2*tb], (RealType) row[2*tb+1]);
                  InType in = Complexm(bc, dp);
                  pol0[tb] += in*v0;
                  pol1[tb] += in*v1;
              }
          }
          
          for(long tb=0; tb<bntime; tb++) {
              out0[t+tb] = (Complex32) (pol0[tb] * inv0);
              out1[t+tb] = (Complex32) (pol1[tb] * inv1);
          }
      }
      
      aligned64_free(pol0);
      aligned64_free(pol1);
  }

  Py_END_ALLOW_THREADS
}


static PyObject *BEngine(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *coeffs, *signals, *validPol0, *validPol1, *beamsF;
    PyArrayObject *cdata=NULL, *data=NULL, *valid0=NULL, *valid1=NULL, *beams=NULL;
    int coeffType = NPY_COMPLEX64;

    long nStand, nChan, nTime;

    char const* kwlist[] = {"coeffs", "signals", "valid_pol0", "valid_pol1", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOOO", const_cast<char **>(kwlist), &coeffs, &signals, &validPol0, &validPol1)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }

    // Preserve the coeffs precision (complex64 or complex128) for both cdata and beams
    if( PyArray_Check(coeffs) && PyArray_TYPE((PyArrayObject *) coeffs) == NPY_COMPLEX128 ) {
        coeffType = NPY_COMPLEX128;
    }

    // Bring the data into C and make it usable
    cdata = (PyArrayObject *) PyArray_ContiguousFromObject(coeffs, coeffType, 2, 2);
    data = (PyArrayObject *) PyArray_ContiguousFromObject(signals,
                                                          PyArray_TYPE((PyArrayObject *) signals),
                                                          3, 4);
    valid0 = (PyArrayObject *) PyArray_ContiguousFromObject(validPol0, NPY_UINT8, 1, 1);
    valid1 = (PyArrayObject *) PyArray_ContiguousFromObject(validPol1, NPY_UINT8, 1, 1);
    if( cdata == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input coeffs array to 2-D complex");
        goto fail;
    }
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array as a 3-D array");
        goto fail;
    }
    if( valid0 == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input valid_pol0 array to 1-D uint8");
        goto fail;
    }
    if( valid1 == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input valid_pol1 array to 1-D uint8");
        goto fail;
    }
    
    // Check data dimensions
    if(PyArray_DIM(cdata, 0) != PyArray_DIM(data, 0)) {
        PyErr_Format(PyExc_RuntimeError, "coeff and signals have different stand counts");
        goto fail;
    }
    if(PyArray_DIM(cdata, 0) != PyArray_DIM(valid0, 0)) {
        PyErr_Format(PyExc_RuntimeError, "coeff and valid_pol0 have different stand counts");
        goto fail;
    }
    if(PyArray_DIM(cdata, 0) != PyArray_DIM(valid1, 0)) {
        PyErr_Format(PyExc_RuntimeError, "coeff and valid_pol1 have different stand counts");
        goto fail;
    }
    
    if(PyArray_DIM(cdata, 1) != PyArray_DIM(data, 1)) {
        PyErr_Format(PyExc_RuntimeError, "coeff and signals have different channel counts");
        goto fail;
    }
    
    // Get the properties of the data
    nStand = PyArray_DIM(cdata, 0);
    nChan = PyArray_DIM(cdata, 1);
    nTime = PyArray_DIM(data, 2);
    
    // Find out how large the output array needs to be and initialize it
    npy_intp dims[3];
    dims[0] = (npy_intp) 2;
    dims[1] = (npy_intp) nChan;
    dims[2] = (npy_intp) nTime;
    beams = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
    if(beams == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }

#define LAUNCH_PHASE_AND_SUM_R(CoeffType, IterType) \
        phase_and_sum_real<CoeffType,IterType>(nStand, nChan, nTime,\
                                        (CoeffType*) PyArray_DATA(cdata), \
                                        (IterType*) PyArray_DATA(data), \
                                        (unsigned char*) PyArray_DATA(valid0), \
                                        (unsigned char*) PyArray_DATA(valid1),\
                                        (Complex32*) PyArray_DATA(beams))

#define LAUNCH_PHASE_AND_SUM_C(CoeffType, IterType) \
        phase_and_sum_complex<CoeffType,IterType>(nStand, nChan, nTime,\
                                        (CoeffType*) PyArray_DATA(cdata), \
                                        (IterType*) PyArray_DATA(data), \
                                        (unsigned char*) PyArray_DATA(valid0), \
                                        (unsigned char*) PyArray_DATA(valid1),\
                                        (Complex32*) PyArray_DATA(beams))

    if( coeffType == NPY_COMPLEX64 ) {
        switch( PyArray_LSL_TYPE(data, 3) ){
            case( NPY_INT8       ): LAUNCH_PHASE_AND_SUM_R(Complex32, int8_t); break;
            case( NPY_INT16      ): LAUNCH_PHASE_AND_SUM_R(Complex32, int16_t); break;
            case( NPY_INT32      ): LAUNCH_PHASE_AND_SUM_R(Complex32, int32_t); break;
            case( NPY_INT64      ): LAUNCH_PHASE_AND_SUM_R(Complex32, int64_t); break;
            case( NPY_FLOAT32    ): LAUNCH_PHASE_AND_SUM_R(Complex32, float); break;
            case( NPY_FLOAT64    ): LAUNCH_PHASE_AND_SUM_R(Complex32, double); break;
            case( LSL_CI8        ): LAUNCH_PHASE_AND_SUM_C(Complex32, int8_t); break;
            case( NPY_COMPLEX64  ): LAUNCH_PHASE_AND_SUM_C(Complex32, float); break;
            case( NPY_COMPLEX128 ): LAUNCH_PHASE_AND_SUM_C(Complex32, double); break;
            default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
        }
    } else if( coeffType == NPY_COMPLEX128 )  {
        switch( PyArray_LSL_TYPE(data, 3) ){
            case( NPY_INT8       ): LAUNCH_PHASE_AND_SUM_R(Complex64, int8_t); break;
            case( NPY_INT16      ): LAUNCH_PHASE_AND_SUM_R(Complex64, int16_t); break;
            case( NPY_INT32      ): LAUNCH_PHASE_AND_SUM_R(Complex64, int32_t); break;
            case( NPY_INT64      ): LAUNCH_PHASE_AND_SUM_R(Complex64, int64_t); break;
            case( NPY_FLOAT32    ): LAUNCH_PHASE_AND_SUM_R(Complex64, float); break;
            case( NPY_FLOAT64    ): LAUNCH_PHASE_AND_SUM_R(Complex64, double); break;
            case( LSL_CI8        ): LAUNCH_PHASE_AND_SUM_C(Complex64, int8_t); break;
            case( NPY_COMPLEX64  ): LAUNCH_PHASE_AND_SUM_C(Complex64, float); break;
            case( NPY_COMPLEX128 ): LAUNCH_PHASE_AND_SUM_C(Complex64, double); break;
            default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
        }
    } else {
        PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }

#undef LAUNCH_PHASE_AND_SUM_R
#undef LAUNCH_PHASE_AND_SUM_C
    
    beamsF = Py_BuildValue("O", PyArray_Return(beams));
    Py_XDECREF(cdata);
    Py_XDECREF(data);
    Py_XDECREF(valid0);
    Py_XDECREF(valid1);
    Py_XDECREF(beams);
    
    return beamsF;
    
fail:
    Py_XDECREF(cdata);
    Py_XDECREF(data);
    Py_XDECREF(valid0);
    Py_XDECREF(valid1);
    Py_XDECREF(beams);
    
    return NULL;
}

PyDoc_STRVAR(BEngine_doc, \
"Perform polarized phase-and-sum beamforming of complex data\n\
\n\
Input arguments are:\n\
 * coeffs: 2-D numpy.complex64/128 (stands by channels) array of beamforming\n\
   coefficients\n\
 * signals: 3-D numpy (stands by channels by time) array of data\n\
 * valid_pol0: 1-D numpy.uint8 array of which stands contribute to the beam's\n\
   first polarization\n\
 * valid_pol2: 1-D numpy.uint8 array of which stands contribute to the beam's\n\
   second polarization\n\
\n\
Outputs:\n\
 * beam: 3-D numpy.complex64 (pol by channels by time) of beamformed data\n\
");

/*
Module Setup - Function Definitions and Documentation
*/

static PyMethodDef beamformer_methods[] = {
    {"BEngine",   (PyCFunction) BEngine,   METH_VARARGS|METH_KEYWORDS, BEngine_doc  },
    {NULL,        NULL,                    0,                          NULL         }
};

PyDoc_STRVAR(beamformer_doc, \
"C-based phase-and-sum beamforming for complex data.\n\
\n\
The function defined in this module are:\n\
 * BEngine - Main (only) function in this module.\n\
\n\
See the individual functions for more details.");


/*
Module Setup - Initialization
*/

static int beamformer_exec(PyObject *module) {
    import_array1(-1);
    
    // Version and revision information
    PyModule_AddObject(module, "__version__", PyUnicode_FromString("0.1"));
    
    // Function listings
    PyObject* all = PyList_New(0);
    PyList_Append(all, PyUnicode_FromString("BEngine"));
    PyModule_AddObject(module, "__all__", all);
    return 0;
}

static PyModuleDef_Slot beamformer_slots[] = {
    {Py_mod_exec, (void *)&beamformer_exec},
    {0,           NULL}
};

static PyModuleDef beamformer_def = {
    PyModuleDef_HEAD_INIT,    /* m_base */
    "_beamformer",            /* m_name */
    beamformer_doc,           /* m_doc */
    0,                        /* m_size */
    beamformer_methods,       /* m_methods */
    beamformer_slots,         /* m_slots */
    NULL,                     /* m_traverse */
    NULL,                     /* m_clear */
    NULL,                     /* m_free */
};

PyMODINIT_FUNC PyInit__beamformer(void) {
    return PyModuleDef_Init(&beamformer_def);
}
