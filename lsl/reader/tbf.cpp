#include "Python.h"
#include <cmath>
#include <cstdint>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.hpp"

/*
  TBF Reader
*/

typedef struct __attribute__((packed)) {
    uint32_t syncWord;
    union {
        struct {
            uint32_t frame_count:24;
            uint8_t  adp_id:8;
        };
        uint32_t frame_count_word;
    };
    uint32_t second_count;
    int16_t  first_chan;
    int16_t  nstand;
} TBFHeader;


typedef struct __attribute__((packed)) {
    int64_t timetag;
} TBFPayload;


typedef struct __attribute__((packed)) {
    TBFHeader header;
    TBFPayload payload;
} TBFFrame;


PyObject *tbf_method = NULL;
PyObject *tbf_size   = NULL;
PyObject *tbf_dsize  = NULL;


template<typename T, NPY_TYPES N>
PyObject *read_tbf(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output;
    PyArrayObject *data;
    int i, nstand;
    const uint8_t *raw_data;
    TBFFrame cFrame;
    
    if(!PyArg_ParseTuple(args, "O", &ph)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Read from the file - header + timestamp
    if( tbf_method == NULL ) {
        tbf_method = Py_BuildValue("s", "read");
        tbf_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, tbf_method, tbf_size, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        return NULL;
    } else if( PyBytes_GET_SIZE(buffer) != sizeof(cFrame) ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(buffer);
        return NULL;
    }
    memcpy(&cFrame, PyBytes_AS_STRING(buffer), sizeof(cFrame));
    Py_XDECREF(buffer);
    
    // Determine the number of stands in the frame - default to 256 for
    // legacy LWA-SV data
    nstand = 256;
    if( cFrame.header.nstand != 0 ) {
      nstand = __bswap_16(cFrame.header.nstand);
      if( nstand > 256 || nstand < 0 ) {
        cFrame.header.syncWord = 0;
        nstand = 256;
      }
    }
    tbf_dsize = Py_BuildValue("i", 12*nstand*2*1);
    
    // Read from the file - data
    buffer = PyObject_CallMethodObjArgs(ph, tbf_method, tbf_dsize, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        return NULL;
    } else if( PyBytes_GET_SIZE(buffer) != 12*nstand*2*1 ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(buffer);
        return NULL;
    }
    Py_buffer view;
    if (PyObject_GetBuffer(buffer, &view, PyBUF_SIMPLE) == -1) {
        PyErr_Format(PyExc_RuntimeError, "Cannot create buffer from read return value");
        Py_XDECREF(data);
        Py_XDECREF(buffer);
        return NULL;
    }
    raw_data = (const uint8_t*) view.buf;
    
    // Create the output data array
    npy_intp dims[3];
    // 4+4-bit Data -> 6144 samples in the data section as 12 channels, 256 stands, and 2 pols.
    dims[0] = (npy_intp) 12;
    dims[1] = (npy_intp) nstand;
    dims[2] = (npy_intp) 2;
    if(N == NPY_INT8) dims[2] *= 2;
    data = (PyArrayObject*) PyArray_ZEROS(3, dims, N, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    // Swap the bits around
    cFrame.header.frame_count_word = __bswap_32(cFrame.header.frame_count_word);
    cFrame.header.second_count = __bswap_32(cFrame.header.second_count);
    cFrame.header.first_chan = __bswap_16(cFrame.header.first_chan);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    
    // Fill the data array
    const int8_t *fp;
    T *a;
    a = (T *) PyArray_DATA(data);
    for(i=0; i<12*nstand*2*1; ) {   // Sets of 12 channels for nstand stands @ 4+4 bit
        fp = tbfLUT[ raw_data[i++] ];
        *a++ = fp[0]; *a++ = fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    PyBuffer_Release(&view);
    Py_XDECREF(buffer);
    
    // Validate
    if( !validSync5C(cFrame.header.syncWord) ) {
        buffer = PyObject_CallMethod(ph, "seek", "ii", -sizeof(cFrame)-12*nstand*2*1, 1);
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to a tuple
    output = Py_BuildValue("(BIIHHKO)", cFrame.header.adp_id, cFrame.header.frame_count, cFrame.header.second_count, cFrame.header.first_chan, nstand, cFrame.payload.timetag, PyArray_Return(data));
    
    Py_XDECREF(data);
    
    return output;
}


PyObject *read_tbf_cf32(PyObject *self, PyObject *args) {
    return read_tbf<float,NPY_COMPLEX64>(self, args);
}

char read_tbf_cf32_doc[] = PyDoc_STR(\
"Function to read in a single TBW frame (header+data) and store the contents\n\
as a Frame object.\n\
\n\
.. versionadded:: 1.2.0\n\
");


PyObject *read_tbf_ci8(PyObject *self, PyObject *args) {
    return read_tbf<int8_t,NPY_INT8>(self, args);
}

char read_tbf_ci8_doc[] = PyDoc_STR(\
"Function to read in a single TBW frame (header+data) and store the contents\n\
as a Frame object.\n\
\n\
.. note::\n\
    This function differs from `read_tbf` in that it returns a\n\
    `lsl.tbf.FramePayload` containing a 4-D numpy.int8 array (channels by\n\
    stands by polarizations by real/complex) rather than a 1-D numpy.complex64\n\
    array.\n\
\n\
.. versionadded:: 2.1.3\n\
");
