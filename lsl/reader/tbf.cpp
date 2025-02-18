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
    uint8_t bytes[6144];       // The maximum possible size
} TBFPayload;


typedef struct __attribute__((packed)) {
    TBFHeader header;
    TBFPayload payload;
} TBFFrame;

/*
  State variables to keep track of the last read_tbf_impl used
*/

static int cached_nstand = 256;

PyObject *tbf_method = NULL;
PyObject *tbf_size = NULL;

template<int NSTAND, typename T, NPY_TYPES N>
PyObject *read_tbf_impl(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i, nstand;
    TBFFrame cFrame;
    
    static constexpr int frameSize = sizeof(TBFHeader)+8+12*NSTAND*2*1;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Read from the file - header + timestamp + NSTAND-sized buffer
    if( tbf_method == NULL ) {
        tbf_method = Py_BuildValue("s", "read");
    }
    if( tbf_size == NULL ) {
        tbf_size = Py_BuildValue("i", frameSize);
    }
    buffer = PyObject_CallMethodObjArgs(ph, tbf_method, tbf_size, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        return NULL;
    } else if( PyBytes_GET_SIZE(buffer) != frameSize ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(buffer);
        return NULL;
    }
    memcpy(&cFrame, PyBytes_AS_STRING(buffer), frameSize);
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
    
    // If nstand is not what we expect, update the cache and retry with correct size
    if( nstand != NSTAND ) {
        cached_nstand = nstand; // Update for next time
        Py_XDECREF(tbf_size);   // Update for next time
        tbf_size = NULL;        // Force NULL since Py_XDECREF isn't guaranteed to change tbf_size
        PyObject_CallMethod(ph, "seek", "ii", -frameSize, 1);
        
        switch(nstand) {
            case 64:  
                return read_tbf_impl<64, T, N>(self, args);
            case 256: 
                return read_tbf_impl<256, T, N>(self, args);
            default:
                PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
                return NULL;
        }
    }
    
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
    for(i=0; i<12*NSTAND*2*1; i+=4) {
        fp = tbfLUT[ cFrame.payload.bytes[i] ];
        *(a + 2*i + 0) = fp[0];
        *(a + 2*i + 1) = fp[1];
        fp = tbfLUT[ cFrame.payload.bytes[i+1] ];
        *(a + 2*i + 2) = fp[0];
        *(a + 2*i + 3) = fp[1];
        fp = tbfLUT[ cFrame.payload.bytes[i+2] ];
        *(a + 2*i + 4) = fp[0];
        *(a + 2*i + 5) = fp[1];
        fp = tbfLUT[ cFrame.payload.bytes[i+3] ];
        *(a + 2*i + 6) = fp[0];
        *(a + 2*i + 7) = fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.syncWord) ) {
        buffer = PyObject_CallMethod(ph, "seek", "ii", -frameSize, 1);
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to the frame object
    // 1.  Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = Py_BuildValue("B", cFrame.header.adp_id);
    PyObject_SetAttrString(fHeader, "adp_id", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.frame_count);
    PyObject_SetAttrString(fHeader, "frame_count", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.second_count);
    PyObject_SetAttrString(fHeader, "second_count", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("h", cFrame.header.first_chan);
    PyObject_SetAttrString(fHeader, "first_chan", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("h", nstand);
    PyObject_SetAttrString(fHeader, "nstand", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromLongLong(cFrame.payload.timetag);
    PyObject_SetAttrString(fPayload, "timetag", temp);
    Py_XDECREF(temp);
    
    PyObject_SetAttrString(fPayload, "_data", PyArray_Return(data));
    
    // 3. Frame
    PyObject_SetAttrString(frame, "header", fHeader);
    PyObject_SetAttrString(frame, "payload", fPayload);
    output = Py_BuildValue("O", frame);
    
    Py_XDECREF(fHeader);
    Py_XDECREF(fPayload);
    Py_XDECREF(data);
    
    return output;
}


template<typename T, NPY_TYPES N>
PyObject *read_tbf(PyObject *self, PyObject *args) {
    // Try fast path first with cached nstand
    switch(cached_nstand) {
        case 64:  
            return read_tbf_impl<64, T, N>(self, args);
        case 256: 
            return read_tbf_impl<256, T, N>(self, args);
        default:
            return read_tbf_impl<256, T, N>(self, args);
    }
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
