#include "Python.h"
#include <cmath>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.hpp"

/*
  DRX Reader
*/

typedef struct __attribute__((packed)) {
    uint32_t sync_word;
    union {
        struct {
            uint32_t frame_count:24;
            uint8_t id:8;
        };
        /* Also: 
            struct {
                uint32_t frame_count:24;
                uint8_t  beam:3;
                uint8_t  tune:3;
                uint8_t  reserved:1;
                uint8_t  pol:1;
            };
        */
        uint32_t frame_count_word;
    };
    uint32_t second_count;
    uint16_t decimation;
    uint16_t time_offset;
} DRXHeader;


typedef struct __attribute__((packed)) {
    uint64_t timetag;
    uint32_t tuning_word;
    uint32_t flags;
    uint8_t  bytes[4096];
} DRXPayload;


typedef struct __attribute__((packed)) {
    DRXHeader header;
    DRXPayload payload;
} DRXFrame;


PyObject *drx_method = NULL;
PyObject *drx_size   = NULL;


template<typename T, NPY_TYPES N>
PyObject *read_drx(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i;
    DRXFrame cFrame;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[1];
    dims[0] = (npy_intp) 4096;
    if( N == NPY_INT8) dims[0] *= 2;
    data = (PyArrayObject*) PyArray_ZEROS(1, dims, N, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Read from the file
    if( drx_method == NULL ) {
        drx_method = Py_BuildValue("s", "read");
        drx_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, drx_method, drx_size, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        Py_XDECREF(data);
        return NULL;
    } else if( PyBytes_GET_SIZE(buffer) != sizeof(cFrame) ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(data);
        Py_XDECREF(buffer);
        return NULL;
    }
    memcpy(&cFrame, PyBytes_AS_STRING(buffer), sizeof(cFrame));
    Py_XDECREF(buffer);
        
    Py_BEGIN_ALLOW_THREADS
    
    // Swap the bits around
    cFrame.header.frame_count_word = __bswap_32(cFrame.header.frame_count_word);
    cFrame.header.second_count = __bswap_32(cFrame.header.second_count);
    cFrame.header.decimation = __bswap_16(cFrame.header.decimation);
    cFrame.header.time_offset = __bswap_16(cFrame.header.time_offset);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    cFrame.payload.tuning_word = __bswap_32(cFrame.payload.tuning_word);
    cFrame.payload.flags = __bswap_32(cFrame.payload.flags);
    
    // Fill the data array
    const int8_t *fp;
    T *a;
    a = (T *) PyArray_DATA(data);
    for(i=0; i<4096; i+=4) {
        fp = drxLUT[ cFrame.payload.bytes[i] ];
        *(a + 2*i + 0) = fp[0];
        *(a + 2*i + 1) = fp[1];
        fp = drxLUT[ cFrame.payload.bytes[i+1] ];
        *(a + 2*i + 2) = fp[0];
        *(a + 2*i + 3) = fp[1];
        fp = drxLUT[ cFrame.payload.bytes[i+2] ];
        *(a + 2*i + 4) = fp[0];
        *(a + 2*i + 5) = fp[1];
        fp = drxLUT[ cFrame.payload.bytes[i+3] ];
        *(a + 2*i + 6) = fp[0];
        *(a + 2*i + 7) = fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.sync_word) ) {
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to the frame object
    // 1. Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = PyLong_FromUnsignedLong(cFrame.header.frame_count);
    PyObject_SetAttrString(fHeader, "frame_count", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("B", cFrame.header.id);
    PyObject_SetAttrString(fHeader, "drx_id", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.second_count);
    PyObject_SetAttrString(fHeader, "second_count", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.decimation);
    PyObject_SetAttrString(fHeader, "decimation", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.time_offset);
    PyObject_SetAttrString(fHeader, "time_offset", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromUnsignedLongLong(cFrame.payload.timetag);
    PyObject_SetAttrString(fPayload, "timetag", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.payload.tuning_word);
    PyObject_SetAttrString(fPayload, "tuning_word", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.payload.flags);
    PyObject_SetAttrString(fPayload, "flags", temp);
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


PyObject *read_drx_cf32(PyObject *self, PyObject *args) {
    return read_drx<float,NPY_COMPLEX64>(self, args);
}

char read_drx_cf32_doc[] = PyDoc_STR(\
"Function to read in a single DRX frame (header+payload) and store the contents\n\
as a Frame object.\n\
");


PyObject *read_drx_ci8(PyObject *self, PyObject *args) {
    return read_drx<int8_t,NPY_INT8>(self, args);
}

char read_drx_ci8_doc[] = PyDoc_STR(\
"Function to read in a single DRX frame (header+payload) and store the contents\n\
as a Frame object.\n\
\n\
.. note::\n\
    This function differs from `read_drx` in that it returns a\n\
    `lsl.drx.FramePayload` containing a 2-D numpy.int8 array (samples by real/\n\
    complex) rather than a 1-D numpy.complex64 array.\n\
\n\
.. versionadded:: 2.1.3\n\
");
