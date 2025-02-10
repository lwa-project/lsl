#include "Python.h"
#include <cmath>
#include <cstdint>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.hpp"

/*
  DRX8 Reader
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
} DRX8Header;


typedef struct __attribute__((packed)) {
    uint64_t timetag;
    uint32_t tuning_word;
    uint32_t flags;
    uint8_t  bytes[8192];
} DRX8Payload;


typedef struct __attribute__((packed)) {
    DRX8Header header;
    DRX8Payload payload;
} DRX8Frame;


PyObject *drx8_method = NULL;
PyObject *drx8_size   = NULL;


template<typename T, NPY_TYPES N>
PyObject *read_drx8(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i;
    DRX8Frame cFrame;
    
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
    if( drx8_method == NULL ) {
        drx8_method = Py_BuildValue("s", "read");
        drx8_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, drx8_method, drx8_size, NULL);
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
        *(a + 2*i + 0) = drx8LUT[ cFrame.payload.bytes[2*i+0] ];
        *(a + 2*i + 1) = drx8LUT[ cFrame.payload.bytes[2*i+1] ];
        *(a + 2*i + 2) = drx8LUT[ cFrame.payload.bytes[2*i+2] ];
        *(a + 2*i + 3) = drx8LUT[ cFrame.payload.bytes[2*i+3] ];
        *(a + 2*i + 4) = drx8LUT[ cFrame.payload.bytes[2*i+4] ];
        *(a + 2*i + 5) = drx8LUT[ cFrame.payload.bytes[2*i+5] ];
        *(a + 2*i + 6) = drx8LUT[ cFrame.payload.bytes[2*i+6] ];
        *(a + 2*i + 7) = drx8LUT[ cFrame.payload.bytes[2*i+7] ];
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

PyObject *read_drx8_cf32(PyObject *self, PyObject *args) {
    return read_drx8<float,NPY_COMPLEX64>(self, args);
}

char read_drx8_cf32_doc[] = PyDoc_STR(\
"Function to read in a single DRX8 frame (header+payload) and store the contents\n\
as a Frame object.\n\
");


PyObject *read_drx8_ci8(PyObject *self, PyObject *args) {
    return read_drx8<int8_t,NPY_INT8>(self, args);
}

char read_drx8_ci8_doc[] = PyDoc_STR(\
"Function to read in a single DRX8 frame (header+payload) and store the contents\n\
as a Frame object.\n\
\n\
.. note::\n\
    This function differs from `read_drx8` in that it returns a\n\
    `lsl.drx8.FramePayload` containing a 2-D numpy.int8 array (samples by real/\n\
    complex) rather than a 1-D numpy.complex64 array.\n\
\n\
");
