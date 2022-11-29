#include "Python.h"
#include <cmath>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"

/*
  TBN Reader
*/

typedef struct __attribute__((packed)) {
    uint32_t sync_word;
    union {
        struct {
            uint32_t frame_count:24;
            uint8_t  id:8;
        };
        uint32_t frame_count_word;
    };
    uint32_t tuning_word;
    union {
        uint16_t tbn_id;
        struct {
            uint16_t stand:14;
            uint8_t  reserved:1;
            uint8_t  is_tbw:1;
        };
    };
    uint16_t gain;
} TBNHeader;


typedef struct __attribute__((packed)) {
    uint64_t timetag;
    uint8_t  bytes[1024];
} TBNPayload;


typedef struct __attribute__((packed)) {
    TBNHeader header;
    TBNPayload payload;
} TBNFrame;


PyObject *tbn_method = NULL;
PyObject *tbn_size   = NULL;


template<typename T, NPY_TYPES N>
PyObject *read_tbn(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i;
    TBNFrame cFrame;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[1];
    dims[0] = (npy_intp) 512;
    if( N == NPY_INT8) dims[0] *= 2;
    data = (PyArrayObject*) PyArray_ZEROS(1, dims, N, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Read from the file
    if( tbn_method == NULL ) {
        tbn_method = Py_BuildValue("s", "read");
        tbn_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, tbn_method, tbn_size, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        Py_XDECREF(data);
        return NULL;
    } else if( PyString_GET_SIZE(buffer) != sizeof(cFrame) ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(data);
        Py_XDECREF(buffer);
        return NULL;
    }
    memcpy(&cFrame, PyString_AS_STRING(buffer), sizeof(cFrame));
    Py_XDECREF(buffer);
    
    Py_BEGIN_ALLOW_THREADS
    
    // Swap the bits around
    cFrame.header.frame_count_word = __bswap_32(cFrame.header.frame_count_word);
    cFrame.header.tuning_word = __bswap_32(cFrame.header.tuning_word);
    cFrame.header.tbn_id = __bswap_16(cFrame.header.tbn_id);
    cFrame.header.gain= __bswap_16(cFrame.header.gain);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    
    // Fill the data array
    T *a;
    a = (T *) PyArray_DATA(data);
    for(i=0; i<512; i++) {
        *(a + 2*i + 0) = tbnLUT[ cFrame.payload.bytes[2*i+0] ];
        *(a + 2*i + 1) = tbnLUT[ cFrame.payload.bytes[2*i+1] ];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.sync_word) ) {
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to the frame object
    // 1.  Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = PyLong_FromUnsignedLong(cFrame.header.frame_count);
    PyObject_SetAttrString(fHeader, "frame_count", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.tuning_word);
    PyObject_SetAttrString(fHeader, "tuning_word", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.tbn_id);
    PyObject_SetAttrString(fHeader, "tbn_id", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.gain);
    PyObject_SetAttrString(fHeader, "gain", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromUnsignedLongLong(cFrame.payload.timetag);
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


PyObject *read_tbn_cf32(PyObject *self, PyObject *args) {
    return read_tbn<float,NPY_COMPLEX64>(self, args);
}

char read_tbn_cf32_doc[] = PyDoc_STR(\
"Function to read in a single TBN frame (header+payload) and store the contents\n\
as a Frame object.\n\
");


PyObject *read_tbn_ci8(PyObject *self, PyObject *args) {
    return read_tbn<int8_t,NPY_INT8>(self, args);
}

char read_tbn_ci8_doc[] = PyDoc_STR(\
"Function to read in a single TBN frame (header+payload) and store the contents\n\
as a Frame object.\n\
\n\
.. note::\n\
    This function differs from `read_tbn` in that it returns a\n\
    `lsl.tbn.FramePayload` containing a 2-D numpy.int8 array (samples by real/\n\
    complex) rather than a 1-D numpy.complex64 array.\n\
\n\
.. versionadded:: 2.1.3\n\
");
