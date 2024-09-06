#include "Python.h"
#include <cmath>
#include <cstdint>
#include <complex>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.hpp"


/*
  COR Reader
*/

typedef struct __attribute__((packed)) {
    uint32_t sync_word;
    union {
        struct {
            uint32_t frame_count:24;
            unsigned char adp_id:8;
        };
        uint32_t frame_count_word;
    };
    uint32_t second_count;
    int16_t  first_chan;
    int16_t  gain;
} CORHeader;


typedef struct __attribute__((packed)) {
    int64_t timetag;
    int32_t navg;
    int16_t stand0;
    int16_t stand1;
    float vis[COR_NCHAN*4*2];   // Really std::complex<float> vis[COR_NCHAN*4]
} CORPayload;


typedef struct __attribute__((packed)) {
    CORHeader header;
    CORPayload payload;
} CORFrame;


PyObject *cor_method = NULL;
PyObject *cor_size   = NULL;


PyObject *read_cor(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data=NULL;
    CORFrame cFrame;
    
    int i;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }
    
    // Create the output data array
    npy_intp dims[3];
    // 32+32-bit Data -> as COR_NCHAN channels, 2 pols @ stand 0, 2 pols @  stand 1
    dims[0] = COR_NCHAN;
    dims[1] = 2;
    dims[2] = 2;
    data = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
    // Read from the file
    if( cor_method == NULL ) {
        cor_method = Py_BuildValue("s", "read");
        cor_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, cor_method, cor_size, NULL);
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
    cFrame.header.first_chan = __bswap_16(cFrame.header.first_chan);
    cFrame.header.gain = __bswap_16(cFrame.header.gain);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    cFrame.payload.navg = __bswap_32(cFrame.payload.navg);
    cFrame.payload.stand0 = __bswap_16(cFrame.payload.stand0);
    cFrame.payload.stand1 = __bswap_16(cFrame.payload.stand1);
    
    // Fill the data array
    std::complex<float> *a;
    a = (std::complex<float> *) PyArray_DATA(data);
    memcpy(a, &cFrame.payload.vis, sizeof(std::complex<float>)*COR_NCHAN*4);
    if( cFrame.payload.stand0 == cFrame.payload.stand1 ) {
        // Deal with the edge of the triangular matrix that ADP outputs
        // so that we do not get strange values in the output.  These are
        // all auto-correlations and we can just use conjgation to get YX
        // from XY.
        for(i=0; i<COR_NCHAN; i++) {
            a[4*i + 1*2 + 0] = std::conj(a[4*i + 0*2 + 1]);
        }
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.sync_word) ) {
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        goto fail;
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
    
    temp = Py_BuildValue("H", cFrame.header.first_chan);
    PyObject_SetAttrString(fHeader, "first_chan", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.gain);
    PyObject_SetAttrString(fHeader, "gain", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromLongLong(cFrame.payload.timetag);
    PyObject_SetAttrString(fPayload, "timetag", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("i", cFrame.payload.navg);
    PyObject_SetAttrString(fPayload, "navg", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("h", cFrame.payload.stand0);
    PyObject_SetAttrString(fPayload, "stand0", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("h", cFrame.payload.stand1);
    PyObject_SetAttrString(fPayload, "stand1", temp);
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
    
fail:
    Py_XDECREF(data);
    
    return NULL;
}

char read_cor_doc[] = PyDoc_STR(\
"Function to read in a single COR frame (header+payload) and store the contents\n\
as a Frame object.\n\
\n\
.. versionchanged:: 1.2.1\n\
\tUpdated readCOR for the switch over to 72 channels, complex64 data, and no\n\
\tdata weights\n\
\n\
.. versionadded:: 1.2.0\n\
");
