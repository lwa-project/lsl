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
    union {
        float vis[2048];   // Really std::complex<float> vis[256*4]
        uint32_t raw[2048];
    };
} CORPayload;


typedef struct __attribute__((packed)) {
    CORHeader header;
    CORPayload payload;
} CORFrame;


static int cached_nchan = 192;

PyObject *cor_method = NULL;
PyObject *cor_size   = NULL;


template<int NCHAN>
PyObject *read_cor_impl(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data=NULL;
    int i, nchan;
    CORFrame cFrame;
    
    static constexpr int frameSize = sizeof(CORHeader)+8+4+2+2+NCHAN*8*4;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }
    
    // Create the output data array
    npy_intp dims[3];
    // 32+32-bit Data -> as COR_NCHAN channels, 2 pols @ stand 0, 2 pols @  stand 1
    dims[0] = NCHAN;
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
    }
    if( cor_size == NULL ) {
        cor_size = Py_BuildValue("i", frameSize);
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
    } else if( PyBytes_GET_SIZE(buffer) != frameSize ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(data);
        Py_XDECREF(buffer);
        return NULL;
    }
    memcpy(&cFrame, PyBytes_AS_STRING(buffer), frameSize);
    Py_XDECREF(buffer);
    
    // Determine the number of channels in the frame
    if( (NCHAN > 33) && validSync5C((uint32_t) cFrame.payload.raw[33*4*2]) ) {
        nchan = 33;
    } else if( (NCHAN > 72) && validSync5C((uint32_t) cFrame.payload.raw[72*4*2]) ) {
        nchan = 72;
    } else if( (NCHAN > 112) && validSync5C((uint32_t) cFrame.payload.raw[112*4*2]) ) {
        nchan = 112;
    } else {
        nchan = NCHAN;
    }
    
    // If nchan is not what we expect, update the cache and retry with correct size
    if( nchan != NCHAN ) {
        cached_nchan = nchan;   // Update for next time
        Py_XDECREF(cor_size);   // Update for next time
        cor_size = NULL;        // Force NULL since Py_XDECREF isn't guaranteed to change cor_size
        PyObject_CallMethod(ph, "seek", "ii", -frameSize, 1);
        
        switch(nchan) {
            case 33:
                return read_cor_impl<33>(self,args);
            case 72:
                return read_cor_impl<72>(self,args);
            case 112:
                return read_cor_impl<112>(self,args);
            case 192:
                return read_cor_impl<192>(self,args);
            default:
                PyErr_Format(SyncError, "Mark 5C sync word differs from expected during mode change");
                return NULL;
        }
    }
    
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
    memcpy(a, &cFrame.payload.vis, sizeof(std::complex<float>)*NCHAN*4);
    if( cFrame.payload.stand0 == cFrame.payload.stand1 ) {
        // Deal with the edge of the triangular matrix that ADP outputs
        // so that we do not get strange values in the output.  These are
        // all auto-correlations and we can just use conjgation to get YX
        // from XY.
        for(i=0; i<NCHAN; i++) {
            a[4*i + 1*2 + 0] = std::conj(a[4*i + 0*2 + 1]);
        }
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.sync_word) ) {
        cached_nchan = 192;   // Update for next time
        Py_XDECREF(cor_size);   // Update for next time
        cor_size = NULL;
        PyObject_CallMethod(ph, "seek", "ii", -frameSize, 1);
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

PyObject *read_cor(PyObject *self, PyObject *args) {
    // Try fast path first with cached nchan
    switch(cached_nchan) {
        case 33:
            return read_cor_impl<33>(self,args);
        case 72:
            return read_cor_impl<72>(self,args);
        case 112:
            return read_cor_impl<112>(self,args);
        case 192:
        default:
            return read_cor_impl<192>(self,args);
    }
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
