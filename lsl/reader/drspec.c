#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"

/*
  DR Spectrometer Reader
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
    unsigned int MAGIC1; 			// must always equal 0xC0DEC0DE
    unsigned long long timeTag0;		// time tag of first frame in ``block''
    unsigned short int timeOffset;		// time offset reported by DP
    unsigned short int decFactor;		// decimation factor
    unsigned int freqCode[2]; 		// DP frequency codes for each tuning
                        //   indexing: 0..1 = Tuning 1..2
    unsigned int fills[4]; 			// fills for each pol/tuning combination
                        //   indexing: 0..3 = X0, Y0 X1, Y1
    unsigned char errors[4]; 		// error flag for each pol/tuning combo
                        //   indexing: 0..3 = X0, Y0 X1, Y1
    unsigned char beam;			// beam number
    unsigned char stokes_format; 		// ouptut format
    unsigned char spec_version;		// version of the spectrometer data file
    unsigned char flags;			// flag bit-field
    unsigned int nFreqs;			// <Transform Length>
    unsigned int nInts;			// <Integration Count>
    unsigned int satCount[4];		// saturation count for each pol/tuning combo
                        //   indexing: 0..3 = X0, Y0 X1, Y1
    unsigned int MAGIC2;			// must always equal 0xED0CED0C
} DRSpecHeader;
#pragma pack(pop)


PyObject *drspec_method   = NULL;
PyObject *drspec_size_hdr = NULL;
PyObject *drspec_size_dat = NULL;


static void parse_linear_single(DRSpecHeader header, float *data, float *S0, float *S1, int isX, int isY) {
    int i;
    float norm0, norm1;
    
    // Spectra normalization factors
    if( isX == 1 && isY == 0) {
        // XX*
        norm0 = header.nFreqs * header.fills[0];
        norm1 = header.nFreqs * header.fills[2];
    } else if( isX == 0 && isY == 1 ) {
        // YY*
        norm0 = header.nFreqs * header.fills[1];
        norm1 = header.nFreqs * header.fills[3];
    } else {
        // real(XY*) or imag(XY*)
        norm0 = header.nFreqs * min(header.fills[0], header.fills[1]);
        norm1 = header.nFreqs * min(header.fills[2], header.fills[3]);
    }
    
    // Sort out the data
    for(i=0; i<header.nFreqs; i++) {
        // XX*/real(XY*)/imag(XY*)/YY only
        *(S0 + i) = *(data + 0*header.nFreqs + i) / norm0;
        *(S1 + i) = *(data + 1*header.nFreqs + i) / norm1;
    }
}


static void parse_linear_half(DRSpecHeader header, float *data, float *XX0, float *XX1, float *YY0, float *YY1) {
    int i;
    float normXX0, normXX1, normYY0, normYY1;
    
    // Spectra normalization factors
    normXX0 = header.nFreqs * header.fills[0];
    normXX1 = header.nFreqs * header.fills[2];
    normYY0 = header.nFreqs * header.fills[1];
    normYY1 = header.nFreqs * header.fills[3];
    
    // Sort out the data
    for(i=0; i<header.nFreqs; i++) {
        // XX*
        *(XX0 + i) = *(data + 0*header.nFreqs + 2*i + 0) / normXX0;
        *(XX1 + i) = *(data + 2*header.nFreqs + 2*i + 0) / normXX1;
        
        // YY*
        *(YY0 + i) = *(data + 0*header.nFreqs + 2*i + 1) / normYY0;
        *(YY1 + i) = *(data + 2*header.nFreqs + 2*i + 1) / normYY1;
    }
}


static void parse_linear_other_half(DRSpecHeader header, float *data, float *CR0, float *CR1, float *CI0, float *CI1) {
    int i;
    float normCH0, normCH1;
    
    // Spectra normalization factors
    normCH0 = header.nFreqs * min(header.fills[0], header.fills[1]);
    normCH1 = header.nFreqs * min(header.fills[2], header.fills[3]);
    
    // Sort out the data
    for(i=0; i<header.nFreqs; i++) {
        // real(XY*)
        *(CR0 + i) = *(data + 0*header.nFreqs + 2*i + 0) / normCH0;
        *(CR1 + i) = *(data + 2*header.nFreqs + 2*i + 0) / normCH1;
        
        // imag(XY*)
        *(CI0 + i) = *(data + 0*header.nFreqs + 2*i + 1) / normCH0;
        *(CI1 + i) = *(data + 2*header.nFreqs + 2*i + 1) / normCH1;
    }
}


static void parse_linear_full(DRSpecHeader header, float *data, float *XX0, float *XX1, float *CR0, float *CR1, float *CI0, float *CI1, float *YY0, float *YY1) {
    int i;
    float normXX0, normXX1, normCH0, normCH1, normYY0, normYY1;
    
    // Spectra normalization factors
    normXX0 = header.nFreqs * header.fills[0];
    normXX1 = header.nFreqs * header.fills[2];
    normCH0 = header.nFreqs * min(header.fills[0], header.fills[1]);
    normCH1 = header.nFreqs * min(header.fills[2], header.fills[3]);
    normYY0 = header.nFreqs * header.fills[1];
    normYY1 = header.nFreqs * header.fills[3];

    // Sort out the data
    for(i=0; i<header.nFreqs; i++) {
        // XX*
        *(XX0 + i) = *(data + 0*header.nFreqs + 4*i + 0) / normXX0;
        *(XX1 + i) = *(data + 4*header.nFreqs + 4*i + 0) / normXX1;
        
        // real(XY*)
        *(CR0 + i) = *(data + 0*header.nFreqs + 4*i + 1) / normCH0;
        *(CR1 + i) = *(data + 4*header.nFreqs + 4*i + 1) / normCH1;
        
        // imag(XY*)
        *(CI0 + i) = *(data + 0*header.nFreqs + 4*i + 2) / normCH0;
        *(CI1 + i) = *(data + 4*header.nFreqs + 4*i + 2) / normCH1;
        
        // YY*
        *(YY0 + i) = *(data + 0*header.nFreqs + 4*i + 3) / normYY0;
        *(YY1 + i) = *(data + 4*header.nFreqs + 4*i + 3) / normYY1;
    }
}


static void parse_stokes_single(DRSpecHeader header, float *data, float *S0, float *S1) {
    int i;
    float norm0, norm1;
    
    // Spectra normalization factors
    norm0 = header.nFreqs * min(header.fills[0], header.fills[1]);
    norm1 = header.nFreqs * min(header.fills[2], header.fills[3]);
    
    // Sort out the data
    for(i=0; i<header.nFreqs; i++) {
        // I/Q/U/V only
        *(S0 + i) = *(data + 0*header.nFreqs + i) / norm0;
        *(S1 + i) = *(data + 1*header.nFreqs + i) / norm1;
    }
}


static void parse_stokes_half(DRSpecHeader header, float *data, float *I0, float *I1, float *V0, float *V1) {
    int i;
    float norm0, norm1;
    
    // Spectra normalization factors
    norm0 = header.nFreqs * min(header.fills[0], header.fills[1]);
    norm1 = header.nFreqs * min(header.fills[2], header.fills[3]);
    
    // Sort out the data
    for(i=0; i<header.nFreqs; i++) {
        // I
        *(I0 + i) = *(data + 0*header.nFreqs + 2*i + 0) / norm0;
        *(I1 + i) = *(data + 2*header.nFreqs + 2*i + 0) / norm1;
        
        // V
        *(V0 + i) = *(data + 0*header.nFreqs + 2*i + 1) / norm0;
        *(V1 + i) = *(data + 2*header.nFreqs + 2*i + 1) / norm1;
    }
}


static void parse_stokes_full(DRSpecHeader header, float *data, float *I0, float *I1, float *Q0, float *Q1, float *U0, float *U1, float *V0, float *V1) {
    int i;
    float norm0, norm1;
    
    // Spectra normalization factors
    norm0 = header.nFreqs * min(header.fills[0], header.fills[1]);
    norm1 = header.nFreqs * min(header.fills[2], header.fills[3]);
    
    // Sort out the data
    for(i=0; i<header.nFreqs; i++) {
        // I
        *(I0 + i) = *(data + 0*header.nFreqs + 4*i + 0) / norm0;
        *(I1 + i) = *(data + 4*header.nFreqs + 4*i + 0) / norm1;
        
        // Q
        *(Q0 + i) = *(data + 0*header.nFreqs + 4*i + 1) / norm0;
        *(Q1 + i) = *(data + 4*header.nFreqs + 4*i + 1) / norm1;
        
        // U
        *(U0 + i) = *(data + 0*header.nFreqs + 4*i + 2) / norm0;
        *(U1 + i) = *(data + 4*header.nFreqs + 4*i + 2) / norm1;
        
        // V
        *(V0 + i) = *(data + 0*header.nFreqs + 4*i + 3) / norm0;
        *(V1 + i) = *(data + 4*header.nFreqs + 4*i + 3) / norm1;
    }
}


PyObject *read_drspec(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyObject *tuningWords, *fills, *errors, *saturations;
    PyArrayObject *dataA0=NULL, *dataB0=NULL, *dataC0=NULL, *dataD0=NULL;
    PyArrayObject *dataA1=NULL, *dataB1=NULL, *dataC1=NULL, *dataD1=NULL;
    int i, nSets;
    DRSpecHeader header;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Read in a single header from the file
    if( drspec_method == NULL ) {
        drspec_method = Py_BuildValue("s", "read");
        drspec_size_hdr = Py_BuildValue("i", sizeof(header));
    }
    buffer = PyObject_CallMethodObjArgs(ph, drspec_method, drspec_size_hdr, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() drspec_method");
        }
        goto fail;
    } else if( PyString_GET_SIZE(buffer) != sizeof(DRSpecHeader) ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        goto fail;
    }
    memcpy(&header, PyString_AS_STRING(buffer), sizeof(header));
    Py_XDECREF(buffer);
    
    // Check the header's magic numbers
    if( header.MAGIC1 != 0xC0DEC0DE || header.MAGIC2 != 0xED0CED0C) {
        buffer = PyObject_CallMethod(ph, "seek", "ii", -sizeof(DRSpecHeader), 1);
        PyErr_Format(SyncError, "Sync word differs from expected");
        goto fail;
    }
    
    // Get the data format
    if( header.stokes_format < 0x10 ) {
        // Linear
        if( header.stokes_format < 0x09 && header.stokes_format != 0x06 ) {
            nSets = 2;
        } else if( header.stokes_format == 0x06 || header.stokes_format == 0x09 ) {
            nSets = 4;
        } else {
            nSets = 8;
        }
    } else {
        // Stokes
        if( header.stokes_format < 0x90 ) {
            nSets = 2;
        } else if( header.stokes_format == 0x90 ) {
            nSets = 4;
        } else {
            nSets = 8;
        }
    }
    
    // Create the output data arrays
    npy_intp dims[1];
    dims[0] = (npy_intp) header.nFreqs;
    
    dataA0 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataA0 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - XX0/I0");
        goto fail;
    }
    
    dataB0 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataB0 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - CR0/Q0");
        goto fail;
    }
    
    dataC0 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataC0 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - CI0/U0");
        goto fail;
    }
    
    dataD0 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataD0 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - YY0/V0");
        goto fail;
    }
    
    dataA1 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataA1 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - XX1/I1");
        goto fail;
    }
    
    dataB1 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataB1 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - CR1/Q1");
        goto fail;
    }
    
    dataC1 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataC1 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - CI1/U1");
        goto fail;
    }
    
    dataD1 = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataD1 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - YY1/V1");
        goto fail;
    }
    
    // Get ready to read in the data
    float *data;
    data = (float *) malloc(sizeof(float)*nSets*header.nFreqs);
    
    // Read in the data section
    if( drspec_size_dat == NULL ) {
        drspec_size_dat = Py_BuildValue("i", sizeof(float)*nSets*header.nFreqs);
    } else if( PyInt_AsLong(drspec_size_dat) != sizeof(float)*nSets*header.nFreqs ) {
        Py_XDECREF(drspec_size_dat);
        drspec_size_dat = Py_BuildValue("i", sizeof(float)*nSets*header.nFreqs);
    }
    buffer = PyObject_CallMethodObjArgs(ph, drspec_method, drspec_size_dat, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() drspec_method");
        }
        goto fail;
    } else if( PyString_GET_SIZE(buffer) != sizeof(float)*nSets*header.nFreqs ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        goto fail;
    }
    memcpy(data, PyString_AS_STRING(buffer), sizeof(float)*nSets*header.nFreqs);
    Py_XDECREF(buffer);
    
    Py_BEGIN_ALLOW_THREADS
    
    // Fill the data arrays
    float *a0, *b0, *c0, *d0, *a1, *b1, *c1, *d1;
    a0 = (float *) PyArray_DATA(dataA0);
    b0 = (float *) PyArray_DATA(dataB0);
    c0 = (float *) PyArray_DATA(dataC0);
    d0 = (float *) PyArray_DATA(dataD0);
    a1 = (float *) PyArray_DATA(dataA1);
    b1 = (float *) PyArray_DATA(dataB1);
    c1 = (float *) PyArray_DATA(dataC1);
    d1 = (float *) PyArray_DATA(dataD1);
    if( header.stokes_format < 0x10 ) {
        // Linear
        if( header.stokes_format == 0x01 ) {
            // XX* only
            parse_linear_single(header, data, a0, a1, 1, 0);
        } else if( header.stokes_format == 0x02 ) {
            // real(XY*) only
            parse_linear_single(header, data, b0, b1, 1, 1);
        } else if( header.stokes_format == 0x04 ) {
            // imag(XY*) only
            parse_linear_single(header, data, c0, c1, 1, 1);
        } else if( header.stokes_format == 0x08 ) {
            // YY* only
            parse_linear_single(header, data, d0, d1, 0, 1);
        } else if( header.stokes_format == 0x09 ) {
            // XX* and YY*
            parse_linear_half(header, data, a0, a1, d0, d1);
        } else if( header.stokes_format == 0x06 ) {
            // real(XY*) and imag(XY*)
            parse_linear_other_half(header, data, b0, b1, c0, c1);
        } else {
            // XX*, real(XY*), imag(XY*), and YY*
            parse_linear_full(header, data, a0, a1, b0, b1, c0, c1, d0, d1);
        }
    } else {
        // Stokes
        if( header.stokes_format == 0x10 ) {
            // I only
            parse_stokes_single(header, data, a0, a1);
        } else if( header.stokes_format == 0x20 ) {
            // Q only
            parse_stokes_single(header, data, b0, b1);
        } else if( header.stokes_format == 0x40 ) {
            // U only
            parse_stokes_single(header, data, c0, c1);
        } else if( header.stokes_format == 0x80 ) {
            // V only
            parse_stokes_single(header, data, d0, d1);
        } else if( header.stokes_format == 0x90 ) {
            // I and V
            parse_stokes_half(header, data, a0, a1, d0, d1);
        } else {
            // I, Q, U, and V
            parse_stokes_full(header, data, a0, a1, b0, b1, c0, c1, d0, d1);
        }
    }
    
    // Clean up the data
    free(data);
    
    Py_END_ALLOW_THREADS
    
    // Fill in the multi-value fields (tuningWords, fills, errors)
    tuningWords = PyList_New(2);
    if(tuningWords == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output list - tuningWords");
        Py_XDECREF(tuningWords);
        goto fail;
    }
    for(i=0; i<2; i++) {
        temp = Py_BuildValue("I", header.freqCode[i]);
        PyList_SetItem(tuningWords, i, temp);
    }
    
    fills = PyList_New(4);
    if(fills == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output list - fills");
        Py_XDECREF(tuningWords);
        Py_XDECREF(fills);
        goto fail;
    }
    for(i=0; i<4; i++) {
        temp = Py_BuildValue("I", header.fills[i]);
        PyList_SetItem(fills, i, temp);
    }
    
    errors = PyList_New(4);
    if(errors == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output list - flags");
        Py_XDECREF(tuningWords);
        Py_XDECREF(fills);
        Py_XDECREF(errors);
        goto fail;
    }
    for(i=0; i<4; i++) {
        temp = Py_BuildValue("B", header.errors[i]);
        PyList_SetItem(errors, i, temp);
    }
    
    saturations = PyList_New(4);
    if(saturations == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output list - saturations");
        Py_XDECREF(tuningWords);
        Py_XDECREF(fills);
        Py_XDECREF(errors);
        Py_XDECREF(saturations);
        goto fail;
    }
    for(i=0; i<4; i++) {
        temp = Py_BuildValue("I", header.satCount[i]);
        PyList_SetItem(saturations, i, temp);
    }
    
    // Save the data to the frame object
    // 1. Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = Py_BuildValue("B", header.beam);
    PyObject_SetAttrString(fHeader, "beam", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("B", header.stokes_format);
    PyObject_SetAttrString(fHeader, "format", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", header.decFactor);
    PyObject_SetAttrString(fHeader, "decimation", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", header.timeOffset);
    PyObject_SetAttrString(fHeader, "time_offset", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", header.nInts);
    PyObject_SetAttrString(fHeader, "nints", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromUnsignedLongLong(header.timeTag0);
    PyObject_SetAttrString(fPayload, "timetag", temp);
    Py_XDECREF(temp);
    
    PyObject_SetAttrString(fPayload, "tuning_words", tuningWords);
    
    PyObject_SetAttrString(fPayload, "fills", fills);
    
    PyObject_SetAttrString(fPayload, "errors", errors);
    
    PyObject_SetAttrString(fPayload, "saturations", saturations);
    
    // Linear
    if( header.stokes_format & 0x01 ) {
        PyObject_SetAttrString(fPayload, "XX0", PyArray_Return(dataA0));
        PyObject_SetAttrString(fPayload, "XX1", PyArray_Return(dataA1));
    }
    if( header.stokes_format & 0x02 ) {
        PyObject_SetAttrString(fPayload, "XY_real0", PyArray_Return(dataB0));
        PyObject_SetAttrString(fPayload, "XY_real1", PyArray_Return(dataB1));
    }
    if( header.stokes_format & 0x04 ) {
        PyObject_SetAttrString(fPayload, "XY_imag0", PyArray_Return(dataC0));
        PyObject_SetAttrString(fPayload, "XY_imag1", PyArray_Return(dataC1));
    }
    if( header.stokes_format & 0x08 ) {
        PyObject_SetAttrString(fPayload, "YY0", PyArray_Return(dataD0));
        PyObject_SetAttrString(fPayload, "YY1", PyArray_Return(dataD1));
    }
    
    // Stokes
    if( header.stokes_format & 0x10 ) {
        PyObject_SetAttrString(fPayload, "I0", PyArray_Return(dataA0));
        PyObject_SetAttrString(fPayload, "I1", PyArray_Return(dataA1));
    }
    if( header.stokes_format & 0x20 ) {
        PyObject_SetAttrString(fPayload, "Q0", PyArray_Return(dataB0));
        PyObject_SetAttrString(fPayload, "Q1", PyArray_Return(dataB1));
    }
    if( header.stokes_format & 0x40 ) {
        PyObject_SetAttrString(fPayload, "U0", PyArray_Return(dataC0));
        PyObject_SetAttrString(fPayload, "U1", PyArray_Return(dataC1));
    }
    if( header.stokes_format & 0x80 ) {
        PyObject_SetAttrString(fPayload, "V0", PyArray_Return(dataD0));
        PyObject_SetAttrString(fPayload, "V1", PyArray_Return(dataD1));
    }
    
    // 3. Frame
    PyObject_SetAttrString(frame, "header", fHeader);
    PyObject_SetAttrString(frame, "payload", fPayload);
    output = Py_BuildValue("O", frame);
    
    Py_XDECREF(fHeader);
    Py_XDECREF(tuningWords);
    Py_XDECREF(fills);
    Py_XDECREF(errors);
    Py_XDECREF(saturations);
    Py_XDECREF(fPayload);
    Py_XDECREF(dataA0);
    Py_XDECREF(dataB0);
    Py_XDECREF(dataC0);
    Py_XDECREF(dataD0);
    Py_XDECREF(dataA1);
    Py_XDECREF(dataB1);
    Py_XDECREF(dataC1);
    Py_XDECREF(dataD1);
    
    return output;
    
fail:
    Py_XDECREF(buffer);
    Py_XDECREF(dataA0);
    Py_XDECREF(dataB0);
    Py_XDECREF(dataC0);
    Py_XDECREF(dataD0);
    Py_XDECREF(dataA1);
    Py_XDECREF(dataB1);
    Py_XDECREF(dataC1);
    Py_XDECREF(dataD1);
    
    return NULL;
}

char read_drspec_doc[] = PyDoc_STR(\
"Function to read in a DR spectrometer header structure and data and return\n\
a drspec.Frame instance.\n\
\n\
.. note::\n\
\tThis function normalizes the spectra by the number of relevant fills.  For\n\
\tproducts that are a function of more than one primary input, i.e., real(XY*) or\n\
\tI, the minimum fill of X and Y are used for normalization.\n\
");
