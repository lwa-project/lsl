#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#if defined(__linux__)
/* Linux */
#include <byteswap.h>
#if defined(__APPLE__) && defined(__MACH__)
/* OSX */
#include <libkern/OSByteOrder.h>
#define __bswap_16 OSSwapInt16
#define __bswap_32 OSSwapInt32
#define __bswap_64 OSSwapInt64
#else
/* BSD */
#include <sys/endian.h>
#define __bswap_16 __bswap16_var
#define __bswap_32 __bswap32_var
#define __bswap_64 __bswap64_var
#endif

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


static void parseLinearSingle(DRSpecHeader header, float *data, float *S0, float *S1, int isX, int isY) {
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
		// XY* or YX*
		norm0 = header.nFreqs * min(header.fills[0], header.fills[1]);
		norm1 = header.nFreqs * min(header.fills[2], header.fills[3]);
	}
	
	// Sort out the data
	for(i=0; i<header.nFreqs; i++) {
		// XX*/XY*/YX*/YY only
		*(S0 + i) = *(data + 0*header.nFreqs + i) / norm0;
		*(S1 + i) = *(data + 1*header.nFreqs + i) / norm1;
	}
}


static void parseLinearHalf(DRSpecHeader header, float *data, float *XX0, float *XX1, float *YY0, float *YY1) {
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
		*(XX0 + i) = *(data + 0*header.nFreqs + i + 0) / normXX0;
		*(XX1 + i) = *(data + 2*header.nFreqs + i + 0) / normXX1;
		
		// YY*
		*(YY0 + i) = *(data + 0*header.nFreqs + i + 1) / normYY0;
		*(YY1 + i) = *(data + 2*header.nFreqs + i + 1) / normYY1;
	}
}


static void parseLinearFull(DRSpecHeader header, float *data, float *XX0, float *XX1, float *XY0, float *XY1, float *YX0, float *YX1, float *YY0, float *YY1) {
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
		*(XX0 + i) = *(data + 0*header.nFreqs + i + 0) / normXX0;
		*(XX1 + i) = *(data + 4*header.nFreqs + i + 0) / normXX1;
		
		// XY*
		*(XY0 + i) = *(data + 0*header.nFreqs + i + 1) / normCH0;
		*(XY1 + i) = *(data + 4*header.nFreqs + i + 1) / normCH1;
		
		// YX*
		*(YX0 + i) = *(data + 0*header.nFreqs + i + 2) / normCH0;
		*(YX1 + i) = *(data + 4*header.nFreqs + i + 2) / normCH1;
		
		// YY*
		*(YY0 + i) = *(data + 0*header.nFreqs + i + 3) / normYY0;
		*(YY1 + i) = *(data + 4*header.nFreqs + i + 3) / normYY1;
	}
}


static void parseStokesSingle(DRSpecHeader header, float *data, float *S0, float *S1) {
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


static void parseStokesHalf(DRSpecHeader header, float *data, float *I0, float *I1, float *V0, float *V1) {
	int i;
	float norm0, norm1;
	
	// Spectra normalization factors
	norm0 = header.nFreqs * min(header.fills[0], header.fills[1]);
	norm1 = header.nFreqs * min(header.fills[2], header.fills[3]);
	
	// Sort out the data
	for(i=0; i<header.nFreqs; i++) {
		// I
		*(I0 + i) = *(data + 0*header.nFreqs + i + 0) / norm0;
		*(I1 + i) = *(data + 2*header.nFreqs + i + 0) / norm1;
		
		// V
		*(V0 + i) = *(data + 0*header.nFreqs + i + 1) / norm0;
		*(V1 + i) = *(data + 2*header.nFreqs + i + 1) / norm1;
	}
}


static void parseStokesFull(DRSpecHeader header, float *data, float *I0, float *I1, float *Q0, float *Q1, float *U0, float *U1, float *V0, float *V1) {
	int i;
	float norm0, norm1;
	
	// Spectra normalization factors
	norm0 = header.nFreqs * min(header.fills[0], header.fills[1]);
	norm1 = header.nFreqs * min(header.fills[2], header.fills[3]);
	
	// Sort out the data
	for(i=0; i<header.nFreqs; i++) {
		// I
		*(I0 + i) = *(data + 0*header.nFreqs + i + 0) / norm0;
		*(I1 + i) = *(data + 4*header.nFreqs + i + 3) / norm1;
		
		// Q
		*(Q0 + i) = *(data + 0*header.nFreqs + i + 1) / norm0;
		*(Q1 + i) = *(data + 4*header.nFreqs + i + 1) / norm1;
		
		// U
		*(U0 + i) = *(data + 0*header.nFreqs + i + 2) / norm0;
		*(U1 + i) = *(data + 4*header.nFreqs + i + 2) / norm1;
		
		// V
		*(V0 + i) = *(data + 0*header.nFreqs + i + 3) / norm0;
		*(V1 + i) = *(data + 4*header.nFreqs + i + 3) / norm1;
	}
}


PyObject *readDRSpec(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyObject *tuningWords, *fills, *errors, *saturations;
	PyArrayObject *dataA0, *dataB0, *dataC0, *dataD0;
	PyArrayObject *dataA1, *dataB1, *dataC1, *dataD1;
	
	int i, nSets;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Read in a single header
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	DRSpecHeader header;
	i = fread(&header, sizeof(DRSpecHeader), 1, fh);
	if(ferror(fh)) {
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
		return NULL;
	}
	if(feof(fh)) {
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(eofError, "End of file encountered during filehandle read");
		return NULL;
	}
	
	// Check the header's magic numbers
	if( header.MAGIC1 != 0xC0DEC0DE || header.MAGIC2 != 0xED0CED0C) {
		fseek(fh, -sizeof(DRSpecHeader), SEEK_CUR);
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(syncError, "Sync word differs from expected");
		return NULL;
	}
	
	// Get the data format
	if( header.stokes_format < 0x10 ) {
		// Linear
		if( header.stokes_format < 0x09 ) {
			nSets = 2;
		} else if( header.stokes_format == 0x09 ) {
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
	
	dataA0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataA0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XX0/I0");
		Py_XDECREF(dataA0);
		return NULL;
	}
	
	dataB0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataB0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XY0/Q0");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		return NULL;
	}
	
	dataC0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataC0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YX0/U0");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		return NULL;
	}
	
	dataD0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataD0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YY0/V0");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		return NULL;
	}
	
	dataA1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataA1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XX1/I1");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		return NULL;
	}
	
	dataB1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataB1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XY1/Q1");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		return NULL;
	}
	
	dataC1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataC1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YX1/U1");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		return NULL;
	}
	
	dataD1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataD1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YY1/V1");
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
	
	// Get ready to read in the data
	float *data;
	data = (float *) malloc(sizeof(float)*nSets*header.nFreqs);
	
	// Read in the data section
	i = fread(data, sizeof(float), nSets*header.nFreqs, fh);
	if(ferror(fh)) {
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
		return NULL;
	}
	if(feof(fh)) {
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(eofError, "End of file encountered during filehandle read");
		return NULL;
	}
	PyFile_DecUseCount((PyFileObject *) ph);
	
	// Fill the data arrays
	float *a0, *b0, *c0, *d0, *a1, *b1, *c1, *d1;
	a0 = (float *) dataA0->data;
	b0 = (float *) dataB0->data;
	c0 = (float *) dataC0->data;
	d0 = (float *) dataD0->data;
	a1 = (float *) dataA1->data;
	b1 = (float *) dataB1->data;
	c1 = (float *) dataC1->data;
	d1 = (float *) dataD1->data;
	if( header.stokes_format < 0x10 ) {
		// Linear
		if( header.stokes_format == 0x01 ) {
			// XX* only
			parseLinearSingle(header, data, a0, a1, 1, 0);
		} else if( header.stokes_format == 0x02 ) {
			// XY* only
			parseLinearSingle(header, data, b0, b1, 1, 1);
		} else if( header.stokes_format == 0x04 ) {
			// YX* only
			parseLinearSingle(header, data, c0, c1, 1, 1);
		} else if( header.stokes_format == 0x08 ) {
			// YY* only
			parseLinearSingle(header, data, d0, d1, 0, 1);
		} else if( header.stokes_format == 0x09 ) {
			// XX* and YY*
			parseLinearHalf(header, data, a0, a1, d0, d1);
		} else {
			// XX*, XY*, YX*, and YY*
			parseLinearFull(header, data, a0, a1, b0, b1, c0, c1, d0, d1);
		}
	} else {
		// Stokes
		if( header.stokes_format == 0x10 ) {
			// I only
			parseStokesSingle(header, data, a0, a1);
		} else if( header.stokes_format == 0x20 ) {
			// Q only
			parseStokesSingle(header, data, b0, b1);
		} else if( header.stokes_format == 0x40 ) {
			// U only
			parseStokesSingle(header, data, c0, c1);
		} else if( header.stokes_format == 0x80 ) {
			// V only
			parseStokesSingle(header, data, d0, d1);
		} else if( header.stokes_format == 0x90 ) {
			// I and V
			parseStokesHalf(header, data, a0, a1, d0, d1);
		} else {
			// I, Q, U, and V
			parseStokesFull(header, data, a0, a1, b0, b1, c0, c1, d0, d1);
		}
	}
	
	// Clean up the data
	free(data);
	
	// Fill in the multi-value fields (tuningWords, fills, errors)
	tuningWords = PyList_New(2);
	if(tuningWords == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - tuningWords");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		return NULL;
	}
	for(i=0; i<2; i++) {
		temp = Py_BuildValue("I", header.freqCode[i]);
		PyList_SetItem(tuningWords, i, temp);
	}
	
	fills = PyList_New(4);
	if(fills == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - fills");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		return NULL;
	}
	for(i=0; i<4; i++) {
		temp = Py_BuildValue("I", header.fills[i]);
		PyList_SetItem(fills, i, temp);
	}
	
	errors = PyList_New(4);
	if(errors == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - flags");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		Py_XDECREF(errors);
		return NULL;
	}
	for(i=0; i<4; i++) {
		temp = Py_BuildValue("B", header.errors[i]);
		PyList_SetItem(errors, i, temp);
	}
	
	saturations = PyList_New(4);
	if(saturations == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - saturations");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		Py_XDECREF(errors);
		Py_XDECREF(saturations);
		return NULL;
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
	PyObject_SetAttrString(fHeader, "timeOffset", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", header.nInts);
	PyObject_SetAttrString(fHeader, "nInts", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(header.timeTag0);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);
	
	PyObject_SetAttrString(fData, "tuningWords", tuningWords);
	
	PyObject_SetAttrString(fData, "fills", fills);
	
	PyObject_SetAttrString(fData, "errors", errors);
	
	PyObject_SetAttrString(fData, "saturations", saturations);
	
	// Linear
	if( header.stokes_format & 0x01 ) {
		PyObject_SetAttrString(fData, "XX0", PyArray_Return(dataA0));
		PyObject_SetAttrString(fData, "XX1", PyArray_Return(dataA1));
	}
	if( header.stokes_format & 0x02 ) {
		PyObject_SetAttrString(fData, "XY0", PyArray_Return(dataB0));
		PyObject_SetAttrString(fData, "XY1", PyArray_Return(dataB1));
	}
	if( header.stokes_format & 0x04 ) {
		PyObject_SetAttrString(fData, "YX0", PyArray_Return(dataC0));
		PyObject_SetAttrString(fData, "YX1", PyArray_Return(dataC1));
	}
	if( header.stokes_format & 0x08 ) {
		PyObject_SetAttrString(fData, "YY0", PyArray_Return(dataD0));
		PyObject_SetAttrString(fData, "YY1", PyArray_Return(dataD1));
	}
	
	// Stokes
	if( header.stokes_format & 0x10 ) {
		PyObject_SetAttrString(fData, "I0", PyArray_Return(dataA0));
		PyObject_SetAttrString(fData, "I1", PyArray_Return(dataA1));
	}
	if( header.stokes_format & 0x20 ) {
		PyObject_SetAttrString(fData, "Q0", PyArray_Return(dataB0));
		PyObject_SetAttrString(fData, "Q1", PyArray_Return(dataB1));
	}
	if( header.stokes_format & 0x40 ) {
		PyObject_SetAttrString(fData, "U0", PyArray_Return(dataC0));
		PyObject_SetAttrString(fData, "U1", PyArray_Return(dataC1));
	}
	if( header.stokes_format & 0x80 ) {
		PyObject_SetAttrString(fData, "V0", PyArray_Return(dataD0));
		PyObject_SetAttrString(fData, "V1", PyArray_Return(dataD1));
	}
	
	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(tuningWords);
	Py_XDECREF(fills);
	Py_XDECREF(errors);
	Py_XDECREF(saturations);
	Py_XDECREF(fData);
	Py_XDECREF(dataA0);
	Py_XDECREF(dataB0);
	Py_XDECREF(dataC0);
	Py_XDECREF(dataD0);
	Py_XDECREF(dataA1);
	Py_XDECREF(dataB1);
	Py_XDECREF(dataC1);
	Py_XDECREF(dataD1);
	
	output = Py_BuildValue("O", frame);
	return output;
}

char readDRSpec_doc[] = PyDoc_STR(\
"Function to read in a DR spectrometer header structure and data and return\n\
a drspec.Frame instance.\n\
\n\
.. note::\n\
\tThis function normalizes the spectra by the number of relevant fills.  For\n\
\tproducts that are a function of more than one primary input, i.e., XY* or\n\
\tI, the minimum fill of X and Y are used for normalization.\n\
");

