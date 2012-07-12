#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "numpy/arrayobject.h"

#define PI 3.1415926535898
#define imaginary _Complex_I

/*
  Exceptions for the Go Fast! Readers
*/

static PyObject *syncError;
static PyObject *eofError;


/*
  Validate a collection of Mark 5C sync words.  Return 1 is all are
  valid.
*/

int validSync(unsigned char w1, unsigned char w2, unsigned char w3, unsigned char w4) {
	int valid = 1;
	
	if( w1 != 92 || w2 != 222 || w3 != 192 || w4 != 222 ) {
		valid = 0;
	}
	
	return valid;
}


/*
  DR Spectrometer header format
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
	unsigned int MAGIC1; 			// must always equal 0xC0DEC0DE
	unsigned long long timeTag0; 		// time tag of first frame in ``block''
	unsigned short int timeOffset;	// time offset reported by DP
	unsigned short int decFactor; 	// decimation factor
	unsigned int freqCode[2]; 		// DP frequency codes for each tuning
								//   indexing: 0..1 = Tuning 1..2
	unsigned int fills[4]; 			// fills for each pol/tuning combination
								//   indexing: 0..3 = X0, Y0 X1, Y1
	unsigned char errors[4]; 		// error flag for each pol/tuning combo
								//   indexing: 0..3 = X0, Y0 X1, Y1
	unsigned char beam;				// beam number
	unsigned int nFreqs;			// <Transform Length>
	unsigned int nInts;				// <Integration Count>
	unsigned int satCount[4];		// saturation count for each pol/tuning combo
								//   indexing: 0..3 = X0, Y0 X1, Y1
	unsigned int MAGIC2;			// must always equal 0xED0CED0C
} DRSpecHeader;
#pragma pack(pop)


/*
  TBW Reader (12-bit and 4-bit samples)
*/

static PyObject *readTBW(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data;
	int i;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Read in a single 1224 byte frame
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	unsigned char bytes[1224];
	i = fread(bytes, 1, sizeof(bytes), fh);	
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

	// Decode the header
	unsigned char sync1, sync2, sync3, sync4;
	sync1 = bytes[3];
	sync2 = bytes[2];
	sync3 = bytes[1];
	sync4 = bytes[0];
	if( !validSync(sync1, sync2, sync3, sync4) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	unsigned long int frameCount;
	frameCount = bytes[5]<<16 | bytes[6]<<8 | bytes[7];
	unsigned long int secondsCount;
	secondsCount = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
	unsigned short int tbwID;
	tbwID = bytes[12]<<8 | bytes[13];

	int bits = (tbwID>>14)&1;
	unsigned long long timeTag;
	timeTag = ((unsigned long long) bytes[16])<<56 | \
			((unsigned long long) bytes[17])<<48 | \
			((unsigned long long) bytes[18])<<40 | \
			((unsigned long long) bytes[19])<<32 | \
			((unsigned long long) bytes[20])<<24 | \
			((unsigned long long) bytes[21])<<16 | \
			((unsigned long long) bytes[22])<<8 | \
			bytes[23];

	// Create the output data array
	npy_intp dims[2];
	if(bits == 0) {
		// 12-bit Data -> 400 samples in the data section
		dims[0] = 2;
		dims[1] = 400;
		data = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_INT16);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		short int tempR;
		short int *a;
		a = (short int *) data->data;
		for(i=0; i<400; i++) {
			tempR = (bytes[24+3*i]<<4) | ((bytes[24+3*i+1]>>4)&15);
			tempR -= ((tempR&2048)<<1);
			*(a + i) = (short int) tempR;

			tempR = ((bytes[24+3*i+1]&15)<<8) | bytes[24+3*i+2];
			tempR -= ((tempR&2048)<<1);
			*(a + 400 + i) = (short int) tempR;
		}
	} else {
		// 4-bit Data -> 1200 samples in the data section
		dims[0] = 2;
		dims[1] = 1200;
		data = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_INT16);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		short int tempR;
		short int *a;
		a = (short int *) data->data;
		for(i=0; i<1200; i++) {
			tempR = (bytes[i+24]>>4)&15;
			tempR -= ((tempR&8)<<1);
			*(a + i) = (short int) tempR;

			tempR = bytes[i+24]&15;
			tempR -= ((tempR&8)<<1);
			*(a + 1200 + i) = (short int) tempR;
		}

	}

	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");

	temp = PyLong_FromUnsignedLong(frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);

	temp = PyLong_FromUnsignedLong(secondsCount);
	PyObject_SetAttrString(fHeader, "secondsCount", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("H", tbwID);
	PyObject_SetAttrString(fHeader, "tbwID", temp);
	Py_XDECREF(temp);

	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");

	temp = PyLong_FromUnsignedLongLong(timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);

	PyObject_SetAttrString(fData, "xy", PyArray_Return(data));

	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);

	output = Py_BuildValue("O", frame);
	return output;
}

PyDoc_STRVAR(readTBW_doc, \
"Function to read in a single TBW frame (header+data) and store the contents\n\
as a Frame object.  This function serves as a replacement for the pure python\n\
reader lsl.reader.tbw.readFrame.\n\
\n\
In order to use this reader in place of lsl.reader.tbw.readFrame change:\n\
\n\
\t>>> import lsl.reader.tbw as tbw\n\
\t>>> fh = open('some-tbw-file.dat', 'rb')\n\
\t>>> frame = tbw.readFrame(fh)\n\
\n\
to:\n\
\n\
\t>>> import lsl.reader.tbw as tbw\n\
\t>>> from lsl.reader._gofast import ReadTBW, syncError, eofError\n\
\t>>> fh = open('some-tbw-file.dat', 'rb')\n\
\t>>> frame = readTBW(fh, tbw.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.tbw`\n\
\tmodule.\n\
");


/*
  TBN Reader
*/

static PyObject *readTBN(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data;
	int i;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Read in a single 1048 byte frame
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	unsigned char bytes[1048];
	i = fread(bytes, 1, sizeof(bytes), fh);	
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

	// Decode the header
	unsigned char sync1, sync2, sync3, sync4;
	sync1 = bytes[3];
	sync2 = bytes[2];
	sync3 = bytes[1];
	sync4 = bytes[0];
	if( !validSync(sync1, sync2, sync3, sync4) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	unsigned long int frameCount;
	frameCount = bytes[5]<<16 | bytes[6]<<8 | bytes[7];
	unsigned long int tuningWord;
	tuningWord = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
	unsigned short int tbnID;
	tbnID = bytes[12]<<8 | bytes[13];
	unsigned short int gain;
	gain = bytes[14]<<8 | bytes[15];

	unsigned long long timeTag;
	timeTag = ((unsigned long long) bytes[16])<<56 | \
			((unsigned long long) bytes[17])<<48 | \
			((unsigned long long) bytes[18])<<40 | \
			((unsigned long long) bytes[19])<<32 | \
			((unsigned long long) bytes[20])<<24 | \
			((unsigned long long) bytes[21])<<16 | \
			((unsigned long long) bytes[22])<<8 | \
			bytes[23];

	// Create the output data array
	npy_intp dims[1];
	short int tempR, tempI;
	dims[0] = 512;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_COMPLEX64);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Fill the data array
	float complex *a;
	a = (float complex *) data->data;
	for(i=0; i<512; i++) {
		tempR = bytes[24+2*i];
		tempR -= ((tempR&128)<<1);
		tempI = bytes[24+2*i+1];
		tempI -= ((tempI&128)<<1);
		*(a + i) = (float) tempR + imaginary * (float) tempI;
	}

	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");

	temp = PyLong_FromUnsignedLong(frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);

	temp = PyLong_FromUnsignedLong(tuningWord);
	PyObject_SetAttrString(fHeader, "tuningWord", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("H", tbnID);
	PyObject_SetAttrString(fHeader, "tbnID", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("H", gain);
	PyObject_SetAttrString(fHeader, "gain", temp);
	Py_XDECREF(temp);

	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");

	temp = PyLong_FromUnsignedLongLong(timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);

	PyObject_SetAttrString(fData, "iq", PyArray_Return(data));

	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);

	output = Py_BuildValue("O", frame);
	return output;
}

PyDoc_STRVAR(readTBN_doc, \
"Function to read in a single TBN frame (header+data) and store the contents\n\
as a Frame object.  This function serves as a replacement for the pure python\n\
reader lsl.reader.tbn.readFrame.\n\
\n\
In order to use this reader in place of lsl.reader.tbn.readFrame change:\n\
\n\
\t>>> import lsl.reader.tbn as tbn\n\
\t>>> fh = open('some-tbn-file.dat', 'rb')\n\
\t>>> frame = tbn.readFrame(fh)\n\
\n\
to:\n\
\n\
\t>>> import lsl.reader.tbn as tbn\n\
\t>>> from lsl.reader._gofast import ReadTBN, syncError, eofError\n\
\t>>> fh = open('some-tbn-file.dat', 'rb')\n\
\t>> frame = readTBN(fh, tbn.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.tbn`\n\
\tmodule.\n\
");


/*
  DRX Reader
*/

static PyObject *readDRX(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data;
	int i;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Read in a single 4128 byte frame
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	unsigned char bytes[4128];
	i = fread(bytes, 1, sizeof(bytes), fh);	
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

	// Decode the header
	unsigned char sync1, sync2, sync3, sync4;
	sync1 = bytes[3];
	sync2 = bytes[2];
	sync3 = bytes[1];
	sync4 = bytes[0];
	if( !validSync(sync1, sync2, sync3, sync4) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	unsigned char drxID;
	drxID = bytes[4];
	unsigned long int frameCount;
	frameCount = bytes[5]<<16 | bytes[6]<<8 | bytes[7];
	unsigned long int secondsCount;
	secondsCount  = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
	unsigned short int decimation;
	decimation = bytes[12]<<8 | bytes[13];
	unsigned short int timeOffset;
	timeOffset = bytes[14]<<8 | bytes[15];

	unsigned long long timeTag;
	timeTag = ((unsigned long long) bytes[16])<<56 | \
			((unsigned long long) bytes[17])<<48 | \
			((unsigned long long) bytes[18])<<40 | \
			((unsigned long long) bytes[19])<<32 | \
			((unsigned long long) bytes[20])<<24 | \
			((unsigned long long) bytes[21])<<16 | \
			((unsigned long long) bytes[22])<<8 | \
			bytes[23];
	unsigned long tuningWord;
	tuningWord = ((unsigned long) bytes[24])<<24 | \
			   ((unsigned long) bytes[25])<<16 | \
			   ((unsigned long) bytes[26])<<8 | \
			   bytes[27];
	unsigned long flags;
	flags  =  ((unsigned long) bytes[28])<<24 | \
			((unsigned long) bytes[29])<<16 | \
			((unsigned long) bytes[30])<<8 | \
			bytes[31];
	
	// Create the output data array
	npy_intp dims[1];
	dims[0] = 4096;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_COMPLEX64);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}

	// Fill the data array
	short int tempR, tempI;
	float complex *a;
	a = (float complex *) data->data;
	for(i=0; i<4096; i++) {
		tempR = (bytes[i+32]>>4)&15;
		tempR -= ((tempR&8)<<1);
		tempI = bytes[i+32]&15;
		tempI -= ((tempI&8)<<1);
		*(a + i) = (float) tempR + imaginary * (float) tempI;
	}

	// Save the data to the frame object
	// 1. Header
	fHeader = PyObject_GetAttrString(frame, "header");

	temp = PyLong_FromUnsignedLong(frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("B", drxID);
	PyObject_SetAttrString(fHeader, "drxID", temp);
	Py_XDECREF(temp);

	temp = PyLong_FromUnsignedLong(secondsCount);
	PyObject_SetAttrString(fHeader, "secondsCount", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("H", decimation);
	PyObject_SetAttrString(fHeader, "decimation", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", timeOffset);
	PyObject_SetAttrString(fHeader, "timeOffset", temp);
	Py_XDECREF(temp);

	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");

	temp = PyLong_FromUnsignedLongLong(timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);

	temp = PyLong_FromUnsignedLong(tuningWord);
	PyObject_SetAttrString(fData, "tuningWord", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(flags);
	PyObject_SetAttrString(fData, "flags", temp);
	Py_XDECREF(temp);

	PyObject_SetAttrString(fData, "iq", PyArray_Return(data));

	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);

	output = Py_BuildValue("O", frame);
	return output;
}

PyDoc_STRVAR(readDRX_doc, \
"Function to read in a single DRX frame (header+data) and store the contents\n\
as a Frame object.  This function serves as a replacement for the pure python\n\
reader lsl.reader.drx.readFrame.\n\
\n\
In order to use this reader in place of lsl.reader.drx.readFrame change:\n\
\n\
\t>>> import lsl.reader.tbn as drx\n\
\t>>> fh = open('some-drx-file.dat', 'rb')\n\
\t>>> frame = drx.readFrame(fh)\n\
\n\
to:\n\
\n\
\t>>> import lsl.reader.drx as drx\n\
\t>>> from lsl.reader._gofast import ReadDRX, syncError, eofError\n\
\t>>> fh = open('some-drx-file.dat', 'rb')\n\
\t>>> frame = readDRX(fh, tbn.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.drx`\n\
\tmodule.\n\
");


/*
  DR Spectrometer Reader
*/

static PyObject *readDRSpec(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyObject *tuningWords, *fills, *errors, *saturations;
	PyArrayObject *dataX0, *dataY0, *dataX1, *dataY1;
	int i;
	
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
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(syncError, "Sync word differs from expected");
		return NULL;
	}
	
	// Read in the data section
	float XY[4*header.nFreqs];
	i = fread(XY, sizeof(float), 4*header.nFreqs, fh);
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
	
	// Create the output data arrays
	npy_intp dims[1];
	dims[0] = header.nFreqs;
	dataX0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataX0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - X0");
		Py_XDECREF(dataX0);
		return NULL;
	}
	dataY0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataY0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - Y0");
		Py_XDECREF(dataX0);
		Py_XDECREF(dataY0);
		return NULL;
	}
	dataX1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataX1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - X1");
		Py_XDECREF(dataX0);
		Py_XDECREF(dataY0);
		Py_XDECREF(dataX1);
		return NULL;
	}
	dataY1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataY1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - Y1");
		Py_XDECREF(dataX0);
		Py_XDECREF(dataY0);
		Py_XDECREF(dataX1);
		Py_XDECREF(dataY1);
		return NULL;
	}

	// Fill the data arrays
	float *a, *b, *c, *d;
	a = (float *) dataX0->data;
	b = (float *) dataY0->data;
	c = (float *) dataX1->data;
	d = (float *) dataY1->data;
	for(i=0; i<header.nFreqs; i++) {
		*(a + i) = XY[0*header.nFreqs + i];
		*(b + i) = XY[1*header.nFreqs + i];
		*(c + i) = XY[2*header.nFreqs + i];
		*(d + i) = XY[3*header.nFreqs + i];
	}
	
	// Fill in the multi-value fields (tuningWords, fills, errors)
	tuningWords = PyList_New(2);
	if(tuningWords == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - tuningWords");
		Py_XDECREF(dataX0);
		Py_XDECREF(dataY0);
		Py_XDECREF(dataX1);
		Py_XDECREF(dataY1);
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
		Py_XDECREF(dataX0);
		Py_XDECREF(dataY0);
		Py_XDECREF(dataX1);
		Py_XDECREF(dataY1);
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
		Py_XDECREF(dataX0);
		Py_XDECREF(dataY0);
		Py_XDECREF(dataX1);
		Py_XDECREF(dataY1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		Py_XDECREF(errors);
		return NULL;
	}
	for(i=0; i<4; i++) {
		temp = Py_BuildValue("H", header.errors[i]);
		PyList_SetItem(errors, i, temp);
	}
	
	saturations = PyList_New(4);
	if(saturations == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - saturations");
		Py_XDECREF(dataX0);
		Py_XDECREF(dataY0);
		Py_XDECREF(dataX1);
		Py_XDECREF(dataY1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		Py_XDECREF(errors);
		Py_XDECREF(saturations);
		return NULL;
	}
	for(i=0; i<4; i++) {
		temp = Py_BuildValue("H", header.satCount[i]);
		PyList_SetItem(saturations, i, temp);
	}
	
	// Save the data to the frame object
	// 1. Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = Py_BuildValue("H", header.beam);
	PyObject_SetAttrString(fHeader, "beam", temp);
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
	
	PyObject_SetAttrString(fData, "flags", errors);
	
	PyObject_SetAttrString(fData, "saturations", saturations);

	PyObject_SetAttrString(fData, "X0", PyArray_Return(dataX0));
	PyObject_SetAttrString(fData, "Y0", PyArray_Return(dataY0));
	PyObject_SetAttrString(fData, "X1", PyArray_Return(dataX1));
	PyObject_SetAttrString(fData, "Y1", PyArray_Return(dataY1));

	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(tuningWords);
	Py_XDECREF(fills);
	Py_XDECREF(errors);
	Py_XDECREF(saturations);
	Py_XDECREF(fData);
	Py_XDECREF(dataX0);
	Py_XDECREF(dataY0);
	Py_XDECREF(dataX1);
	Py_XDECREF(dataY1);

	output = Py_BuildValue("O", frame);
	return output;
}

PyDoc_STRVAR(readDRSpec_doc, \
"Function to read in a DR spectrometer header structure and data and retrun\n\
a drspec.Frame instance.");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef GoFastMethods[] = {
	{"readTBW", readTBW, METH_VARARGS, readTBW_doc}, 
	{"readTBN", readTBN, METH_VARARGS, readTBN_doc}, 
	{"readDRX", readDRX, METH_VARARGS, readDRX_doc}, 
	{"readDRSpec", readDRSpec, METH_VARARGS, readDRSpec_doc},
	{NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(GoFast_doc, "Go Fast! (TM) - TBW, TBN, and DRX readers written in C");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_gofast(void) {
	PyObject *m, *dict1, *dict2;

	// Module definitions and functions
	m = Py_InitModule3("_gofast", GoFastMethods, GoFast_doc);
	import_array();

	// Exceptions
	
	//   1.  syncError -> similar to lsl.reader.errors.syncError
	dict1 = (PyObject *) PyDict_New();
	if(dict1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
		Py_XDECREF(dict1);
		Py_XDECREF(m);
		return NULL;
	}
	PyDict_SetItemString(dict1, "__doc__", \
		PyString_FromString("Exception raised when a reader encounters an error with one or more of the four sync. words."));
	syncError = PyErr_NewException("_gofast.syncError", PyExc_IOError, dict1);
	Py_INCREF(syncError);
	PyModule_AddObject(m, "syncError", syncError);
	
	//    2. eofError -> similar to lsl.reader.errors.eofError
	dict2 = (PyObject *) PyDict_New();
	if(dict2 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
		Py_XDECREF(dict1);
		Py_XDECREF(syncError);
		Py_XDECREF(dict2);
		Py_XDECREF(m);
		return NULL;
	}
	PyDict_SetItemString(dict2, "__doc__", \
		PyString_FromString("Exception raised when a reader encounters the end-of-file while reading."));
	eofError = PyErr_NewException("_gofast.eofError", PyExc_IOError, dict2);
	Py_INCREF(eofError);
	PyModule_AddObject(m, "eofError", eofError);
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.4"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
}
