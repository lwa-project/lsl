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
  TBW Reader (12-bit and 4-bit samples)
*/

static PyObject *readTBW(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData;
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
	short int temp;
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
		short int *a;
		a = (short int *) data->data;
		for(i=0; i<400; i++) {
			temp = (bytes[24+3*i]<<4) | ((bytes[24+3*i+1]>>4)&15);
			temp -= ((temp&2048)<<1);
			*(a + i) = (short int) temp;

			temp = ((bytes[24+3*i+1]&15)<<8) | bytes[24+3*i+2];
			temp -= ((temp&2048)<<1);
			*(a + 400 + i) = (short int) temp;
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
		short int *a;
		a = (short int *) data->data;
		for(i=0; i<1200; i++) {
			temp = (bytes[i+24]>>4)&15;
			temp -= ((temp&8)<<1);
			*(a + i) = (short int) temp;

			temp = bytes[i+24]&15;
			temp -= ((temp&8)<<1);
			*(a + 400 + i) = (short int) temp;
		}

	}

	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");
	PyObject_SetAttrString(fHeader, "frameCount", Py_BuildValue("l", frameCount));
	PyObject_SetAttrString(fHeader, "secondsCount", Py_BuildValue("l", secondsCount));
	PyObject_SetAttrString(fHeader, "tbwID", Py_BuildValue("i", tbwID));
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	PyObject_SetAttrString(fData, "timeTag", Py_BuildValue("l", timeTag));
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
In order to use this\n\
reader in place of lsl.reader.tbw.readFrame change:\n\
  >>> import lsl.reader.tbw as tbw\n\
  >>> fh = open('some-tbw-file.dat', 'rb')\n\
  >>> frame = tbw.readFrame(fh)\n\
to:\n\
  >>> import lsl.reader.tbw as tbw\n\
  >>> from lsl.reader._gofast import ReadTBW, syncError, eofError\n\
  >>> fh = open('some-tbw-file.dat', 'rb')\n\
  >>> frame = readTBW(fh, tbw.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classs of IOError.\n\
");


/*
  TBN Reader
*/

static PyObject *readTBN(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData;
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
	unsigned long int secondsCount;
	secondsCount = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
	unsigned short int tbnID;
	tbnID = bytes[12]<<8 | bytes[13];

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
	npy_intp *fLoc;
	fLoc = PyDimMem_NEW(1);
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
		fLoc[0] = (npy_intp) i;
		tempR = bytes[24+2*i];
		tempR -= ((tempR&128)<<1);
		tempI = bytes[24+2*i+1];
		tempI -= ((tempI&128)<<1);
		*(a + i) = (float) tempR + imaginary * (float) tempI;
	}
	PyDimMem_FREE(fLoc);

	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");
	PyObject_SetAttrString(fHeader, "frameCount", Py_BuildValue("l", frameCount));
	PyObject_SetAttrString(fHeader, "secondsCount", Py_BuildValue("l", secondsCount));
	PyObject_SetAttrString(fHeader, "tbnID", Py_BuildValue("i", tbnID));
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	PyObject_SetAttrString(fData, "timeTag", Py_BuildValue("l", timeTag));
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
In order to use this\n\
reader in place of lsl.reader.tbn.readFrame change:\n\
  >>> import lsl.reader.tbn as tbn\n\
  >>> fh = open('some-tbn-file.dat', 'rb')\n\
  >>> frame = tbn.readFrame(fh)\n\
to:\n\
  >>> import lsl.reader.tbn as tbn\n\
  >>> from lsl.reader._gofast import ReadTBN, syncError, eofError\n\
  >>> fh = open('some-tbn-file.dat', 'rb')\n\
  >>> frame = readTBN(fh, tbn.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classs of IOError.\n\
");


/*
  DRX Reader
*/

static PyObject *readDRX(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData;
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
	timeTag  = ((unsigned long long) bytes[16])<<56;
	timeTag |= ((unsigned long long) bytes[17])<<48; 
	timeTag |= ((unsigned long long) bytes[18])<<40;
	timeTag |= ((unsigned long long) bytes[19])<<32;
	timeTag |= ((unsigned long long) bytes[20])<<24;
	timeTag |= ((unsigned long long) bytes[21])<<16;
	timeTag |= ((unsigned long long) bytes[22])<<8;
	timeTag |= bytes[23];
	unsigned long long flags;
	flags  = ((unsigned long long) bytes[24])<<56;
	flags |= ((unsigned long long) bytes[25])<<48;
	flags |= ((unsigned long long) bytes[26])<<40;
	flags |= ((unsigned long long) bytes[27])<<32;
	flags |= ((unsigned long long) bytes[28])<<24;
	flags |= ((unsigned long long) bytes[29])<<16;
	flags |= ((unsigned long long) bytes[30])<<8;
	flags |= bytes[31];
	
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
	npy_intp *fLoc;
	fLoc = PyDimMem_NEW(1);
	short int tempR, tempI;
	float complex *a;
	a = (float complex *) data->data;
	for(i=0; i<4096; i++) {
		fLoc[0] = (npy_intp) i;
		tempR = (bytes[i+32]>>4)&15;
		tempR -= ((tempR&8)<<1);
		tempI = bytes[i+32]&15;
		tempI -= ((tempI&8)<<1);
		*(a + i) = (float) tempR + imaginary * (float) tempI;
	}
	PyDimMem_FREE(fLoc);

	// Save the data to the frame object
	// 1. Header
	fHeader = PyObject_GetAttrString(frame, "header");
	PyObject_SetAttrString(fHeader, "frameCount", Py_BuildValue("l", frameCount));
	PyObject_SetAttrString(fHeader, "drxID", Py_BuildValue("i", drxID));
	PyObject_SetAttrString(fHeader, "secondsCount", Py_BuildValue("l", secondsCount));
	PyObject_SetAttrString(fHeader, "decimation", Py_BuildValue("i", decimation));
	PyObject_SetAttrString(fHeader, "timeOffset", Py_BuildValue("i", timeOffset));
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	PyObject_SetAttrString(fData, "timeTag", Py_BuildValue("l", timeTag));
	PyObject_SetAttrString(fData, "flags", Py_BuildValue("l", flags));
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
In order to use this\n\
reader in place of lsl.reader.drx.readFrame change:\n\
  >>> import lsl.reader.tbn as drx\n\
  >>> fh = open('some-drx-file.dat', 'rb')\n\
  >>> frame = drx.readFrame(fh)\n\
to:\n\
  >>> import lsl.reader.drx as drx\n\
  >>> from lsl.reader._gofast import ReadDRX, syncError, eofError\n\
  >>> fh = open('some-drx-file.dat', 'rb')\n\
  >>> frame = readDRX(fh, tbn.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classs of IOError.\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef GoFastMethods[] = {
	{"readTBW", readTBW, METH_VARARGS, readTBW_doc}, 
	{"readTBN", readTBN, METH_VARARGS, readTBN_doc}, 
	{"readDRX", readDRX, METH_VARARGS, readDRX_doc}, 
	{NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(GoFast_doc, "Go Fast! (TM) - TBW, TBN, and DRX readers written in C");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_gofast(void) {
	PyObject *m;

	// Module definitions and functions
	m = Py_InitModule3("_gofast", GoFastMethods, GoFast_doc);
	import_array();

	// Exceptions
	//   1.  syncError -> similar to lsl.reader.errors.syncError
	syncError = PyErr_NewException("_gofast.syncError", PyExc_IOError, NULL);
	Py_INCREF(syncError);
	PyModule_AddObject(m, "syncError", syncError);
	//    2. eofError -> similar to lsl.reader.errors.eofError
	eofError = PyErr_NewException("_gofast.eofError", PyExc_IOError, NULL);
	Py_INCREF(eofError);
	PyModule_AddObject(m, "eofError", eofError);
}
