#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "numpy/arrayobject.h"

#define PI 3.1415926535898
#define imaginary _Complex_I

static PyObject *syncError;
static PyObject *eofError;

/*
  Validate a collection of Mark 5C sync words.  Return 1 is all are
  valid.
*/

int validSync(unsigned char w1, unsigned char w2, unsigned char w3, unsigned char w4) {
	int valid = 1;
	
	if( w1 != 92 || w2 != 222 || w3 != 192 || w4 != 222 ) {
		valid = 1;
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
		fclose(fh);
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_SetString(eofError, "End of file encountered during filehandle read");
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
		PyErr_SetString(syncError, "Mark 5C sync word differs from expected");
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
	npy_intp *fLoc;
	fLoc = PyDimMem_NEW(2);
	short int temp;
	if(bits == 0) {
		dims[0] = 2;
		dims[1] = 400;
		data = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_INT16);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		for(i=0; i<400; i++) {
			fLoc[0] = (npy_intp) 0;
			fLoc[1] = (npy_intp) i;
			temp = bytes[24+3*i]<<4 | (bytes[24+3*i+1]>>4)&15;
			temp -= 4096*((temp>>11)&1);
			*(short int *) PyArray_GetPtr(data, fLoc) = temp;

			fLoc[0] = (npy_intp) 1;
			fLoc[1] = (npy_intp) i;
			temp = (bytes[24+3*i+1]&15)<<8 | bytes[24+3*i+2];
			temp -= 4096*((temp>>11)&1);
			*(short int *) PyArray_GetPtr(data, fLoc) = temp;
		}
	} else {
		dims[0] = 2;
		dims[1] = 1200;
		data = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_UINT16);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		for(i=0; i<1200; i++) {
			fLoc[0] = 0;
			fLoc[1] = (npy_intp) i;
			temp = (bytes[i+24]>>4)&15;
			temp -= 16*((temp>>3)&1);
			*(short int *) PyArray_GetPtr(data, fLoc) = temp;

			fLoc[0] = 1;
			fLoc[1] = (npy_intp) i;
			temp = bytes[i+24]&15;
			temp -= 16*((temp>>3)&1);
			*(short int *) PyArray_GetPtr(data, fLoc) = temp;
		}

	}
	PyDimMem_FREE(fLoc);

	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttr(frame, Py_BuildValue("s", "header"));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "frameCount"), Py_BuildValue("l", frameCount));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "secondsCount"), Py_BuildValue("l", secondsCount));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "tbwID"), Py_BuildValue("i", tbwID));
	// 2. Data
	fData = PyObject_GetAttr(frame, Py_BuildValue("s", "data"));
	PyObject_SetAttr(fData, Py_BuildValue("s", "timeTag"), Py_BuildValue("l", timeTag));
	PyObject_SetAttr(fData, Py_BuildValue("s", "xy"), PyArray_Return(data));
	// 3. Frame
	PyObject_SetAttr(frame, Py_BuildValue("s", "header"), fHeader);
	PyObject_SetAttr(frame, Py_BuildValue("s", "data"), fData);
	Py_DECREF(fHeader);
	Py_DECREF(fData);
	Py_XDECREF(data);

	return Py_BuildValue("O", frame);
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
		fclose(fh);
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_SetString(eofError, "End of file encountered during filehandle read");
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
		PyErr_SetString(syncError, "Mark 5C sync word differs from expected");
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
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, PyArray_CDOUBLE);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Fill the data array
	for(i=0; i<512; i++) {
		fLoc[0] = (npy_intp) i;
		tempR = bytes[24+2*i];
		tempR -= 256*((tempR>>7)&1);
		tempI = bytes[24+2*i+1];
		tempI -= 256*((tempI>>7)&1);
		*(double complex *) PyArray_GetPtr(data, fLoc) = (double) tempR + imaginary * (double) tempI; 
	}
	PyDimMem_FREE(fLoc);

	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttr(frame, Py_BuildValue("s", "header"));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "frameCount"), Py_BuildValue("l", frameCount));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "secondsCount"), Py_BuildValue("l", secondsCount));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "tbnID"), Py_BuildValue("i", tbnID));
	// 2. Data
	fData = PyObject_GetAttr(frame, Py_BuildValue("s", "data"));
	PyObject_SetAttr(fData, Py_BuildValue("s", "timeTag"), Py_BuildValue("l", timeTag));
	PyObject_SetAttr(fData, Py_BuildValue("s", "iq"), PyArray_Return(data));
	// 3. Frame
	PyObject_SetAttr(frame, Py_BuildValue("s", "header"), fHeader);
	PyObject_SetAttr(frame, Py_BuildValue("s", "data"), fData);
	Py_DECREF(fHeader);
	Py_DECREF(fData);
	Py_XDECREF(data);

	return Py_BuildValue("O", frame);
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
		fclose(fh);
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_SetString(eofError, "End of file encountered during filehandle read");
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
		PyErr_SetString(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	unsigned char drxID;
	drxID = bytes[4];
	unsigned long int frameCount;
	frameCount = bytes[5]<<16 | bytes[6]<<8 | bytes[7];
	unsigned long int secondsCount;
	secondsCount = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
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
	unsigned long long flags;
	flags = ((unsigned long long) bytes[24])<<56 | \
			((unsigned long long) bytes[25])<<48 | \
			((unsigned long long) bytes[26])<<40 | \
			((unsigned long long) bytes[27])<<32 | \
			((unsigned long long) bytes[28])<<24 | \
			((unsigned long long) bytes[29])<<16 | \
			((unsigned long long) bytes[30])<<8 | \
			bytes[31];
	// Create the output data array
	npy_intp dims[2];
	npy_intp *fLoc;
	fLoc = PyDimMem_NEW(1);
	
	short int tempR, tempI;
	dims[0] = 4096;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, PyArray_CDOUBLE);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}

	// Fill the data array
	for(i=0; i<4096; i++) {
		fLoc[0] = (npy_intp) i;
		tempR = (bytes[i+32]>>4)&15;
		tempR -= 16*((tempR>>3)&1);
		tempI = bytes[i+24]&15;
		tempI -= 16*((tempI>>3)&1);
		*(double complex *) PyArray_GetPtr(data, fLoc) = (double) tempR + imaginary * (double) tempI;
	}
	PyDimMem_FREE(fLoc);

	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttr(frame, Py_BuildValue("s", "header"));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "frameCount"), Py_BuildValue("l", frameCount));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "drxID"), Py_BuildValue("i", drxID));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "secondsCount"), Py_BuildValue("l", secondsCount));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "decimation"), Py_BuildValue("i", decimation));
	PyObject_SetAttr(fHeader, Py_BuildValue("s", "timeOffset"), Py_BuildValue("i", timeOffset));
	// 2. Data
	fData = PyObject_GetAttr(frame, Py_BuildValue("s", "data"));
	PyObject_SetAttr(fData, Py_BuildValue("s", "timeTag"), Py_BuildValue("l", timeTag));
	PyObject_SetAttr(fData, Py_BuildValue("s", "flags"), Py_BuildValue("l", flags));
	PyObject_SetAttr(fData, Py_BuildValue("s", "iq"), PyArray_Return(data));
	// 3. Frame
	PyObject_SetAttr(frame, Py_BuildValue("s", "header"), fHeader);
	PyObject_SetAttr(frame, Py_BuildValue("s", "data"), fData);
	Py_DECREF(fHeader);
	Py_DECREF(fData);
	Py_XDECREF(data);

	return Py_BuildValue("O", frame);
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
