#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#if defined(__linux__)
/* Linux */
#include <byteswap.h>
#elif defined(__APPLE__) && defined(__MACH__)
/* OSX */
#include <libkern/OSByteOrder.h>
#define __bswap_16 OSSwapInt16
#define __bswap_32 OSSwapInt32
#define __bswap_64 OSSwapInt64
#endif

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"

/*
  TBW Reader (12-bit and 4-bit samples)
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
	unsigned int syncWord;
	union {
		struct {
			unsigned int frameCount:24;
			unsigned char id:8;
		};
		unsigned int frameCountWord;
	};
	unsigned int secondsCount;
	union {
		unsigned short int tbwID;
		struct {
			unsigned short int stand:14;
			unsigned char bits:1;
			unsigned char isTBW:1;
		};
	};
	unsigned short int unassigned;
} TBWHeader;


typedef struct {
	unsigned long long timeTag;
	unsigned char bytes[1200];
} TBWPayload;


typedef struct {
	TBWHeader header;
	TBWPayload data;
} TBWFrame;
#pragma pack(pop)


PyObject *readTBW(PyObject *self, PyObject *args) {
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
	TBWFrame cFrame;
	i = fread(&cFrame, sizeof(cFrame), 1, fh);
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

	// Validate
	if( !validSync5C(cFrame.header.syncWord) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	
	// Swap the bits around
	cFrame.header.frameCountWord = __bswap_32(cFrame.header.frameCountWord);
	cFrame.header.secondsCount = __bswap_32(cFrame.header.secondsCount);
	cFrame.header.tbwID = __bswap_16(cFrame.header.tbwID);
	cFrame.data.timeTag = __bswap_64(cFrame.data.timeTag);
	
	// Create the output data array
	npy_intp dims[2];
	if(cFrame.header.bits == 0) {
		// 12-bit Data -> 400 samples in the data section
		dims[0] = (npy_intp) 2;
		dims[1] = (npy_intp) 400;
		data = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_INT16, 0);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		short int tempR;
		short int *a;
		a = (short int *) PyArray_DATA(data);
		for(i=0; i<400; i++) {
			tempR = (cFrame.data.bytes[3*i]<<4) | ((cFrame.data.bytes[3*i+1]>>4)&15);
			tempR -= ((tempR&2048)<<1);
			*(a + i) = (short int) tempR;

			tempR = ((cFrame.data.bytes[3*i+1]&15)<<8) | cFrame.data.bytes[3*i+2];
			tempR -= ((tempR&2048)<<1);
			*(a + 400 + i) = (short int) tempR;
		}
	} else {
		// 4-bit Data -> 1200 samples in the data section
		dims[0] = (npy_intp) 2;
		dims[1] = (npy_intp) 1200;
		data = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_INT16, 0);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		const short int *fp;
		short int *a;
		a = (short int *) PyArray_DATA(data);
		for(i=0; i<1200; i++) {
			fp = tbw4LUT[ cFrame.data.bytes[i] ];
			*(a + i) = (short int) fp[0];
			*(a + 1200 + i) = (short int) fp[1];
		}
		
	}
	
	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = PyLong_FromUnsignedLong(cFrame.header.frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.header.secondsCount);
	PyObject_SetAttrString(fHeader, "secondsCount", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", cFrame.header.tbwID);
	PyObject_SetAttrString(fHeader, "tbwID", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(cFrame.data.timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);
	
	PyObject_SetAttrString(fData, "xy", PyArray_Return(data));
	
	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	output = Py_BuildValue("O", frame);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);
	
	return output;
}

char readTBW_doc[] = PyDoc_STR(\
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
