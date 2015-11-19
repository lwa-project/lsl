#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"


/*
  Exceptions for the Go Fast! Readers
*/

PyObject *syncError;
PyObject *eofError;


/*
  Validate a collection of Mark 5C sync words.  Return 1 is all are
  valid.
*/

int validSync5C(unsigned int syncWord) {
	int valid = 0;
	
	if( syncWord == 0x5CDEC0DE ) {
		valid = 1;
	}
	
	return valid;
}


/* 
  Look-up Tables
*/
short int tbw4LUT[256][2];
float tbnLUT[256];
float drxLUT[256][2];
float tbfLUT[256][2];
float drx8LUT[256];

static void initLWALUTs(void) {
	// Look-up table inialization function from the VDIFIO library
	
	int i,j;
	short int t;
	
	//TBW - 4 bit
	for(i=0; i<256; i++) {
		for(j=0; j<2; j++) {
			t = (i >> 4*(1-j)) & 15;
			tbw4LUT[i][j] = t;
			tbw4LUT[i][j] -= ((t&8)<<1);
		}
	}
	
	// TBN & DRX8
	for(i=0; i<256; i++) {
		tbnLUT[i] = i;
		tbnLUT[i] -= ((i&128)<<1);
		
		drx8LUT[i] = i;
		drx8LUT[i] -= ((i&128)<<1);
	}
	
	// DRX & TBF
	for(i=0; i<256; i++) {
		for(j=0; j<2; j++) {
			t = (i >> 4*(1-j)) & 15;
			drxLUT[i][j] = t;
			drxLUT[i][j] -= ((t&8)<<1);
			
			tbfLUT[i][j] = t;
			tbfLUT[i][j] -= ((t&8)<<1);
		}
	}
}


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef GoFastMethods[] = {
	{"readTBW",    (PyCFunction) readTBW,    METH_VARARGS,               readTBW_doc   }, 
	{"readTBN",    (PyCFunction) readTBN,    METH_VARARGS,               readTBN_doc   }, 
	{"readDRX",    (PyCFunction) readDRX,    METH_VARARGS,               readDRX_doc   }, 
	{"readDRSpec", (PyCFunction) readDRSpec, METH_VARARGS,               readDRSpec_doc},
	{"readVDIF",   (PyCFunction) readVDIF,   METH_VARARGS|METH_KEYWORDS, readVDIF_doc  }, 
	{"readTBF",    (PyCFunction) readTBF,    METH_VARARGS,               readTBF_doc   }, 
	{"readCOR",    (PyCFunction) readCOR,    METH_VARARGS,               readCOR_doc   }, 
	{"readDRX8",   (PyCFunction) readDRX8,   METH_VARARGS,               readDRX8_doc  }, 
	{NULL,         NULL,                     0,                          NULL          }
};

PyDoc_STRVAR(GoFast_doc, "Go Fast! (TM) - TBW, TBN, DRX, DR Spectrometer, and VDIF readers written in C");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_gofast(void) {
	PyObject *m, *dict1, *dict2;
	
	// Initialize the look-up tables
	initLWALUTs();
	initVDIFLUTs();
	
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
	}
	PyDict_SetItemString(dict2, "__doc__", \
		PyString_FromString("Exception raised when a reader encounters the end-of-file while reading."));
	eofError = PyErr_NewException("_gofast.eofError", PyExc_IOError, dict2);
	Py_INCREF(eofError);
	PyModule_AddObject(m, "eofError", eofError);
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.8"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
}
