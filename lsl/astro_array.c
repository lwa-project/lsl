/* 
  astro_array.c
  
  Array-based astronomical calculations implemented in C using numpy and libnova.

  $Id: astro_array.c 85 2010-05-19 22:53:54Z dwood $ 
*/


#include "Python.h"
#include "numpy/arrayobject.h"

#include <libnova/transform.h>
#include <libnova/julian_day.h>
#include <libnova/utility.h>
#include <libnova/sidereal_time.h>
#include <math.h>
#include <stdio.h>



/*  function hrz_from_equ() */
/*  Python/numpy wrapper to compute az/el from ra/dec for an array of points using libnova */
/*    This is much faster than iterating in python */

/* Written by Johnathan York, 5/2006, University of Texas Applied Research Laboratories */


static int 
hrz_from_equ(int n,double lng, double lat, double JD, double *ra, double *dec, double *alt, double *az)
{
        struct ln_equ_posn object;
        struct ln_hrz_posn hrz;
        struct ln_lnlat_posn observer;
        int i;

        observer.lng=lng;
        observer.lat=lat;

        for(i=0;i<n;i++) {
          object.ra=*(ra++);
          object.dec=*(dec++);

           ln_get_hrz_from_equ (&object, &observer, JD, &hrz);

		 /* Note:  libnova works with an azimuth definition that 0 is south, not 
		    the standard definition that 0 is north.  This convention difference 
		    was not noted in the pervious version of astro_array.c but this should
		    take care of it.
		 */
           *(alt++)=hrz.alt;
           *(az++)=fmod((hrz.az+180.0), 360.0);
        }

 return n;
}


static PyObject *
Py_hrz_from_equ(PyObject *obj, PyObject *args)
{
        double lng,lat,JD;
        PyObject   *ora, *odec, *ret;
        PyArrayObject *ra = NULL, *dec = NULL, *alt = NULL, *az = NULL;

        if (!PyArg_ParseTuple(args, "dddOO", &lng, &lat, &JD,&ora, &odec)) {
                PyErr_Format(PyExc_TypeError, "Invalid parameters.");
                goto _fail;
        }

       /* Align, Byteswap, Contiguous, Typeconvert */
        ra  = (PyArrayObject*)  PyArray_ContiguousFromObject(ora, PyArray_DOUBLE,1,1);
        dec  = (PyArrayObject*) PyArray_ContiguousFromObject(odec, PyArray_DOUBLE,1,1);

        if (!ra || !dec) {
                PyErr_Format(PyExc_TypeError, "error converting array inputs.");
                goto _fail;
        }

        if ((ra->nd != 1) || (dec->nd != 1)) {
                PyErr_Format(PyExc_TypeError,"arrays must have 1 dimension.");
                goto _fail;
        }

        if (ra->dimensions[0]!=dec->dimensions[0]) {
                PyErr_Format(PyExc_TypeError,"input data arrays need identitcal shapes.");
                goto _fail;
        }

        alt = (PyArrayObject*) PyArray_SimpleNew(1,ra->dimensions,PyArray_DOUBLE);
        az  = (PyArrayObject*) PyArray_SimpleNew(1,ra->dimensions,PyArray_DOUBLE);

        if (!alt || !az) {
                PyErr_Format(PyExc_TypeError, "error creating array outputs.");
                goto _fail;
        }

        hrz_from_equ(ra->dimensions[0],lng,lat,JD,(double*)ra->data,(double*)dec->data,
                   (double*)alt->data, (double*)az->data);

        Py_XDECREF(ra);
        Py_XDECREF(dec);


        /* Align, Byteswap, Contiguous, Typeconvert */
        ret=Py_BuildValue("OO",PyArray_Return(alt),PyArray_Return(az));

        Py_XDECREF(alt);
        Py_XDECREF(az);

        return ret;

_fail:
        Py_XDECREF(ra);
        Py_XDECREF(dec);
        Py_XDECREF(alt);
        Py_XDECREF(az);
        return NULL;
}


PyDoc_STRVAR(hrz_from_equ_doc,
    "Compute alt/az from ra/dec\nhrz_from_equ(lng, lat, jD, ra, dec) -> (alt,az)");


static PyMethodDef AstroArrayMethods[] = {
    {"hrz_from_equ",  Py_hrz_from_equ, METH_VARARGS, hrz_from_equ_doc},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyDoc_STRVAR(astro_array_doc,
    "Array-based astronomical calculations implemented in C\nusing numpy and libnova");

PyMODINIT_FUNC
initastro_array(void)
{
  PyObject *m;

  // Module definitions and functions
  m = Py_InitModule3("astro_array", AstroArrayMethods, astro_array_doc);
  import_array();
  
  // Version and revision information
  PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
  PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
  
  // Author and maintainer information
  PyModule_AddObject(m, "__author__", PyString_FromString("J. York"));
  PyModule_AddObject(m, "__maintainer__", PyString_FromString("Jayce Dowell"));
  
}
