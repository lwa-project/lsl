#include "complex_int8.h"

typedef struct {
    PyObject_HEAD
    complex_int8 obval;
} PyComplexInt8;

static PyTypeObject PyComplexInt8_Type;

PyArray_Descr *complex_int8_descr;

static NPY_INLINE int PyComplexInt8_Check(PyObject *o) {
    return PyObject_IsInstance(o, (PyObject*) &PyComplexInt8_Type);
}

static PyObject* PyComplexInt8_FromComplexInt8(complex_int8 c) {
    PyComplexInt8 *p = (PyComplexInt8 *) PyComplexInt8_Type.tp_alloc(&PyComplexInt8_Type, 0);
    if( p != NULL ) {
        p->obval = c;
    }
    return (PyObject *) p;
}

#define PyComplexInt8_AsComplexInt8(c, o)                      \
    if( PyComplexInt8_Check(o) ) {                             \
        c = ((PyComplexInt8 *) o)->obval;                      \
    } else {                                                   \
        PyErr_SetString(PyExc_TypeError,                       \
                        "Input object is not a complex_int8"); \
        return NULL;                                           \
    }                                                          \

#define PyComplexInt8_AsComplexInt8Pointer(c, o)               \
    if( PyComplexInt8_Check(o) ) {                             \
        c = &((PyComplexInt8 *) o)->obval;                     \
    } else {                                                   \
        PyErr_SetString(PyExc_TypeError,                       \
                        "Input object is not a complex_int8"); \
        return NULL;                                           \
    }    

static PyObject* pycomplexint8_new(PyTypeObject *type, PyObject *NPY_UNUSED(args), PyObject *NPY_UNUSED(kwds)) {
    PyComplexInt8 *self;
    self = (PyComplexInt8 *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

static int pycomplexint8_init(PyObject *self, PyObject *args, PyObject *kwds) {
    Py_ssize_t size = PyTuple_Size(args);
    complex_int8 *c;
    Py_complex cmplx;
    long real, imag;
    
    PyObject *C = {0};
    c = &(((PyComplexInt8 *) self)->obval);
    
    if( kwds != NULL && PyDict_Size(kwds) ) {
        PyErr_SetString(PyExc_TypeError,
                        "complex_int8 constructor takes no keyword arguments");
        return -1;
    }
    
    c->real_imag = 0;
    if( size == 0 ) {
        return 0;
    } else if( size == 1) {
        if( PyArg_ParseTuple(args, "O", &C) && PyComplexInt8_Check(C) ) {
            c->real_imag = ((PyComplexInt8 *) C)->obval.real_imag;
            return 0;
        } else if( PyArg_ParseTuple(args, "D", &cmplx) ) {
            inplace_pack_ci8(cmplx.real, cmplx.imag, c);
            return 0;
        } else if( PyArg_ParseTuple(args, "i", &real) ) {
            inplace_pack_ci8((signed char) real, 0, c);
            return 0;
        }
    } else if( size == 2 ) {
        if( PyArg_ParseTuple(args, "ii", &real, &imag) ) {
            inplace_pack_ci8((signed char) real, imag, c);
            return 0;
        }
    }
    
    PyErr_SetString(PyExc_TypeError,
                    "complex_int8 constructor takes zero or two int8 values, or a single argument of a complex_int8 or Python complex");
    return -1;
}

#define UNARY_BOOL_RETURNER_CI8(name)                                   \
    static PyObject*                                                    \
    pycomplexint8_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {        \
        complex_int8 c = {0};                                           \
        PyComplexInt8_AsComplexInt8(c, a);                              \
        return PyBool_FromLong(complex_int8_##name(c));                 \
    }
UNARY_BOOL_RETURNER_CI8(nonzero)
UNARY_BOOL_RETURNER_CI8(isnan)
UNARY_BOOL_RETURNER_CI8(isinf)
UNARY_BOOL_RETURNER_CI8(isfinite)

#define BINARY_BOOL_RETURNER_CI8(name)                                  \
    static PyObject*                                                    \
    pycomplexint8_##name(PyObject* a, PyObject* b) {                    \
        complex_int8 p = {0};                                           \
        complex_int8 q = {0};                                           \
        PyComplexInt8_AsComplexInt8(p, a);                              \
        PyComplexInt8_AsComplexInt8(q, b);                              \
        return PyBool_FromLong(complex_int8_##name(p,q));               \
    }
BINARY_BOOL_RETURNER_CI8(equal)
BINARY_BOOL_RETURNER_CI8(not_equal)
BINARY_BOOL_RETURNER_CI8(less)
BINARY_BOOL_RETURNER_CI8(greater)
BINARY_BOOL_RETURNER_CI8(less_equal)
BINARY_BOOL_RETURNER_CI8(greater_equal)

#define UNARY_INT_RETURNER_CI8(name)                                    \
    static PyObject*                                                    \
        pycomplexint8_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int8 q = {0};                                           \
        PyComplexInt8_AsComplexInt8(q, a);                              \
        return PyInt_FromLong(complex_int8_##name(q));                  \
    }
UNARY_INT_RETURNER_CI8(real)
UNARY_INT_RETURNER_CI8(imag)

#define UNARY_FLOAT_RETURNER_CI8(name)                                  \
    static PyObject*                                                    \
        pycomplexint8_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int8 q = {0};                                           \
        PyComplexInt8_AsComplexInt8(q, a);                              \
        return PyFloat_FromDouble(complex_int8_##name(q));              \
    }
UNARY_FLOAT_RETURNER_CI8(absolute)

#define UNARY_COMPLEX_INT8_RETURNER_CI8(name)                           \
    static PyObject*                                                    \
        pycomplexint8_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int8 q = {0};                                           \
        PyComplexInt8_AsComplexInt8(q, a);                              \
        return PyComplexInt8_FromComplexInt8(complex_int8_##name(q));   \
    }
UNARY_COMPLEX_INT8_RETURNER_CI8(negative)
UNARY_COMPLEX_INT8_RETURNER_CI8(conjugate)

static PyObject* pycomplexint8_positive(PyObject* self, PyObject* NPY_UNUSED(b)) {
    Py_INCREF(self);
    return self;
}

#define CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(fake_name, name)        \
    static PyObject*                                                          \
    pycomplexint8_##fake_name##_array_operator(PyObject* a, PyObject* b) {    \
        NpyIter *iter;                                                        \
        NpyIter_IterNextFunc *iternext;                                       \
        PyArrayObject *op[2];                                                 \
        PyObject *ret;                                                        \
        npy_uint32 flags;                                                     \
        npy_uint32 op_flags[2];                                               \
        PyArray_Descr *op_dtypes[2];                                          \
        npy_intp itemsize, *innersizeptr, innerstride;                        \
        char **dataptrarray;                                                  \
        char *src, *dst;                                                      \
        complex_int8 p = {0};                                                 \
        PyComplexInt8_AsComplexInt8(p, a);                                    \
        flags = NPY_ITER_EXTERNAL_LOOP;                                       \
        op[0] = (PyArrayObject *) b;                                          \
        op[1] = NULL;                                                         \
        op_flags[0] = NPY_ITER_READONLY;                                      \
        op_flags[1] = NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE;                 \
        op_dtypes[0] = PyArray_DESCR((PyArrayObject*) b);                     \
        op_dtypes[1] = complex_int8_descr;                                    \
        iter = NpyIter_MultiNew(2, op, flags, NPY_KEEPORDER, NPY_NO_CASTING, op_flags, op_dtypes); \
        if (iter == NULL) {                                                   \
            return NULL;                                                      \
        }                                                                     \
        iternext = NpyIter_GetIterNext(iter, NULL);                           \
        innerstride = NpyIter_GetInnerStrideArray(iter)[0];                   \
        itemsize = NpyIter_GetDescrArray(iter)[1]->elsize;                    \
        innersizeptr = NpyIter_GetInnerLoopSizePtr(iter);                     \
        dataptrarray = NpyIter_GetDataPtrArray(iter);                         \
        if(PyArray_EquivTypes(PyArray_DESCR((PyArrayObject*) b), complex_int8_descr)) { \
            npy_intp i;                                                       \
            do {                                                              \
                npy_intp size = *innersizeptr;                                \
                src = dataptrarray[0];                                        \
                dst = dataptrarray[1];                                        \
                for(i = 0; i < size; i++, src += innerstride, dst += itemsize) {  \
                    *((complex_int8 *) dst) = complex_int8_##name(p, *((complex_int8 *) src)); \
                }                                                             \
            } while (iternext(iter));                                         \
        } else if(PyArray_ISINTEGER((PyArrayObject*) b)) {                    \
            npy_intp i;                                                       \
            do {                                                              \
                npy_intp size = *innersizeptr;                                \
                src = dataptrarray[0];                                        \
                dst = dataptrarray[1];                                        \
                for(i = 0; i < size; i++, src += innerstride, dst += itemsize) {  \
                    *(complex_int8 *) dst = complex_int8_##name##_scalar(p, *((npy_long *) src)); \
                }                                                             \
            } while (iternext(iter));                                         \
        } else {                                                              \
            NpyIter_Deallocate(iter);                                         \
            return NULL;                                                      \
        }                                                                     \
        ret = (PyObject *) NpyIter_GetOperandArray(iter)[1];                  \
        Py_INCREF(ret);                                                       \
        if (NpyIter_Deallocate(iter) != NPY_SUCCEED) {                        \
            Py_DECREF(ret);                                                   \
            return NULL;                                                      \
        }                                                                     \
        return ret;                                                           \
    }                                                                         \
    static PyObject*                                                          \
    pycomplexint8_##fake_name(PyObject* a, PyObject* b) {                     \
        /* PyObject *a_type, *a_repr, *b_type, *b_repr, *a_repr2, *b_repr2;   \ */ \
        /* char* a_char, b_char, a_char2, b_char2;                            \ */ \
        npy_int64 val64;                                                      \
        npy_int32 val32;                                                      \
        complex_int8 p = {0};                                                 \
        if(PyArray_Check(b)) { return pycomplexint8_##fake_name##_array_operator(a, b); } \
        if(PyInt_Check(a) && PyComplexInt8_Check(b)) {                        \
            return PyComplexInt8_FromComplexInt8(complex_int8_scalar_##name(PyInt_AsLong(a), ((PyComplexInt8*)b)->obval));                                                 \
        }                                                                     \
        PyComplexInt8_AsComplexInt8(p, a);                                    \
        if(PyComplexInt8_Check(b)) {                                          \
            return PyComplexInt8_FromComplexInt8(complex_int8_##name(p,((PyComplexInt8*)b)->obval)); \
        } else if(PyInt_Check(b)) {                                           \
            return PyComplexInt8_FromComplexInt8(complex_int8_##name##_scalar(p,PyInt_AsLong(b))); \
        }                                                                     \
        PyErr_SetString(PyExc_TypeError, "Binary operation involving complex_int8 and neither integer nor complex_int8.");                                                              \
        return NULL;                                                          \
    }
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(add, add)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(subtract, subtract)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(multiply, multiply)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(divide, divide)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(true_divide, divide)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(floor_divide, divide)
/* CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_RETURNER(power, power) */

static PyObject* pycomplexint8_kludgy_arctan2_array_operator(PyObject* a, PyObject* b) {
    NpyIter *iter;
    NpyIter_IterNextFunc *iternext;
    PyArrayObject *op[2];
    PyObject *ret;
    npy_uint32 flags;
    npy_uint32 op_flags[2];
    PyArray_Descr *op_dtypes[2];
    npy_intp itemsize, *innersizeptr, innerstride;
    char **dataptrarray;
    char *src, *dst;
    long p = 0;
    p = PyInt_AsLong(a);
    flags = NPY_ITER_EXTERNAL_LOOP;
    op[0] = (PyArrayObject *) b;
    op[1] = NULL;
    op_flags[0] = NPY_ITER_READONLY;
    op_flags[1] = NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE;
    op_dtypes[0] = PyArray_DESCR((PyArrayObject*) b);
    op_dtypes[1] = complex_int8_descr;
    iter = NpyIter_MultiNew(2, op, flags, NPY_KEEPORDER, NPY_NO_CASTING, op_flags, op_dtypes);
    if (iter == NULL) {
        return NULL;
    }
    iternext = NpyIter_GetIterNext(iter, NULL);
    innerstride = NpyIter_GetInnerStrideArray(iter)[0];
    itemsize = NpyIter_GetDescrArray(iter)[1]->elsize;
    innersizeptr = NpyIter_GetInnerLoopSizePtr(iter);
    dataptrarray = NpyIter_GetDataPtrArray(iter);
    if(PyArray_EquivTypes(PyArray_DESCR((PyArrayObject*) b), complex_int8_descr)) {
        npy_intp i;
        do {
            npy_intp size = *innersizeptr;
            src = dataptrarray[0];
            dst = dataptrarray[1];
            for(i = 0; i < size; i++, src += innerstride, dst += itemsize) {
                *((double *) dst) = complex_int8_kludgy_arctan2(p, *((complex_int8 *) src));
            }
        } while (iternext(iter));
    } else {
        NpyIter_Deallocate(iter);
        return NULL;
    }
    ret = (PyObject *) NpyIter_GetOperandArray(iter)[1];
    Py_INCREF(ret);
    if (NpyIter_Deallocate(iter) != NPY_SUCCEED) {
        Py_DECREF(ret);
        return NULL;
    }
    return ret;
}

static PyObject* pycomplexint8_kludgy_arctan2(PyObject* a, PyObject* b) {
    /* PyObject *a_type, *a_repr, *b_type, *b_repr, *a_repr2, *b_repr2;   \ */
    /* char* a_char, b_char, a_char2, b_char2;                            \ */
    npy_int64 val64;
    npy_int32 val32;
    long p = 0;
    if(PyArray_Check(b)) { return pycomplexint8_kludgy_arctan2_array_operator(a, b); }
    p = PyInt_AsLong(a);
    if(PyComplexInt8_Check(b)) {
        return PyFloat_FromDouble(complex_int8_kludgy_arctan2(p,((PyComplexInt8*)b)->obval));
    }
    PyErr_SetString(PyExc_TypeError, "Binary operation involving complex_int8 and neither integer nor complex_int8.");
    return NULL;
}

#define CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_INPLACE(fake_name, name)        \
    static PyObject*                                                         \
    pycomplexint8_inplace_##fake_name(PyObject* a, PyObject* b) {            \
        complex_int8* p = {0};                                               \
        if(PyFloat_Check(a) || PyInt_Check(a)) {                             \
            PyErr_SetString(PyExc_TypeError, "Cannot in-place "#fake_name" a scalar by a complex_int8; should be handled by python.");                                                     \
            return NULL;                                                     \
        }                                                                    \
        PyComplexInt8_AsComplexInt8Pointer(p, a);                            \
        if(PyComplexInt8_Check(b)) {                                         \
            complex_int8_inplace_##name(p,((PyComplexInt8*)b)->obval);       \
            Py_INCREF(a);                                                    \
            return a;                                                        \
        } else if(PyInt_Check(b)) {                                          \
            complex_int8_inplace_##name##_scalar(p,PyInt_AsLong(b));         \
            Py_INCREF(a);                                                    \
            return a;                                                        \
        }                                                                    \
        PyErr_SetString(PyExc_TypeError, "Binary in-place operation involving complex_int8 and neither integer nor complex_int8.");                                                 \
        return NULL;                                                         \
    }
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_INPLACE(add, add)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_INPLACE(subtract, subtract)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_INPLACE(multiply, multiply)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_INPLACE(divide, divide)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_INPLACE(true_divide, divide)
CI8CI8_CI8S_SCI8_BINARY_COMPLEX_INT8_INPLACE(floor_divide, divide)

static PyObject* pycomplexint8__reduce(PyComplexInt8* self) {
    return Py_BuildValue("O(O)", Py_TYPE(self), PyInt_FromLong(self->obval.real_imag));
}

static PyObject* pycomplexint8_getstate(PyComplexInt8* self, PyObject* args) {
    if( !PyArg_ParseTuple(args, ":getstate") ) {
        return NULL;
    }
    return Py_BuildValue("O", PyInt_FromLong(self->obval.real_imag));
}

static PyObject* pycomplexint8_setstate(PyComplexInt8* self, PyObject* args) {
    complex_int8* c;
    c = &(self->obval);
    
    if( !PyArg_ParseTuple(args, "b:setstate", &c->real_imag) ) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// This is an array of methods (member functions) that will be
// available to use on the complex_int8 objects in python.  This is
// packaged up here, and will be used in the `tp_methods` field when
// definining the PyComplexInt8_Type below.
PyMethodDef pycomplexint8_methods[] = {
  // Unary bool returners
  {"nonzero", pycomplexint8_nonzero, METH_NOARGS,
   "True if the complex_int8 has all zero components"},
  {"isnan", pycomplexint8_isnan, METH_NOARGS,
   "True if the complex_int8 has any NAN components"},
  {"isinf", pycomplexint8_isinf, METH_NOARGS,
   "True if the complex_int8 has any INF components"},
  {"isfinite", pycomplexint8_isfinite, METH_NOARGS,
   "True if the complex_int8 has all finite components"},

  // Binary bool returners
  {"equal", pycomplexint8_equal, METH_O,
   "True if the complex_int8s are PRECISELY equal"},
  {"not_equal", pycomplexint8_not_equal, METH_O,
   "True if the complex_int8s are not PRECISELY equal"},
  {"less", pycomplexint8_less, METH_O,
   "Strict dictionary ordering"},
  {"greater", pycomplexint8_greater, METH_O,
   "Strict dictionary ordering"},
  {"less_equal", pycomplexint8_less_equal, METH_O,
   "Dictionary ordering"},
  {"greater_equal", pycomplexint8_greater_equal, METH_O,
   "Dictionary ordering"},
  
  // Unary ing returners
  {"real_part", pycomplexint8_real, METH_NOARGS,
   "real part of the complex value"},
  {"imag_part", pycomplexint8_imag, METH_NOARGS,
   "imaginary part of the complex value"},
  
  // Unary float returners
  {"absolute", pycomplexint8_absolute, METH_NOARGS,
   "Absolute value of complex_int8"},
  {"abs", pycomplexint8_absolute, METH_NOARGS,
   "Absolute value (Euclidean norm) of complex_int8"},
  
  // Unary complex_int8 returners
  // {"negative", pycomplexint8_negative, METH_NOARGS,
  //  "Return the negated complex_int8"},
  // {"positive", pycomplexint8_positive, METH_NOARGS,
  //  "Return the complex_int8 itself"},
  {"conjugate", pycomplexint8_conjugate, METH_NOARGS,
   "Return the complex conjugate of the complex_int8"},
  {"conj", pycomplexint8_conjugate, METH_NOARGS,
   "Return the complex conjugate of the complex_int8"},
  
  // complex_int8-complex_int8 or complex_int8-scalar binary complex_int8 returners
  // {"multiply", pycomplexint8_multiply, METH_O,
  //  "Standard (geometric) complex_int8 product"},
  // {"divide", pycomplexint8_divide, METH_O,
  //  "Standard (geometric) complex_int8 division"},
  // {"power", pycomplexint8_power, METH_O,
  //  "q.power(p) = (q.log() * p).exp()"},

  {"__reduce__", (PyCFunction)pycomplexint8__reduce, METH_NOARGS,
   "Return state information for pickling."},
  {"__getstate__", (PyCFunction)pycomplexint8_getstate, METH_VARARGS,
   "Return state information for pickling."},
  {"__setstate__", (PyCFunction)pycomplexint8_setstate, METH_VARARGS,
   "Reconstruct state information from pickle."},

  {NULL, NULL, 0, NULL}
};

static PyObject* pycomplexint8_num_negative(PyObject* a) { return pycomplexint8_negative(a,NULL); }
static PyObject* pycomplexint8_num_positive(PyObject* a) { return pycomplexint8_positive(a,NULL); }
static PyObject* pycomplexint8_num_absolute(PyObject* a) { return pycomplexint8_absolute(a,NULL); }
static int pycomplexint8_num_nonzero(PyObject* a) {
  complex_int8 q = ((PyComplexInt8*)a)->obval;
  return complex_int8_nonzero(q);
}
#define CANNOT_CONVERT_CI8(target)                                           \
  static PyObject* pycomplexint8_convert_##target(PyObject* a) {             \
    PyErr_SetString(PyExc_TypeError, "Cannot convert complex_int8 to " #target); \
    return NULL;                                                             \
  }
CANNOT_CONVERT_CI8(int)
CANNOT_CONVERT_CI8(float)
#if PY_MAJOR_VERSION < 3
CANNOT_CONVERT_CI8(long)
CANNOT_CONVERT_CI8(oct)
CANNOT_CONVERT_CI8(hex)
#endif

static PyNumberMethods pycomplexint8_as_number = {
  pycomplexint8_add,              // nb_add
  pycomplexint8_subtract,         // nb_subtract
  pycomplexint8_multiply,         // nb_multiply
  #if PY_MAJOR_VERSION < 3
  pycomplexint8_divide,           // nb_divide
  #endif
  0,                              // nb_remainder
  0,                              // nb_divmod
  0,                              // nb_power
  pycomplexint8_num_negative,     // nb_negative
  pycomplexint8_num_positive,     // nb_positive
  pycomplexint8_num_absolute,     // nb_absolute
  pycomplexint8_num_nonzero,      // nb_nonzero
  0,                              // nb_invert
  0,                              // nb_lshift
  0,                              // nb_rshift
  0,                              // nb_and
  0,                              // nb_xor
  0,                              // nb_or
  #if PY_MAJOR_VERSION < 3
  0,                              // nb_coerce
  #endif
  pycomplexint8_convert_int,      // nb_int
  #if PY_MAJOR_VERSION >= 3
  0,                              // nb_reserved
  #else
  pycomplexint8_convert_long,     // nb_long
  #endif
  pycomplexint8_convert_float,    // nb_float
  #if PY_MAJOR_VERSION < 3
  pycomplexint8_convert_oct,      // nb_oct
  pycomplexint8_convert_hex,      // nb_hex
  #endif
  pycomplexint8_inplace_add,      // nb_inplace_add
  pycomplexint8_inplace_subtract, // nb_inplace_subtract
  pycomplexint8_inplace_multiply, // nb_inplace_multiply
  #if PY_MAJOR_VERSION < 3
  pycomplexint8_inplace_divide,   // nb_inplace_divide
  #endif
  0,                              // nb_inplace_remainder
  0,                              // nb_inplace_power
  0,                              // nb_inplace_lshift
  0,                              // nb_inplace_rshift
  0,                              // nb_inplace_and
  0,                              // nb_inplace_xor
  0,                              // nb_inplace_or
  pycomplexint8_divide,           // nb_floor_divide
  pycomplexint8_divide,           // nb_true_divide
  pycomplexint8_inplace_divide,   // nb_inplace_floor_divide
  pycomplexint8_inplace_divide,   // nb_inplace_true_divide
  0,                              // nb_index
  #if PY_MAJOR_VERSION >= 3
  #if PY_MINOR_VERSION >= 5
  0,                              // nb_matrix_multiply
  0,                              // nb_inplace_matrix_multiply
  #endif
  #endif
};

// This is an array of members (member data) that will be available to
// use on the complex_int8 objects in python.  This is packaged up here,
// and will be used in the `tp_members` field when definining the
// PyComplexInt8_Type below.
PyMemberDef pycomplexint8_members[] = {
  {"real_imag", T_BYTE, offsetof(PyComplexInt8, obval.real_imag), 0,
   "The packed real and imaginary component of the complex_int8"},
  {NULL, 0, 0, 0, NULL}
};

static PyObject* pycomplexint8_richcompare(PyObject* a, PyObject* b, int op) {
    complex_int8 x = {0};
    complex_int8 y = {0};
    int result = 0;
    PyComplexInt8_AsComplexInt8(x,a);
    PyComplexInt8_AsComplexInt8(y,b);
    #define COMPARISONOP(py,op) case py: result = complex_int8_##op(x,y); break;
    switch (op) {
        COMPARISONOP(Py_LT,less)
        COMPARISONOP(Py_LE,less_equal)
        COMPARISONOP(Py_EQ,equal)
        COMPARISONOP(Py_NE,not_equal)
        COMPARISONOP(Py_GT,greater)
        COMPARISONOP(Py_GE,greater_equal)
    };
    #undef COMPARISONOP
    return PyBool_FromLong(result);
}

static long pycomplexint8_hash(PyObject *o) {
    complex_int8 q = ((PyComplexInt8 *)o)->obval;
    long value = 0x456789;
    value = (10000004 * value) ^ _Py_HashDouble(q.real_imag);
    if (value == -1) {
        value = -2;
    }
    return value;
}

static PyObject* pycomplexint8_repr(PyObject *o) {
    char str[128];
    complex_int8 c = ((PyComplexInt8 *)o)->obval;
    const signed char* sc = fourBitLUT[c.real_imag];
    sprintf(str, "complex_int8(%u [%i, %i])", c.real_imag, sc[0], sc[1]);
    return PyUString_FromString(str);
}

static PyObject* pycomplexint8_str(PyObject *o) {
    char str[64];
    complex_int8 c = ((PyComplexInt8 *)o)->obval;
    const signed char* sc = fourBitLUT[c.real_imag];
    sprintf(str, "%i%+ij", sc[0], sc[1]);
    return PyUString_FromString(str);
}

// This establishes the complex_int8 as a python object (not yet a numpy
// scalar type).  The name may be a little counterintuitive; the idea
// is that this will be a type that can be used as an array dtype.
// Note that many of the slots below will be filled later, after the
// corresponding functions are defined.
static PyTypeObject PyComplexInt8_Type = {
#if PY_MAJOR_VERSION >= 3
  PyVarObject_HEAD_INIT(NULL, 0)
#else
  PyObject_HEAD_INIT(NULL)
  0,                                          // ob_size
#endif
  "complex_int8.complex_int8",                // tp_name
  sizeof(PyComplexInt8),                      // tp_basicsize
  0,                                          // tp_itemsize
  0,                                          // tp_dealloc
  0,                                          // tp_print
  0,                                          // tp_getattr
  0,                                          // tp_setattr
#if PY_MAJOR_VERSION >= 3
  0,                                          // tp_reserved
#else
  0,                                          // tp_compare
#endif
  pycomplexint8_repr,                         // tp_repr
  &pycomplexint8_as_number,                   // tp_as_number
  0,                                          // tp_as_sequence
  0,                                          // tp_as_mapping
  pycomplexint8_hash,                         // tp_hash
  0,                                          // tp_call
  pycomplexint8_str,                          // tp_str
  0,                                          // tp_getattro
  0,                                          // tp_setattro
  0,                                          // tp_as_buffer
#if PY_MAJOR_VERSION >= 3
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   // tp_flags
#else
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES, // tp_flags
#endif
  "Complex integer complex_int8 numbers",     // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  pycomplexint8_richcompare,                  // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  pycomplexint8_methods,                      // tp_methods
  pycomplexint8_members,                      // tp_members
  0,                                          // tp_getset
  0,                                          // tp_base; will be reset to &PyGenericArrType_Type after numpy import
  0,                                          // tp_dict
  0,                                          // tp_descr_get
  0,                                          // tp_descr_set
  0,                                          // tp_dictoffset
  pycomplexint8_init,                         // tp_init
  0,                                          // tp_alloc
  pycomplexint8_new,                          // tp_new
  0,                                          // tp_free
  0,                                          // tp_is_gc
  0,                                          // tp_bases
  0,                                          // tp_mro
  0,                                          // tp_cache
  0,                                          // tp_subclasses
  0,                                          // tp_weaklist
  0,                                          // tp_del
#if PY_VERSION_HEX >= 0x02060000
  0,                                          // tp_version_tag
#endif
#if PY_VERSION_HEX >= 0x030400a1
  0,                                          // tp_finalize
#endif
};

// Functions implementing internal features. Not all of these function
// pointers must be defined for a given type. The required members are
// nonzero, copyswap, copyswapn, setitem, getitem, and cast.
static PyArray_ArrFuncs _PyComplexInt8_ArrFuncs;
PyArray_Descr *complex_int8_descr;

static npy_bool CI8_nonzero (char *ip, PyArrayObject *ap) {
  complex_int8 c;
  complex_int8 zero = {0};
  if (ap == NULL || PyArray_ISBEHAVED_RO(ap)) {
    c = *(complex_int8 *)ip;
  }
  else {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT8);
    descr->f->copyswap(&c.real_imag, ip, !PyArray_ISNOTSWAPPED(ap), NULL);
    Py_DECREF(descr);
  }
  return (npy_bool) !complex_int8_equal(c, zero);
}

static void CI8_copyswap(complex_int8 *dst, complex_int8 *src, int swap, void *NPY_UNUSED(arr)) {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT8);
    descr->f->copyswapn(dst, sizeof(unsigned char), src, sizeof(unsigned char), 1, swap, NULL);
    Py_DECREF(descr);
}

static void CI8_copyswapn(complex_int8 *dst, npy_intp dstride,
                               complex_int8 *src, npy_intp sstride,
                               npy_intp n, int swap, void *NPY_UNUSED(arr)) {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT8);
    descr->f->copyswapn(&dst->real_imag, dstride, &src->real_imag, sstride, n, swap, NULL);
    Py_DECREF(descr);    
}


static int CI8_setitem(PyObject* item, complex_int8* c, void* NPY_UNUSED(ap)) {
    PyObject *element;
    long real, imag;
    if( PyComplexInt8_Check(item) ) {
        memcpy(c, &(((PyComplexInt8 *)item)->obval), sizeof(complex_int8));
    } else if( PyComplex_Check(item) ) {
        inplace_pack_ci8(PyComplex_RealAsDouble(item), \
                         PyComplex_ImagAsDouble(item), \
                         c);
    } else if( PyInt_Check(item) ) {
        inplace_pack_ci8(PyInt_AsLong(item), \
                         0, \
                         c);
    } else if( PySequence_Check(item) && PySequence_Length(item) == 2 ) {
        element = PySequence_GetItem(item, 0);
        if(element == NULL) {
            return -1;
        } /* Not a sequence, or other failure */
        real = PyInt_AsLong(element);
        Py_DECREF(element);
        
        element = PySequence_GetItem(item, 1);
        if(element == NULL) {
            return -1;
        } /* Not a sequence, or other failure */
        imag = PyInt_AsLong(element);
        Py_DECREF(element);
        
        inplace_pack_ci8(real, imag, c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Unknown input to CI8_setitem");
        return -1;
    }
    return 0;
}

// When a numpy array of dtype=complex_int8 is indexed, this function is
// called, returning a new complex_int8 object with a copy of the
// data... sometimes...
static PyObject* CI8_getitem(void* data, void* NPY_UNUSED(arr)) {
    complex_int8 q;
    memcpy(&q, data, sizeof(complex_int8));
    return PyComplexInt8_FromComplexInt8(q);
}

static int CI8_compare(complex_int8 *pa, complex_int8 *pb, PyArrayObject *NPY_UNUSED(ap)) {
    complex_int8 a = *pa, b = *pb;
    npy_bool anan, bnan;
    int ret;
    
    anan = complex_int8_isnan(a);
    bnan = complex_int8_isnan(b);
    
    if( anan ) {
        ret = bnan ? 0 : -1;
    } else if( bnan ) {
        ret = 1;
    } else if( complex_int8_less(a, b) ) {
        ret = -1;
    } else if( complex_int8_less(b, a) ) {
        ret = 1;
    } else {
        ret = 0;
    }
    
    return ret;
}

static int CI8_argmax(complex_int8 *ip, npy_intp n, npy_intp *max_ind, PyArrayObject *NPY_UNUSED(aip)) {
    npy_intp i;
    complex_int8 mp = *ip;
    
    *max_ind = 0;
    
    if( complex_int8_isnan(mp) ) {
        /* nan encountered; it's maximal */
        return 0;
    }

    for(i = 1; i < n; i++) {
        ip++;
        /*
         * Propagate nans, similarly as max() and min()
         */
        if( !(complex_int8_less_equal(*ip, mp)) ) {  /* negated, for correct nan handling */
            mp = *ip;
            *max_ind = i;
            if (complex_int8_isnan(mp)) {
                /* nan encountered, it's maximal */
                break;
            }
        }
    }
    return 0;
}

static void CI8_fillwithscalar(complex_int8 *buffer, npy_intp length, complex_int8 *value, void *NPY_UNUSED(ignored)) {
    npy_intp i;
    complex_int8 val = *value;

    for (i = 0; i < length; ++i) {
        buffer[i] = val;
    }
}

// This is a macro (followed by applications of the macro) that cast
// the input types to standard complex_int8 with only a nonzero scalar
// part.
#define MAKE_T_TO_CI8(TYPE, type)                                      \
static void TYPE ## _to_complex_int8(type *ip, complex_int8 *op, npy_intp n, \
                                  PyArrayObject *NPY_UNUSED(aip),      \
                                  PyArrayObject *NPY_UNUSED(aop)) {    \
    while (n--) {                                                      \
        op->real_imag = (unsigned char) ((*ip++) * 16);                \
        *op++;                                                         \
    }                                                                  \
}

MAKE_T_TO_CI8(BOOL, npy_bool);
MAKE_T_TO_CI8(BYTE, npy_byte);

// This is a macro (followed by applications of the macro) that cast
// the input complex types from complex_int8.
#define MAKE_CI8_TO_CT(TYPE, type)                                    \
static void complex_int8_to_## TYPE(complex_int8* ip, type *op, npy_intp n, \
                                PyArrayObject *NPY_UNUSED(aip),       \
                                PyArrayObject *NPY_UNUSED(aop)) {     \
    const signed char* sc;                                            \
    while (n--) {                                                     \
        sc = fourBitLUT[ip->real_imag];                               \
        op->real = sc[0];                                             \
        op->imag = sc[1];                                             \
        (*ip++); (*op++);                                             \
    }                                                                 \
}

MAKE_CI8_TO_CT(CFLOAT, npy_cfloat);
MAKE_CI8_TO_CT(CDOUBLE, npy_cdouble);
MAKE_CI8_TO_CT(CLONGDOUBLE, npy_clongdouble);

static void register_cast_function_ci8(int sourceType, int destType, PyArray_VectorUnaryFunc *castfunc) {
    PyArray_Descr *descr = PyArray_DescrFromType(sourceType);
    PyArray_RegisterCastFunc(descr, destType, castfunc);
    PyArray_RegisterCanCast(descr, destType, NPY_NOSCALAR);
    Py_DECREF(descr);
}

// This is a macro that will be used to define the various basic unary
// complex_int8 functions, so that they can be applied quickly to a
// numpy array of complex_int8s.
#define UNARY_GEN_UFUNC_CI8(ufunc_name, func_name, ret_type)                \
    static void complex_int8_##ufunc_name##_ufunc(char** args,              \
                                                  npy_intp* dimensions,     \
                                                  npy_intp* steps,          \
                                                  void* NPY_UNUSED(data)) { \
    char *ip1 = args[0], *op1 = args[1];                                    \
    npy_intp is1 = steps[0], os1 = steps[1];                                \
    npy_intp n = dimensions[0];                                             \
    npy_intp i;                                                             \
    for(i = 0; i < n; i++, ip1 += is1, op1 += os1){                         \
      const complex_int8 in1 = *(complex_int8 *)ip1;                        \
      *((ret_type *)op1) = complex_int8_##func_name(in1);};}
#define UNARY_UFUNC_CI8(name, ret_type) \
  UNARY_GEN_UFUNC_CI8(name, name, ret_type)
// And these all do the work mentioned above, using the macro
UNARY_UFUNC_CI8(isnan, npy_bool)
UNARY_UFUNC_CI8(isinf, npy_bool)
UNARY_UFUNC_CI8(isfinite, npy_bool)
UNARY_UFUNC_CI8(real, npy_long)
UNARY_UFUNC_CI8(imag, npy_long)
UNARY_UFUNC_CI8(absolute, npy_double)
UNARY_UFUNC_CI8(negative, complex_int8)
UNARY_UFUNC_CI8(conjugate, complex_int8)
static void complex_int8_positive_ufunc(char** args, npy_intp* dimensions, npy_intp* steps, void* NPY_UNUSED(data)) {
    char *ip1 = args[0], *op1 = args[1];
    npy_intp is1 = steps[0], os1 = steps[1];
    npy_intp n = dimensions[0];
    npy_intp i;
    for(i = 0; i < n; i++, ip1 += is1, op1 += os1) {
        const complex_int8 in1 = *(complex_int8 *)ip1;
        *((complex_int8 *)op1) = in1;
    }
}

// This is a macro that will be used to define the various basic binary
// complex_int8 functions, so that they can be applied quickly to a
// numpy array of complex_int8s.
#define BINARY_GEN_UFUNC_CI8(ufunc_name, func_name, arg_type1, arg_type2, ret_type) \
  static void complex_int8_##ufunc_name##_ufunc(char** args,                    \
                                               npy_intp* dimensions,            \
                                               npy_intp* steps,                 \
                                               void* NPY_UNUSED(data)) {        \
    char *ip1 = args[0], *ip2 = args[1], *op1 = args[2];                        \
    npy_intp is1 = steps[0], is2 = steps[1], os1 = steps[2];                    \
    npy_intp n = dimensions[0];                                                 \
    npy_intp i;                                                                 \
    for(i = 0; i < n; i++, ip1 += is1, ip2 += is2, op1 += os1) {                \
      const arg_type1 in1 = *(arg_type1 *)ip1;                                  \
      const arg_type2 in2 = *(arg_type2 *)ip2;                                  \
      *((ret_type *)op1) = complex_int8_##func_name(in1, in2);                  \
    };                                                                          \
  };
// A couple special-case versions of the above
#define BINARY_UFUNC_CI8(name, ret_type)                    \
  BINARY_GEN_UFUNC_CI8(name, name, complex_int8, complex_int8, ret_type)
#define BINARY_SCALAR_UFUNC_CI8(name, ret_type)                             \
  BINARY_GEN_UFUNC_CI8(name##_scalar, name##_scalar, complex_int8, npy_long, ret_type) \
  BINARY_GEN_UFUNC_CI8(scalar_##name, scalar_##name, npy_long, complex_int8, ret_type)
// And these all do the work mentioned above, using the macros
BINARY_UFUNC_CI8(add, complex_int8)
BINARY_UFUNC_CI8(subtract, complex_int8)
BINARY_UFUNC_CI8(multiply, complex_int8)
BINARY_UFUNC_CI8(divide, complex_int8)
BINARY_GEN_UFUNC_CI8(true_divide, divide, complex_int8, complex_int8, complex_int8)
BINARY_GEN_UFUNC_CI8(floor_divide, divide, complex_int8, complex_int8, complex_int8)
BINARY_UFUNC_CI8(equal, npy_bool)
BINARY_UFUNC_CI8(not_equal, npy_bool)
BINARY_UFUNC_CI8(less, npy_bool)
BINARY_UFUNC_CI8(less_equal, npy_bool)
BINARY_SCALAR_UFUNC_CI8(add, complex_int8)
BINARY_SCALAR_UFUNC_CI8(subtract, complex_int8)
BINARY_SCALAR_UFUNC_CI8(multiply, complex_int8)
BINARY_SCALAR_UFUNC_CI8(divide, complex_int8)
BINARY_GEN_UFUNC_CI8(true_divide_scalar, divide_scalar, complex_int8, long, complex_int8)
BINARY_GEN_UFUNC_CI8(floor_divide_scalar, divide_scalar, complex_int8, long, complex_int8)
BINARY_GEN_UFUNC_CI8(scalar_true_divide, scalar_divide, long, complex_int8, complex_int8)
BINARY_GEN_UFUNC_CI8(scalar_floor_divide, scalar_divide, long, complex_int8, complex_int8)

BINARY_GEN_UFUNC_CI8(arctan2, kludgy_arctan2, long, complex_int8, npy_double)

static PyObject* complex_int8_arrtype_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    complex_int8 c;
    Py_complex cmplx;

    if( !PyArg_ParseTuple(args, "D", &cmplx) ) {
        return NULL;
    }
    
    c = pack_ci8(cmplx.real, cmplx.imag);
    return PyArray_Scalar(&c, complex_int8_descr, NULL);
}

static PyObject* gentype_richcompare_ci8(PyObject *self, PyObject *other, int cmp_op) {
    PyObject *arr, *ret;

    arr = PyArray_FromScalar(self, NULL);
    if (arr == NULL) {
        return NULL;
    }
    ret = Py_TYPE(arr)->tp_richcompare(arr, other, cmp_op);
    Py_DECREF(arr);
    return ret;
}

static long complex_int8_arrtype_hash(PyObject *o) {
    complex_int8 c = ((PyComplexInt8 *)o)->obval;
    long value = 0x456789;
    value = (10000004 * value) ^ _Py_HashDouble(c.real_imag);
    if( value == -1 ) {
        value = -2;
    }
    return value;
}

static PyObject* complex_int8_arrtype_repr(PyObject *o) {
    char str[64];
    complex_int8 c = ((PyComplexInt8 *)o)->obval;
    const signed char* sc = fourBitLUT[c.real_imag];
    sprintf(str, "complex_int8(%u [%i, %i])", c.real_imag, sc[0], sc[1]);
    return PyUString_FromString(str);
}

static PyObject* complex_int8_arrtype_str(PyObject *o) {
    char str[64];
    complex_int8 c = ((PyComplexInt8 *)o)->obval;
    const signed char* sc = fourBitLUT[c.real_imag];
    sprintf(str, "%i%+ij", sc[0], sc[1]);
    return PyUString_FromString(str);
}

int complex_int8_elsize = sizeof(complex_int8);

typedef struct {
    char c;
    complex_int8 q;
} align_test_ci8;
int complex_int8_alignment = offsetof(align_test_ci8, q);

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//                                                             //
//  Everything above was preparation for the following set up  //
//                                                             //
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

int create_complex_int8(PyObject* m, PyObject* numpy_dict) {
    int complexi8Num;
    PyObject *tmp_ufunc;
    int arg_types[3];
    
    /* Fill the lookup table */
    complex_int8_fillLUT();
    
    /* Register the complex_int8 array base type. */
    PyComplexInt8_Type.tp_base = &PyGenericArrType_Type;
    if (PyType_Ready(&PyComplexInt8_Type) < 0) {
        return -2;
    }
    
    // The array functions, to be used below.  This InitArrFuncs
    // function is a convenient way to set all the fields to zero
    // initially, so we don't get undefined behavior.
    PyArray_InitArrFuncs(&_PyComplexInt8_ArrFuncs);
    _PyComplexInt8_ArrFuncs.nonzero = (PyArray_NonzeroFunc*)CI8_nonzero;
    _PyComplexInt8_ArrFuncs.copyswap = (PyArray_CopySwapFunc*)CI8_copyswap;
    _PyComplexInt8_ArrFuncs.copyswapn = (PyArray_CopySwapNFunc*)CI8_copyswapn;
    _PyComplexInt8_ArrFuncs.setitem = (PyArray_SetItemFunc*)CI8_setitem;
    _PyComplexInt8_ArrFuncs.getitem = (PyArray_GetItemFunc*)CI8_getitem;
    _PyComplexInt8_ArrFuncs.compare = (PyArray_CompareFunc*)CI8_compare;
    _PyComplexInt8_ArrFuncs.argmax = (PyArray_ArgFunc*)CI8_argmax;
    _PyComplexInt8_ArrFuncs.fillwithscalar = (PyArray_FillWithScalarFunc*)CI8_fillwithscalar;
    
    // The complex_int8 array descr
    complex_int8_descr = PyObject_New(PyArray_Descr, &PyArrayDescr_Type);
    complex_int8_descr->typeobj = &PyComplexInt8_Type;
    complex_int8_descr->kind = 'V';
    complex_int8_descr->type = 'i';
    complex_int8_descr->byteorder = '=';
    complex_int8_descr->flags = NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM;
    complex_int8_descr->type_num = 0; // assigned at registration
    complex_int8_descr->elsize = complex_int8_elsize;
    complex_int8_descr->alignment = complex_int8_alignment;
    complex_int8_descr->subarray = NULL;
    complex_int8_descr->fields = NULL;
    complex_int8_descr->names = NULL;
    complex_int8_descr->f = &_PyComplexInt8_ArrFuncs;
    complex_int8_descr->metadata = NULL;
    complex_int8_descr->c_metadata = NULL;
    
    Py_INCREF(&PyComplexInt8_Type);
    complexi8Num = PyArray_RegisterDataType(complex_int8_descr);
    
    if( complexi8Num < 0 || complexi8Num != NPY_COMPLEX_INT8 ) {
        return -1;
    }
    
    register_cast_function_ci8(NPY_BOOL, complexi8Num, (PyArray_VectorUnaryFunc*)BOOL_to_complex_int8);
    register_cast_function_ci8(NPY_BYTE, complexi8Num, (PyArray_VectorUnaryFunc*)BYTE_to_complex_int8);
    
    register_cast_function_ci8(complexi8Num, NPY_CFLOAT, (PyArray_VectorUnaryFunc*)complex_int8_to_CFLOAT);
    register_cast_function_ci8(complexi8Num, NPY_CDOUBLE, (PyArray_VectorUnaryFunc*)complex_int8_to_CDOUBLE);
    register_cast_function_ci8(complexi8Num, NPY_CLONGDOUBLE, (PyArray_VectorUnaryFunc*)complex_int8_to_CLONGDOUBLE);
    
    // These macros will be used below
#define REGISTER_UFUNC_CI8(name)\
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name),\
            complex_int8_descr->type_num, complex_int8_##name##_ufunc, arg_types, NULL)
    
#define REGISTER_SCALAR_UFUNC_CI8(name)\
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name),\
            complex_int8_descr->type_num, complex_int8_scalar_##name##_ufunc, arg_types, NULL)
    
#define REGISTER_UFUNC_SCALAR_CI8(name)                                   \
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name), \
                                complex_int8_descr->type_num, complex_int8_##name##_scalar_ufunc, arg_types, NULL)
    
#define REGISTER_NEW_UFUNC_GENERAL_CI8(pyname, cname, nargin, nargout, doc) \
    tmp_ufunc = PyUFunc_FromFuncAndData(NULL, NULL, NULL, 0, nargin, nargout, \
                                        PyUFunc_None, #pyname, doc, 0); \
    PyUFunc_RegisterLoopForType((PyUFuncObject *)tmp_ufunc,             \
                                complex_int8_descr->type_num, complex_int8_##cname##_ufunc, arg_types, NULL); \
    PyDict_SetItemString(numpy_dict, #pyname, tmp_ufunc);               \
    Py_DECREF(tmp_ufunc)

  #define REGISTER_NEW_UFUNC_CI8(name, nargin, nargout, doc)                \
    REGISTER_NEW_UFUNC_GENERAL_CI8(name, name, nargin, nargout, doc)
    
    /* complex_int8 -> bool */
    arg_types[0] = complex_int8_descr->type_num;
    arg_types[1] = NPY_BOOL;
    
    REGISTER_UFUNC_CI8(isnan);
    REGISTER_UFUNC_CI8(isinf);
    REGISTER_UFUNC_CI8(isfinite);
    
    /* complex_int8 -> long */
    arg_types[1] = NPY_LONG;
    
    REGISTER_NEW_UFUNC_GENERAL_CI8(real_part, real, 1, 1, \
                                   "Return the real part of the complex number\n");
    REGISTER_NEW_UFUNC_GENERAL_CI8(imag_part, imag, 1, 1, \
                                   "Return the imaginary part of the complex number\n");
    
    /* complex_int8 -> double */
    arg_types[1] = NPY_DOUBLE;
    
    REGISTER_UFUNC_CI8(absolute);
    
    /* complex_int8 -> complex_int8 */
    arg_types[1] = complex_int8_descr->type_num;
    
    REGISTER_UFUNC_CI8(negative);
    REGISTER_UFUNC_CI8(conjugate);

    /* complex_int8, complex_int8 -> bool */
    arg_types[2] = NPY_BOOL;

    REGISTER_UFUNC_CI8(equal);
    REGISTER_UFUNC_CI8(not_equal);
    REGISTER_UFUNC_CI8(less);
    REGISTER_UFUNC_CI8(less_equal);
    
    /* complex_int8, complex_int8 -> complex_int8 */
    arg_types[2] = complex_int8_descr->type_num;
    
    REGISTER_UFUNC_CI8(add);
    REGISTER_UFUNC_CI8(subtract);
    REGISTER_UFUNC_CI8(multiply);
    REGISTER_UFUNC_CI8(divide);
    REGISTER_UFUNC_CI8(true_divide);
    REGISTER_UFUNC_CI8(floor_divide);
    
    /* long, complex_int8 -> double */
    arg_types[0] = NPY_LONG;
    arg_types[1] = complex_int8_descr->type_num;
    arg_types[2] = NPY_DOUBLE;
    
    REGISTER_UFUNC_CI8(arctan2);
    
    /* long, complex_int8 -> complex_int8 */
    arg_types[0] = NPY_LONG;
    arg_types[1] = complex_int8_descr->type_num;
    arg_types[2] = complex_int8_descr->type_num;
    
    REGISTER_SCALAR_UFUNC_CI8(add);
    REGISTER_SCALAR_UFUNC_CI8(subtract);
    REGISTER_SCALAR_UFUNC_CI8(multiply);
    REGISTER_SCALAR_UFUNC_CI8(divide);
    REGISTER_SCALAR_UFUNC_CI8(true_divide);
    REGISTER_SCALAR_UFUNC_CI8(floor_divide);
    
    /* complex_int8, long -> complex_int8 */
    arg_types[0] = complex_int8_descr->type_num;
    arg_types[1] = NPY_LONG;
    arg_types[2] = complex_int8_descr->type_num;
    
    REGISTER_UFUNC_SCALAR_CI8(add);
    REGISTER_UFUNC_SCALAR_CI8(subtract);
    REGISTER_UFUNC_SCALAR_CI8(multiply);
    REGISTER_UFUNC_SCALAR_CI8(divide);
    REGISTER_UFUNC_SCALAR_CI8(true_divide);
    REGISTER_UFUNC_SCALAR_CI8(floor_divide);
    
    PyModule_AddObject(m, "complex_int8", (PyObject *)&PyComplexInt8_Type);
    
    return complexi8Num;
}
