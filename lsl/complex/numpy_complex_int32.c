#include "complex_int32.h"

typedef struct {
    PyObject_HEAD
    complex_int32 obval;
} PyComplexInt32;

static PyTypeObject PyComplexInt32_Type;

PyArray_Descr *complex_int32_descr;

static NPY_INLINE int PyComplexInt32_Check(PyObject *o) {
    return PyObject_IsInstance(o, (PyObject*) &PyComplexInt32_Type);
}

static PyObject* PyComplexInt32_FromComplexInt32(complex_int32 c) {
    PyComplexInt32 *p = (PyComplexInt32 *) PyComplexInt32_Type.tp_alloc(&PyComplexInt32_Type, 0);
    if( p != NULL ) {
        p->obval = c;
    }
    return (PyObject *) p;
}

#define PyComplexInt32_AsComplexInt32(c, o)                     \
    if( PyComplexInt32_Check(o) ) {                             \
        c = ((PyComplexInt32 *) o)->obval;                      \
    } else {                                                    \
        PyErr_SetString(PyExc_TypeError,                        \
                        "Input object is not a complex_int32"); \
        return NULL;                                            \
    }                                                           \

#define PyComplexInt32_AsComplexInt32Pointer(c, o)              \
    if( PyComplexInt32_Check(o) ) {                             \
        c = &((PyComplexInt32 *) o)->obval;                     \
    } else {                                                    \
        PyErr_SetString(PyExc_TypeError,                        \
                        "Input object is not a complex_int32"); \
        return NULL;                                            \
    }    

static PyObject* pycomplexint32_new(PyTypeObject *type, PyObject *NPY_UNUSED(args), PyObject *NPY_UNUSED(kwds)) {
    PyComplexInt32 *self;
    self = (PyComplexInt32 *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

static int pycomplexint32_init(PyObject *self, PyObject *args, PyObject *kwds) {
    Py_ssize_t size = PyTuple_Size(args);
    complex_int32 *c;
    Py_complex cmplx;
    long real, imag;
    
    PyObject *C = {0};
    c = &(((PyComplexInt32 *) self)->obval);
    
    if( kwds != NULL && PyDict_Size(kwds) ) {
        PyErr_SetString(PyExc_TypeError,
                        "complex_int32 constructor takes no keyword arguments");
        return -1;
    }
    
    c->real = 0;
    c->imag = 0;
    if( size == 0 ) {
        return 0;
    } else if( size == 1) {
        if( PyArg_ParseTuple(args, "O", &C) && (PyComplexInt8_Check(C) || PyComplexInt16_Check(C) || PyComplexInt32_Check(C)) ) {
            if( PyComplexInt8_Check(C) ) {
                const signed char* sc = fourBitLUT[((PyComplexInt8 *) C)->obval.real_imag];
                c->real = sc[0];
                c->imag = sc[1];
                return 0;
            } else if( PyComplexInt16_Check(C) ) {
                c->real = ((PyComplexInt16 *) C)->obval.real;
                c->imag = ((PyComplexInt16 *) C)->obval.imag;
                return 0;
            } else if( PyComplexInt32_Check(C) ) {
                c->real = ((PyComplexInt32 *) C)->obval.real;
                c->imag = ((PyComplexInt32 *) C)->obval.imag;
                return 0;
            }
        } else if( PyArg_ParseTuple(args, "D", &cmplx) ) {
            c->real = cmplx.real;
            c->imag = cmplx.imag;
            return 0;
        } else if( PyArg_ParseTuple(args, "i", &real) ) {
            c->real = real;
            c->imag = 0;
            return 0;
        }
    } else if( size == 2 ) {
        if( PyArg_ParseTuple(args, "ii", &real, &imag) ) {
            c->real = real;
            c->imag = imag;
            return 0;
        }
    }
    
    PyErr_SetString(PyExc_TypeError,
                    "complex_int32 constructor takes zero or two int32 values, or a single argument of a complex_int32 or Python complex");
    return -1;
}

#define UNARY_BOOL_RETURNER_CI32(name)                                       \
    static PyObject*                                                    \
    pycomplexint32_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {       \
        complex_int32 c = {0, 0};                                          \
        PyComplexInt32_AsComplexInt32(c, a);                            \
        return PyBool_FromLong(complex_int32_##name(c));                \
    }
UNARY_BOOL_RETURNER_CI32(nonzero)
UNARY_BOOL_RETURNER_CI32(isnan)
UNARY_BOOL_RETURNER_CI32(isinf)
UNARY_BOOL_RETURNER_CI32(isfinite)

#define BINARY_BOOL_RETURNER_CI32(name)                                      \
    static PyObject*                                                    \
    pycomplexint32_##name(PyObject* a, PyObject* b) {                   \
        complex_int32 p = {0, 0};                                          \
        complex_int32 q = {0, 0};                                          \
        PyComplexInt32_AsComplexInt32(p, a);                            \
        PyComplexInt32_AsComplexInt32(q, b);                            \
        return PyBool_FromLong(complex_int32_##name(p,q));              \
    }
BINARY_BOOL_RETURNER_CI32(equal)
BINARY_BOOL_RETURNER_CI32(not_equal)
BINARY_BOOL_RETURNER_CI32(less)
BINARY_BOOL_RETURNER_CI32(greater)
BINARY_BOOL_RETURNER_CI32(less_equal)
BINARY_BOOL_RETURNER_CI32(greater_equal)

#define UNARY_INT_RETURNER_CI32(name)                                      \
    static PyObject*                                                    \
        pycomplexint32_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int32 q = {0, 0};                                           \
        PyComplexInt32_AsComplexInt32(q, a);                              \
        return PyInt_FromLong(complex_int32_##name(q));              \
    }
UNARY_INT_RETURNER_CI32(real)
UNARY_INT_RETURNER_CI32(imag)

#define UNARY_FLOAT_RETURNER_CI32(name)                                      \
    static PyObject*                                                    \
        pycomplexint32_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int32 q = {0, 0};                                           \
        PyComplexInt32_AsComplexInt32(q, a);                              \
        return PyFloat_FromDouble(complex_int32_##name(q));              \
    }
UNARY_FLOAT_RETURNER_CI32(absolute)
UNARY_FLOAT_RETURNER_CI32(norm)
UNARY_FLOAT_RETURNER_CI32(angle)

#define UNARY_COMPLEX_INT32_RETURNER_CI32(name)                               \
    static PyObject*                                                    \
        pycomplexint32_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int32 q = {0, 0};                          \
        PyComplexInt32_AsComplexInt32(q, a);                              \
        return PyComplexInt32_FromComplexInt32(complex_int32_##name(q));   \
    }
UNARY_COMPLEX_INT32_RETURNER_CI32(negative)
UNARY_COMPLEX_INT32_RETURNER_CI32(conjugate)

static PyObject* pycomplexint32_positive(PyObject* self, PyObject* NPY_UNUSED(b)) {
    Py_INCREF(self);
    return self;
}

#define CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(fake_name, name)        \
    static PyObject*                                                          \
    pycomplexint32_##fake_name##_array_operator(PyObject* a, PyObject* b) {    \
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
        complex_int32 p = {0, 0};                                             \
        PyComplexInt32_AsComplexInt32(p, a);                                    \
        flags = NPY_ITER_EXTERNAL_LOOP;                                       \
        op[0] = (PyArrayObject *) b;                                          \
        op[1] = NULL;                                                         \
        op_flags[0] = NPY_ITER_READONLY;                                      \
        op_flags[1] = NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE;                 \
        op_dtypes[0] = PyArray_DESCR((PyArrayObject*) b);                     \
        op_dtypes[1] = complex_int32_descr;                                    \
        iter = NpyIter_MultiNew(2, op, flags, NPY_KEEPORDER, NPY_NO_CASTING, op_flags, op_dtypes); \
        if (iter == NULL) {                                                   \
            return NULL;                                                      \
        }                                                                     \
        iternext = NpyIter_GetIterNext(iter, NULL);                           \
        innerstride = NpyIter_GetInnerStrideArray(iter)[0];                   \
        itemsize = NpyIter_GetDescrArray(iter)[1]->elsize;                    \
        innersizeptr = NpyIter_GetInnerLoopSizePtr(iter);                     \
        dataptrarray = NpyIter_GetDataPtrArray(iter);                         \
        if(PyArray_EquivTypes(PyArray_DESCR((PyArrayObject*) b), complex_int32_descr)) { \
            npy_intp i;                                                       \
            do {                                                              \
                npy_intp size = *innersizeptr;                                \
                src = dataptrarray[0];                                        \
                dst = dataptrarray[1];                                        \
                for(i = 0; i < size; i++, src += innerstride, dst += itemsize) {  \
                    *((complex_int32 *) dst) = complex_int32_##name(p, *((complex_int32 *) src)); \
                }                                                             \
            } while (iternext(iter));                                         \
        } else if(PyArray_ISINTEGER((PyArrayObject*) b)) {                      \
            npy_intp i;                                                       \
            do {                                                              \
                npy_intp size = *innersizeptr;                                \
                src = dataptrarray[0];                                        \
                dst = dataptrarray[1];                                        \
                for(i = 0; i < size; i++, src += innerstride, dst += itemsize) {  \
                    *(complex_int32 *) dst = complex_int32_##name##_scalar(p, *((npy_long *) src)); \
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
    pycomplexint32_##fake_name(PyObject* a, PyObject* b) {                     \
        /* PyObject *a_type, *a_repr, *b_type, *b_repr, *a_repr2, *b_repr2;   \ */ \
        /* char* a_char, b_char, a_char2, b_char2;                            \ */ \
        npy_int64 val64;                                                      \
        npy_int32 val32;                                                      \
        complex_int32 p = {0, 0};                                                 \
        if(PyArray_Check(b)) { return pycomplexint32_##fake_name##_array_operator(a, b); } \
        if(PyInt_Check(a) && PyComplexInt32_Check(b)) {                      \
            return PyComplexInt32_FromComplexInt32(complex_int32_scalar_##name(PyInt_AsLong(a), ((PyComplexInt32*)b)->obval));                                     \
        }                                                                  \
        PyComplexInt32_AsComplexInt32(p, a);                                    \
        if(PyComplexInt32_Check(b)) {                                          \
            return PyComplexInt32_FromComplexInt32(complex_int32_##name(p,((PyComplexInt32*)b)->obval)); \
        } else if(PyInt_Check(b)) {                                         \
            return PyComplexInt32_FromComplexInt32(complex_int32_##name##_scalar(p,PyInt_AsLong(b))); \
        }                                                                     \
        PyErr_SetString(PyExc_TypeError, "Binary operation involving complex_int32 and neither integer nor complex_int32.");                                                      \
        return NULL;                                                          \
    }
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(add, add)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(subtract, subtract)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(multiply, multiply)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(divide, divide)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(true_divide, divide)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(floor_divide, divide)
/* CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_RETURNER(power, power) */

#define CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_INPLACE(fake_name, name)        \
    static PyObject*                                                         \
    pycomplexint32_inplace_##fake_name(PyObject* a, PyObject* b) {            \
        complex_int32* p = {0};                                               \
        /* fprintf (stderr, "file %s, line %d, pycomplexint32_inplace_"#fake_name"(PyObject* a, PyObject* b).\n", __FILE__, __LINE__); \ */                                       \
        if(PyFloat_Check(a) || PyInt_Check(a)) {                             \
            PyErr_SetString(PyExc_TypeError, "Cannot in-place "#fake_name" a scalar by a complex_int32; should be handled by python.");                                         \
            return NULL;                                                     \
        }                                                                    \
        PyComplexInt32_AsComplexInt32Pointer(p, a);                            \
        if(PyComplexInt32_Check(b)) {                                         \
            complex_int32_inplace_##name(p,((PyComplexInt32*)b)->obval);       \
            Py_INCREF(a);                                                    \
            return a;                                                        \
        } else if(PyInt_Check(b)) {                                        \
            complex_int32_inplace_##name##_scalar(p,PyInt_AsLong(b));     \
            Py_INCREF(a);                                                    \
            return a;                                                        \
        }                                                                    \
        PyErr_SetString(PyExc_TypeError, "Binary in-place operation involving complex_int32 and neither complex nor complex_int32.");                                                 \
        return NULL;                                                         \
    }
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_INPLACE(add, add)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_INPLACE(subtract, subtract)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_INPLACE(multiply, multiply)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_INPLACE(divide, divide)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_INPLACE(true_divide, divide)
CI32CI32_CI32S_SCI32_BINARY_COMPLEX_INT32_INPLACE(floor_divide, divide)

static PyObject* pycomplexint32__reduce(PyComplexInt32* self) {
    return Py_BuildValue("O(OO)", Py_TYPE(self), PyInt_FromLong(self->obval.real), PyInt_FromLong(self->obval.imag));
}

static PyObject* pycomplexint32_getstate(PyComplexInt32* self, PyObject* args) {
    if( !PyArg_ParseTuple(args, ":getstate") ) {
        return NULL;
    }
    return Py_BuildValue("OO", PyInt_FromLong(self->obval.real),  PyInt_FromLong(self->obval.imag));
}

static PyObject* pycomplexint32_setstate(PyComplexInt32* self, PyObject* args) {
    complex_int32* c;
    c = &(self->obval);
    
    if( !PyArg_ParseTuple(args, "bb:setstate", &c->real, &c->imag) ) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// This is an array of methods (member functions) that will be
// available to use on the complex_int32 objects in python.  This is
// packaged up here, and will be used in the `tp_methods` field when
// definining the PyComplexInt32_Type below.
PyMethodDef pycomplexint32_methods[] = {
  // Unary bool returners
  {"nonzero", pycomplexint32_nonzero, METH_NOARGS,
   "True if the complex_int32 has all zero components"},
  {"isnan", pycomplexint32_isnan, METH_NOARGS,
   "True if the complex_int32 has any NAN components"},
  {"isinf", pycomplexint32_isinf, METH_NOARGS,
   "True if the complex_int32 has any INF components"},
  {"isfinite", pycomplexint32_isfinite, METH_NOARGS,
   "True if the complex_int32 has all finite components"},

  // Binary bool returners
  {"equal", pycomplexint32_equal, METH_O,
   "True if the complex_int32s are PRECISELY equal"},
  {"not_equal", pycomplexint32_not_equal, METH_O,
   "True if the complex_int32s are not PRECISELY equal"},
  {"less", pycomplexint32_less, METH_O,
   "Strict dictionary ordering"},
  {"greater", pycomplexint32_greater, METH_O,
   "Strict dictionary ordering"},
  {"less_equal", pycomplexint32_less_equal, METH_O,
   "Dictionary ordering"},
  {"greater_equal", pycomplexint32_greater_equal, METH_O,
   "Dictionary ordering"},
  
  // Unary ing returners
  {"real_part", pycomplexint32_real, METH_NOARGS,
   "real part of the complex value"},
  {"imag_part", pycomplexint32_imag, METH_NOARGS,
   "imaginary part of the complex value"},
  
  // Unary float returners
  {"absolute", pycomplexint32_absolute, METH_NOARGS,
   "Absolute value of complex_int32"},
  {"abs", pycomplexint32_absolute, METH_NOARGS,
   "Absolute value (Euclidean norm) of complex_int32"},
  {"norm", pycomplexint32_norm, METH_NOARGS,
   "Cayley norm (square of the absolute value) of complex_int32"},
  {"angle", pycomplexint32_angle, METH_NOARGS,
   "Angle through which rotor rotates"},

  // Unary complex_int32 returners
  // {"negative", pycomplexint32_negative, METH_NOARGS,
  //  "Return the negated complex_int32"},
  // {"positive", pycomplexint32_positive, METH_NOARGS,
  //  "Return the complex_int32 itself"},
  {"conjugate", pycomplexint32_conjugate, METH_NOARGS,
   "Return the complex conjugate of the complex_int32"},
  {"conj", pycomplexint32_conjugate, METH_NOARGS,
   "Return the complex conjugate of the complex_int32"},
  
  // complex_int32-complex_int32 or complex_int32-scalar binary complex_int32 returners
  // {"multiply", pycomplexint32_multiply, METH_O,
  //  "Standard (geometric) complex_int32 product"},
  // {"divide", pycomplexint32_divide, METH_O,
  //  "Standard (geometric) complex_int32 division"},
  // {"power", pycomplexint32_power, METH_O,
  //  "q.power(p) = (q.log() * p).exp()"},

  {"__reduce__", (PyCFunction)pycomplexint32__reduce, METH_NOARGS,
   "Return state information for pickling."},
  {"__getstate__", (PyCFunction)pycomplexint32_getstate, METH_VARARGS,
   "Return state information for pickling."},
  {"__setstate__", (PyCFunction)pycomplexint32_setstate, METH_VARARGS,
   "Reconstruct state information from pickle."},

  {NULL, NULL, 0, NULL}
};

static PyObject* pycomplexint32_num_negative(PyObject* a) { return pycomplexint32_negative(a,NULL); }
static PyObject* pycomplexint32_num_positive(PyObject* a) { return pycomplexint32_positive(a,NULL); }
static PyObject* pycomplexint32_num_absolute(PyObject* a) { return pycomplexint32_absolute(a,NULL); }
static int pycomplexint32_num_nonzero(PyObject* a) {
  complex_int32 q = ((PyComplexInt32*)a)->obval;
  return complex_int32_nonzero(q);
}
#define CANNOT_CONVERT_CI32(target)                                          \
  static PyObject* pycomplexint32_convert_##target(PyObject* a) {         \
    PyErr_SetString(PyExc_TypeError, "Cannot convert complex_int32 to " #target); \
    return NULL;                                                        \
  }
CANNOT_CONVERT_CI32(int)
CANNOT_CONVERT_CI32(float)
#if PY_MAJOR_VERSION < 3
CANNOT_CONVERT_CI32(long)
CANNOT_CONVERT_CI32(oct)
CANNOT_CONVERT_CI32(hex)
#endif

static PyNumberMethods pycomplexint32_as_number = {
  pycomplexint32_add,               // nb_add
  pycomplexint32_subtract,          // nb_subtract
  pycomplexint32_multiply,          // nb_multiply
  #if PY_MAJOR_VERSION < 3
  pycomplexint32_divide,            // nb_divide
  #endif
  0,                              // nb_remainder
  0,                              // nb_divmod
  0,                              // nb_power
  pycomplexint32_num_negative,      // nb_negative
  pycomplexint32_num_positive,      // nb_positive
  pycomplexint32_num_absolute,      // nb_absolute
  pycomplexint32_num_nonzero,       // nb_nonzero
  0,                              // nb_invert
  0,                              // nb_lshift
  0,                              // nb_rshift
  0,                              // nb_and
  0,                              // nb_xor
  0,                              // nb_or
  #if PY_MAJOR_VERSION < 3
  0,                              // nb_coerce
  #endif
  pycomplexint32_convert_int,       // nb_int
  #if PY_MAJOR_VERSION >= 3
  0,                              // nb_reserved
  #else
  pycomplexint32_convert_long,      // nb_long
  #endif
  pycomplexint32_convert_float,     // nb_float
  #if PY_MAJOR_VERSION < 3
  pycomplexint32_convert_oct,       // nb_oct
  pycomplexint32_convert_hex,       // nb_hex
  #endif
  pycomplexint32_inplace_add,       // nb_inplace_add
  pycomplexint32_inplace_subtract,  // nb_inplace_subtract
  pycomplexint32_inplace_multiply,  // nb_inplace_multiply
  #if PY_MAJOR_VERSION < 3
  pycomplexint32_inplace_divide,    // nb_inplace_divide
  #endif
  0,                              // nb_inplace_remainder
  0,                              // nb_inplace_power
  0,                              // nb_inplace_lshift
  0,                              // nb_inplace_rshift
  0,                              // nb_inplace_and
  0,                              // nb_inplace_xor
  0,                              // nb_inplace_or
  pycomplexint32_divide,            // nb_floor_divide
  pycomplexint32_divide,            // nb_true_divide
  pycomplexint32_inplace_divide,    // nb_inplace_floor_divide
  pycomplexint32_inplace_divide,    // nb_inplace_true_divide
  0,                              // nb_index
  #if PY_MAJOR_VERSION >= 3
  #if PY_MINOR_VERSION >= 5
  0,                              // nb_matrix_multiply
  0,                              // nb_inplace_matrix_multiply
  #endif
  #endif
};

// This is an array of members (member data) that will be available to
// use on the complex_int32 objects in python.  This is packaged up here,
// and will be used in the `tp_members` field when definining the
// PyComplexInt32_Type below.
PyMemberDef pycomplexint32_members[] = {
  {"real", T_SHORT, offsetof(PyComplexInt32, obval.real), 0,
   "The packed real component of the complex_int32"},
  {"imag", T_SHORT, offsetof(PyComplexInt32, obval.imag), 0,
   "The packed imaginary component of the complex_int32"},
  {NULL, 0, 0, 0, NULL}
};

static PyObject* pycomplexint32_richcompare(PyObject* a, PyObject* b, int op) {
    complex_int32 x = {0, 0};
    complex_int32 y = {0, 0};
    int result = 0;
    PyComplexInt32_AsComplexInt32(x,a);
    PyComplexInt32_AsComplexInt32(y,b);
    #define COMPARISONOP(py,op) case py: result = complex_int32_##op(x,y); break;
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

static long pycomplexint32_hash(PyObject *o) {
    complex_int32 q = ((PyComplexInt32 *)o)->obval;
    long value = 0x4567329;
    value = (10000004 * value) ^ _Py_HashDouble(q.real);
    value = (10000004 * value) ^ _Py_HashDouble(q.imag);
    if (value == -1) {
        value = -2;
    }
    return value;
}

static PyObject* pycomplexint32_repr(PyObject *o) {
    char str[1232];
    complex_int32 c = ((PyComplexInt32 *)o)->obval;
    sprintf(str, "complex_int32(%i, %i)", c.real, c.imag);
    return PyUString_FromString(str);
}

static PyObject* pycomplexint32_str(PyObject *o) {
    char str[64];
    complex_int32 c = ((PyComplexInt32 *)o)->obval;
    sprintf(str, "%i%+ij", c.real, c.imag);
    return PyUString_FromString(str);
}

// This establishes the complex_int32 as a python object (not yet a numpy
// scalar type).  The name may be a little counterintuitive; the idea
// is that this will be a type that can be used as an array dtype.
// Note that many of the slots below will be filled later, after the
// corresponding functions are defined.
static PyTypeObject PyComplexInt32_Type = {
#if PY_MAJOR_VERSION >= 3
  PyVarObject_HEAD_INIT(NULL, 0)
#else
  PyObject_HEAD_INIT(NULL)
  0,                                          // ob_size
#endif
  "complex_int32.complex_int32",                // tp_name
  sizeof(PyComplexInt32),                      // tp_basicsize
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
  pycomplexint32_repr,                         // tp_repr
  &pycomplexint32_as_number,                   // tp_as_number
  0,                                          // tp_as_sequence
  0,                                          // tp_as_mapping
  pycomplexint32_hash,                         // tp_hash
  0,                                          // tp_call
  pycomplexint32_str,                          // tp_str
  0,                                          // tp_getattro
  0,                                          // tp_setattro
  0,                                          // tp_as_buffer
#if PY_MAJOR_VERSION >= 3
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   // tp_flags
#else
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES, // tp_flags
#endif
  "Complex integer complex_int32 numbers",     // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  pycomplexint32_richcompare,                  // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  pycomplexint32_methods,                      // tp_methods
  pycomplexint32_members,                      // tp_members
  0,                                          // tp_getset
  0,                                          // tp_base; will be reset to &PyGenericArrType_Type after numpy import
  0,                                          // tp_dict
  0,                                          // tp_descr_get
  0,                                          // tp_descr_set
  0,                                          // tp_dictoffset
  pycomplexint32_init,                         // tp_init
  0,                                          // tp_alloc
  pycomplexint32_new,                          // tp_new
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
static PyArray_ArrFuncs _PyComplexInt32_ArrFuncs;
PyArray_Descr *complex_int32_descr;

static npy_bool CI32_nonzero (char *ip, PyArrayObject *ap) {
  complex_int32 c;
  complex_int32 zero = {0, 0};
  if (ap == NULL || PyArray_ISBEHAVED_RO(ap)) {
    c = *(complex_int32 *)ip;
  }
  else {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT16);
    descr->f->copyswap(&c.real, ip, !PyArray_ISNOTSWAPPED(ap), NULL);
    descr->f->copyswap(&c.imag, ip+1, !PyArray_ISNOTSWAPPED(ap), NULL);
    Py_DECREF(descr);
  }
  return (npy_bool) !complex_int32_equal(c, zero);
}

static void CI32_copyswap(complex_int32 *dst, complex_int32 *src, int swap, void *NPY_UNUSED(arr)) {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT16);
    descr->f->copyswapn(dst, sizeof(short int), src, sizeof(short int), 2, swap, NULL);
    Py_DECREF(descr);
}

static void CI32_copyswapn(complex_int32 *dst, npy_intp dstride,
                               complex_int32 *src, npy_intp sstride,
                               npy_intp n, int swap, void *NPY_UNUSED(arr)) {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT16);
    descr->f->copyswapn(&dst->real, dstride, &src->real, sstride, n, swap, NULL);
    descr->f->copyswapn(&dst->imag, dstride, &src->imag, sstride, n, swap, NULL);
    Py_DECREF(descr);    
}

static int CI32_setitem(PyObject* item, complex_int32* c, void* NPY_UNUSED(ap)) {
    PyObject *element;
    if( PyComplexInt32_Check(item) ) {
        memcpy(c, &(((PyComplexInt32 *)item)->obval), sizeof(complex_int32));
    } else if( PyComplexInt8_Check(item) ) {
        complex_int8 i8;
        PyComplexInt8_AsComplexInt8(i8, item);
        const signed char* sc = fourBitLUT[i8.real_imag];
        c->real = sc[0];
        c->imag = sc[1];
    } else if( PyComplexInt16_Check(item) ) {
        complex_int16 i16;
        PyComplexInt16_AsComplexInt16(i16, item);
        c->real = i16.real;
        c->imag = i16.imag;
    } else if( PyComplex_Check(item) ) {
        c->real = PyComplex_RealAsDouble(item);
        c->imag = PyComplex_ImagAsDouble(item);
    } else if( PyInt_Check(item) ) {
        c->real = PyInt_AsLong(item);
        c->imag = 0;
    } else if( PySequence_Check(item) && PySequence_Length(item)==2 ) {
        element = PySequence_GetItem(item, 0);
        if(element == NULL) {
            return -1;
        } /* Not a sequence, or other failure */
        c->real = PyInt_AsLong(element);
        Py_DECREF(element);
        
        element = PySequence_GetItem(item, 1);
        if(element == NULL) {
            return -1;
        } /* Not a sequence, or other failure */
        c->imag = PyInt_AsLong(element);
        Py_DECREF(element);
    } else {
        PyErr_SetString(PyExc_TypeError, "Unknown input to CI32_setitem");
        return -1;
    }
    return 0;
}

// When a numpy array of dtype=complex_int32 is indexed, this function is
// called, returning a new complex_int32 object with a copy of the
// data... sometimes...
static PyObject* CI32_getitem(void* data, void* NPY_UNUSED(arr)) {
    complex_int32 q;
    memcpy(&q, data, sizeof(complex_int32));
    return PyComplexInt32_FromComplexInt32(q);
}

static int CI32_compare(complex_int32 *pa, complex_int32 *pb, PyArrayObject *NPY_UNUSED(ap)) {
    complex_int32 a = *pa, b = *pb;
    npy_bool anan, bnan;
    int ret;
    
    anan = complex_int32_isnan(a);
    bnan = complex_int32_isnan(b);
    
    if( anan ) {
        ret = bnan ? 0 : -1;
    } else if( bnan ) {
        ret = 1;
    } else if( complex_int32_less(a, b) ) {
        ret = -1;
    } else if( complex_int32_less(b, a) ) {
        ret = 1;
    } else {
        ret = 0;
    }
    
    return ret;
}

static int CI32_argmax(complex_int32 *ip, npy_intp n, npy_intp *max_ind, PyArrayObject *NPY_UNUSED(aip)) {
    npy_intp i;
    complex_int32 mp = *ip;
    
    *max_ind = 0;
    
    if( complex_int32_isnan(mp) ) {
        /* nan encountered; it's maximal */
        return 0;
    }

    for(i = 1; i < n; i++) {
        ip++;
        /*
         * Propagate nans, similarly as max() and min()
         */
        if( !(complex_int32_less_equal(*ip, mp)) ) {  /* negated, for correct nan handling */
            mp = *ip;
            *max_ind = i;
            if (complex_int32_isnan(mp)) {
                /* nan encountered, it's maximal */
                break;
            }
        }
    }
    return 0;
}

static void CI32_fillwithscalar(complex_int32 *buffer, npy_intp length, complex_int32 *value, void *NPY_UNUSED(ignored)) {
    npy_intp i;
    complex_int32 val = *value;

    for (i = 0; i < length; ++i) {
        buffer[i] = val;
    }
}

// This is a macro (followed by applications of the macro) that cast
// the input types to standard complex_int32 with only a nonzero scalar
// part.
#define MAKE_T_TO_CI32(TYPE, type)                                      \
static void TYPE ## _to_complex_int32(type *ip, complex_int32 *op, npy_intp n, \
                                  PyArrayObject *NPY_UNUSED(aip),      \
                                  PyArrayObject *NPY_UNUSED(aop)) {    \
    while (n--) {                                                      \
        op->real = (signed char) (*ip++);                              \
        op->imag = (signed char) (*ip++);                              \
        *op++;                                                         \
    }                                                                  \
}

MAKE_T_TO_CI32(BOOL, npy_bool);
MAKE_T_TO_CI32(BYTE, npy_byte);
MAKE_T_TO_CI32(SHORT, npy_short);

static void CI8_to_complex_int32(complex_int8* ip, short int *op, npy_intp n, PyArrayObject *NPY_UNUSED(aip), PyArrayObject *NPY_UNUSED(aop)) { 
    const signed char *sc;
    while (n--) {
        sc = fourBitLUT[ip->real_imag];
        *(op++) = (short int) sc[0];
        *(op++) = (short int) sc[1];
        (*ip++);
    }
}

static void CI16_to_complex_int32(complex_int16* ip, short int *op, npy_intp n, PyArrayObject *NPY_UNUSED(aip), PyArrayObject *NPY_UNUSED(aop)) { 
    while (n--) {
        *(op++) = (short int) ip->real;
        *(op++) = (short int) ip->imag;
        (*ip++);
    }
}


// This is a macro (followed by applications of the macro) that cast
// the input complex types from complex_int32.
#define MAKE_CI32_TO_CT(TYPE, type)                                    \
static void complex_int32_to_## TYPE(complex_int32* ip, type *op, npy_intp n, \
                                PyArrayObject *NPY_UNUSED(aip),       \
                                PyArrayObject *NPY_UNUSED(aop)) {     \
    while (n--) {                                                     \
        *(op++) = (type) ip->real;                                    \
        *(op++) = (type) ip->imag;                                    \
        (*ip++);                                                      \
    }                                                                 \
}

MAKE_CI32_TO_CT(CFLOAT, npy_float);
MAKE_CI32_TO_CT(CDOUBLE, npy_double);
MAKE_CI32_TO_CT(CLONGDOUBLE, npy_longdouble);

static void register_cast_function_ci32(int sourceType, int destType, PyArray_VectorUnaryFunc *castfunc) {
    PyArray_Descr *descr = PyArray_DescrFromType(sourceType);
    PyArray_RegisterCastFunc(descr, destType, castfunc);
    if( ( (destType == NPY_CFLOAT) \
         || (destType == NPY_CDOUBLE) \
         || (destType == NPY_CLONGDOUBLE) ) ) {
        PyArray_RegisterCanCast(descr, destType, NPY_COMPLEX_SCALAR);
    } else {
        PyArray_RegisterCanCast(descr, destType, NPY_NOSCALAR);
    }
    Py_DECREF(descr);
}

// This is a macro that will be used to define the various basic unary
// complex_int32 functions, so that they can be applied quickly to a
// numpy array of complex_int32s.
#define UNARY_GEN_UFUNC_CI32(ufunc_name, func_name, ret_type)                \
    static void complex_int32_##ufunc_name##_ufunc(char** args,              \
                                                  npy_intp* dimensions,     \
                                                  npy_intp* steps,          \
                                                  void* NPY_UNUSED(data)) { \
    char *ip1 = args[0], *op1 = args[1];                                    \
    npy_intp is1 = steps[0], os1 = steps[1];                                \
    npy_intp n = dimensions[0];                                             \
    npy_intp i;                                                             \
    for(i = 0; i < n; i++, ip1 += is1, op1 += os1){                         \
      const complex_int32 in1 = *(complex_int32 *)ip1;                        \
      *((ret_type *)op1) = complex_int32_##func_name(in1);};}
#define UNARY_UFUNC_CI32(name, ret_type) \
  UNARY_GEN_UFUNC_CI32(name, name, ret_type)
// And these all do the work mentioned above, using the macro
UNARY_UFUNC_CI32(isnan, npy_bool)
UNARY_UFUNC_CI32(isinf, npy_bool)
UNARY_UFUNC_CI32(isfinite, npy_bool)
UNARY_UFUNC_CI32(real, npy_long)
UNARY_UFUNC_CI32(imag, npy_long)
UNARY_UFUNC_CI32(norm, npy_double)
UNARY_UFUNC_CI32(absolute, npy_double)
UNARY_UFUNC_CI32(angle, npy_double)
UNARY_UFUNC_CI32(negative, complex_int32)
UNARY_UFUNC_CI32(conjugate, complex_int32)
static void complex_int32_positive_ufunc(char** args, npy_intp* dimensions, npy_intp* steps, void* NPY_UNUSED(data)) {
    char *ip1 = args[0], *op1 = args[1];
    npy_intp is1 = steps[0], os1 = steps[1];
    npy_intp n = dimensions[0];
    npy_intp i;
    for(i = 0; i < n; i++, ip1 += is1, op1 += os1) {
        const complex_int32 in1 = *(complex_int32 *)ip1;
        *((complex_int32 *)op1) = in1;
    }
}

// This is a macro that will be used to define the various basic binary
// complex_int32 functions, so that they can be applied quickly to a
// numpy array of complex_int32s.
#define BINARY_GEN_UFUNC_CI32(ufunc_name, func_name, arg_type1, arg_type2, ret_type) \
  static void complex_int32_##ufunc_name##_ufunc(char** args,                     \
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
      *((ret_type *)op1) = complex_int32_##func_name(in1, in2);                   \
    };                                                                          \
  };
// A couple special-case versions of the above
#define BINARY_UFUNC_CI32(name, ret_type)                    \
  BINARY_GEN_UFUNC_CI32(name, name, complex_int32, complex_int32, ret_type)
#define BINARY_SCALAR_UFUNC_CI32(name, ret_type)                             \
  BINARY_GEN_UFUNC_CI32(name##_scalar, name##_scalar, complex_int32, npy_long, ret_type) \
  BINARY_GEN_UFUNC_CI32(scalar_##name, scalar_##name, npy_long, complex_int32, ret_type)
// And these all do the work mentioned above, using the macros
BINARY_UFUNC_CI32(add, complex_int32)
BINARY_UFUNC_CI32(subtract, complex_int32)
BINARY_UFUNC_CI32(multiply, complex_int32)
BINARY_UFUNC_CI32(divide, complex_int32)
BINARY_GEN_UFUNC_CI32(true_divide, divide, complex_int32, complex_int32, complex_int32)
BINARY_GEN_UFUNC_CI32(floor_divide, divide, complex_int32, complex_int32, complex_int32)
BINARY_UFUNC_CI32(equal, npy_bool)
BINARY_UFUNC_CI32(not_equal, npy_bool)
BINARY_UFUNC_CI32(less, npy_bool)
BINARY_UFUNC_CI32(less_equal, npy_bool)
BINARY_SCALAR_UFUNC_CI32(add, complex_int32)
BINARY_SCALAR_UFUNC_CI32(subtract, complex_int32)
BINARY_SCALAR_UFUNC_CI32(multiply, complex_int32)
BINARY_SCALAR_UFUNC_CI32(divide, complex_int32)
BINARY_GEN_UFUNC_CI32(true_divide_scalar, divide_scalar, complex_int32, long, complex_int32)
BINARY_GEN_UFUNC_CI32(floor_divide_scalar, divide_scalar, complex_int32, long, complex_int32)
BINARY_GEN_UFUNC_CI32(scalar_true_divide, scalar_divide, long, complex_int32, complex_int32)
BINARY_GEN_UFUNC_CI32(scalar_floor_divide, scalar_divide, long, complex_int32, complex_int32)

static PyObject* complex_int32_arrtype_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    complex_int32 c;
    Py_complex cmplx;

    if( !PyArg_ParseTuple(args, "D", &cmplx) ) {
        return NULL;
    }
    
    // Ouch
    c.real = cmplx.real;
    c.imag = cmplx.imag;
    return PyArray_Scalar(&c, complex_int32_descr, NULL);
}

static PyObject* gentype_richcompare_ci32(PyObject *self, PyObject *other, int cmp_op) {
    PyObject *arr, *ret;

    arr = PyArray_FromScalar(self, NULL);
    if (arr == NULL) {
        return NULL;
    }
    ret = Py_TYPE(arr)->tp_richcompare(arr, other, cmp_op);
    Py_DECREF(arr);
    return ret;
}

static long complex_int32_arrtype_hash(PyObject *o) {
    complex_int32 c = ((PyComplexInt32 *)o)->obval;
    long value = 0x4567329;
    value = (10000004 * value) ^ _Py_HashDouble(c.real);
    value = (10000004 * value) ^ _Py_HashDouble(c.imag);
    if( value == -1 ) {
        value = -2;
    }
    return value;
}

static PyObject* complex_int32_arrtype_repr(PyObject *o) {
    char str[64];
    complex_int32 c = ((PyComplexInt32 *)o)->obval;
    sprintf(str, "complex_int32(%i, %i)", c.real, c.imag);
    return PyUString_FromString(str);
}

static PyObject* complex_int32_arrtype_str(PyObject *o) {
    char str[64];
    complex_int32 c = ((PyComplexInt32 *)o)->obval;
    sprintf(str, "%i%+ij", c.real, c.imag);
    return PyUString_FromString(str);
}

int complex_int32_elsize = sizeof(complex_int32);

typedef struct {
    char c;
    complex_int32 q;
} align_test_ci32;
int complex_int32_alignment = offsetof(align_test_ci32, q);

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//                                                             //
//  Everything above was preparation for the following set up  //
//                                                             //
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

int create_complex_int32(PyObject* m, PyObject* numpy_dict) {
    int complexi32Num;
    PyObject *tmp_ufunc;
    int arg_types[3];
    
    /* Register the complex_int32 array base type. */
    PyComplexInt32_Type.tp_base = &PyGenericArrType_Type;
    if (PyType_Ready(&PyComplexInt32_Type) < 0) {
        return -2;
    }
    
    // The array functions, to be used below.  This InitArrFuncs
    // function is a convenient way to set all the fields to zero
    // initially, so we don't get undefined behavior.
    PyArray_InitArrFuncs(&_PyComplexInt32_ArrFuncs);
    _PyComplexInt32_ArrFuncs.nonzero = (PyArray_NonzeroFunc*)CI32_nonzero;
    _PyComplexInt32_ArrFuncs.copyswap = (PyArray_CopySwapFunc*)CI32_copyswap;
    _PyComplexInt32_ArrFuncs.copyswapn = (PyArray_CopySwapNFunc*)CI32_copyswapn;
    _PyComplexInt32_ArrFuncs.setitem = (PyArray_SetItemFunc*)CI32_setitem;
    _PyComplexInt32_ArrFuncs.getitem = (PyArray_GetItemFunc*)CI32_getitem;
    _PyComplexInt32_ArrFuncs.compare = (PyArray_CompareFunc*)CI32_compare;
    _PyComplexInt32_ArrFuncs.argmax = (PyArray_ArgFunc*)CI32_argmax;
    _PyComplexInt32_ArrFuncs.fillwithscalar = (PyArray_FillWithScalarFunc*)CI32_fillwithscalar;
    
    // The complex_int32 array descr
    complex_int32_descr = PyObject_New(PyArray_Descr, &PyArrayDescr_Type);
    complex_int32_descr->typeobj = &PyComplexInt32_Type;
    complex_int32_descr->kind = 'V';
    complex_int32_descr->type = 'i';
    complex_int32_descr->byteorder = '=';
    complex_int32_descr->flags = NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM;
    complex_int32_descr->type_num = 0; // assigned at registration
    complex_int32_descr->elsize = complex_int32_elsize;
    complex_int32_descr->alignment = complex_int32_alignment;
    complex_int32_descr->subarray = NULL;
    complex_int32_descr->fields = NULL;
    complex_int32_descr->names = NULL;
    complex_int32_descr->f = &_PyComplexInt32_ArrFuncs;
    complex_int32_descr->metadata = NULL;
    complex_int32_descr->c_metadata = NULL;
    
    Py_INCREF(&PyComplexInt32_Type);
    complexi32Num = PyArray_RegisterDataType(complex_int32_descr);
    
    if( complexi32Num < 0 || complexi32Num != NPY_COMPLEX_INT32 ) {
        return -1;
    }
    
    register_cast_function_ci32(NPY_BOOL, complexi32Num, (PyArray_VectorUnaryFunc*)BOOL_to_complex_int32);
    register_cast_function_ci32(NPY_BYTE, complexi32Num, (PyArray_VectorUnaryFunc*)BYTE_to_complex_int32);
    register_cast_function_ci32(NPY_SHORT, complexi32Num, (PyArray_VectorUnaryFunc*)SHORT_to_complex_int32);
    
    register_cast_function_ci32(NPY_COMPLEX_INT8, complexi32Num, (PyArray_VectorUnaryFunc*)CI8_to_complex_int32);
    register_cast_function_ci32(NPY_COMPLEX_INT16, complexi32Num, (PyArray_VectorUnaryFunc*)CI16_to_complex_int32);
    
    register_cast_function_ci32(complexi32Num, NPY_CFLOAT, (PyArray_VectorUnaryFunc*)complex_int32_to_CFLOAT);
    register_cast_function_ci32(complexi32Num, NPY_CDOUBLE, (PyArray_VectorUnaryFunc*)complex_int32_to_CDOUBLE);
    register_cast_function_ci32(complexi32Num, NPY_CLONGDOUBLE, (PyArray_VectorUnaryFunc*)complex_int32_to_CLONGDOUBLE);
    
    // These macros will be used below
#define REGISTER_UFUNC_CI32(name)\
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name),\
            complex_int32_descr->type_num, complex_int32_##name##_ufunc, arg_types, NULL)
    
#define REGISTER_SCALAR_UFUNC_CI32(name)\
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name),\
            complex_int32_descr->type_num, complex_int32_##name##_scalar_ufunc, arg_types, NULL)
    
#define REGISTER_UFUNC_SCALAR_CI32(name)                                   \
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name), \
                                complex_int32_descr->type_num, complex_int32_##name##_scalar_ufunc, arg_types, NULL)
    
#define REGISTER_NEW_UFUNC_GENERAL_CI32(pyname, cname, nargin, nargout, doc) \
    tmp_ufunc = PyUFunc_FromFuncAndData(NULL, NULL, NULL, 0, nargin, nargout, \
                                        PyUFunc_None, #pyname, doc, 0); \
    PyUFunc_RegisterLoopForType((PyUFuncObject *)tmp_ufunc,             \
                                complex_int32_descr->type_num, complex_int32_##cname##_ufunc, arg_types, NULL); \
    PyDict_SetItemString(numpy_dict, #pyname, tmp_ufunc);               \
    Py_DECREF(tmp_ufunc)

  #define REGISTER_NEW_UFUNC_CI32(name, nargin, nargout, doc)                \
    REGISTER_NEW_UFUNC_GENERAL_CI32(name, name, nargin, nargout, doc)
    
    /* complex_int32 -> bool */
    arg_types[0] = complex_int32_descr->type_num;
    arg_types[1] = NPY_BOOL;
    
    REGISTER_UFUNC_CI32(isnan);
    REGISTER_UFUNC_CI32(isinf);
    REGISTER_UFUNC_CI32(isfinite);
    
    /* complex_int32 -> long */
    arg_types[1] = NPY_LONG;
    
    REGISTER_NEW_UFUNC_GENERAL_CI32(real_part, real, 1, 1, \
                                    "Return the real part of the complex number\n");
    REGISTER_NEW_UFUNC_GENERAL_CI32(imag_part, imag, 1, 1, \
                                    "Return the imaginary part of the complex number\n");
    
    /* complex_int32 -> double */
    arg_types[1] = NPY_DOUBLE;
    
    REGISTER_NEW_UFUNC_CI32(norm, 1, 1, \
                           "Return Cayley norm (square of the absolute value) of each complex_int32\n");
    REGISTER_UFUNC_CI32(absolute);
    REGISTER_NEW_UFUNC_GENERAL_CI32(angle, angle, 1, 1, \
                                   "Return the angle of the complex argument\n");
    
    /* complex_int32 -> complex_int32 */
    arg_types[1] = complex_int32_descr->type_num;
    
    REGISTER_UFUNC_CI32(negative);
    REGISTER_UFUNC_CI32(conjugate);

    /* complex_int32, complex_int32 -> bool */
    arg_types[2] = NPY_BOOL;

    REGISTER_UFUNC_CI32(equal);
    REGISTER_UFUNC_CI32(not_equal);
    REGISTER_UFUNC_CI32(less);
    REGISTER_UFUNC_CI32(less_equal);
    
    /* complex_int32, complex_int32 -> complex_int32 */
    arg_types[2] = complex_int32_descr->type_num;
    
    REGISTER_UFUNC_CI32(add);
    REGISTER_UFUNC_CI32(subtract);
    REGISTER_UFUNC_CI32(multiply);
    REGISTER_UFUNC_CI32(divide);
    REGISTER_UFUNC_CI32(true_divide);
    REGISTER_UFUNC_CI32(floor_divide);
    
    /* long, complex_int32 -> complex_int32 */
    arg_types[0] = NPY_LONG;
    arg_types[1] = complex_int32_descr->type_num;
    arg_types[2] = complex_int32_descr->type_num;
    
    REGISTER_SCALAR_UFUNC_CI32(add);
    REGISTER_SCALAR_UFUNC_CI32(subtract);
    REGISTER_SCALAR_UFUNC_CI32(multiply);
    REGISTER_SCALAR_UFUNC_CI32(divide);
    REGISTER_SCALAR_UFUNC_CI32(true_divide);
    REGISTER_SCALAR_UFUNC_CI32(floor_divide);
    
    /* complex_int32, long -> complex_int32 */
    arg_types[0] = complex_int32_descr->type_num;
    arg_types[1] = NPY_LONG;
    arg_types[2] = complex_int32_descr->type_num;
    
    REGISTER_UFUNC_SCALAR_CI32(add);
    REGISTER_UFUNC_SCALAR_CI32(subtract);
    REGISTER_UFUNC_SCALAR_CI32(multiply);
    REGISTER_UFUNC_SCALAR_CI32(divide);
    REGISTER_UFUNC_SCALAR_CI32(true_divide);
    REGISTER_UFUNC_SCALAR_CI32(floor_divide);
    
    PyModule_AddObject(m, "complex_int32", (PyObject *)&PyComplexInt32_Type);
    
    return complexi32Num;
}
