#include "complex_int16.h"

typedef struct {
    PyObject_HEAD
    complex_int16 obval;
} PyComplexInt16;

static PyTypeObject PyComplexInt16_Type;

PyArray_Descr *complex_int16_descr;

static NPY_INLINE int PyComplexInt16_Check(PyObject *o) {
    return PyObject_IsInstance(o, (PyObject*) &PyComplexInt16_Type);
}

static PyObject* PyComplexInt16_FromComplexInt16(complex_int16 c) {
    PyComplexInt16 *p = (PyComplexInt16 *) PyComplexInt16_Type.tp_alloc(&PyComplexInt16_Type, 0);
    if( p != NULL ) {
        p->obval = c;
    }
    return (PyObject *) p;
}

#define PyComplexInt16_AsComplexInt16(c, o)                     \
    if( PyComplexInt16_Check(o) ) {                             \
        c = ((PyComplexInt16 *) o)->obval;                      \
    } else {                                                    \
        PyErr_SetString(PyExc_TypeError,                        \
                        "Input object is not a complex_int16"); \
        return NULL;                                            \
    }                                                           \

#define PyComplexInt16_AsComplexInt16Pointer(c, o)              \
    if( PyComplexInt16_Check(o) ) {                             \
        c = &((PyComplexInt16 *) o)->obval;                     \
    } else {                                                    \
        PyErr_SetString(PyExc_TypeError,                        \
                        "Input object is not a complex_int16"); \
        return NULL;                                            \
    }    

static PyObject* pycomplexint16_new(PyTypeObject *type, PyObject *NPY_UNUSED(args), PyObject *NPY_UNUSED(kwds)) {
    PyComplexInt16 *self;
    self = (PyComplexInt16 *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

static int pycomplexint16_init(PyObject *self, PyObject *args, PyObject *kwds) {
    Py_ssize_t size = PyTuple_Size(args);
    complex_int16 *c;
    Py_complex cmplx;
    long real, imag;
    
    PyObject *C = {0};
    c = &(((PyComplexInt16 *) self)->obval);
    
    if( kwds != NULL && PyDict_Size(kwds) ) {
        PyErr_SetString(PyExc_TypeError,
                        "complex_int16 constructor takes no keyword arguments");
        return -1;
    }
    
    c->real = 0;
    c->imag = 0;
    if( size == 0 ) {
        return 0;
    } else if( size == 1) {
        if( PyArg_ParseTuple(args, "O", &C) && (PyComplexInt8_Check(C) || PyComplexInt16_Check(C)) ) {
            if( PyComplexInt8_Check(C) ) {
                const signed char* sc = fourBitLUT[((PyComplexInt8 *) C)->obval.real_imag];
                c->real = sc[0];
                c->imag = sc[1];
                return 0;
            } else if( PyComplexInt16_Check(C) ) {
                c->real = ((PyComplexInt16 *) C)->obval.real;
                c->imag = ((PyComplexInt16 *) C)->obval.imag;
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
                    "complex_int16 constructor takes zero or two int16 values, or a single argument of a complex_int16 or Python complex");
    return -1;
}

#define UNARY_BOOL_RETURNER_CI16(name)                                       \
    static PyObject*                                                    \
    pycomplexint16_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {       \
        complex_int16 c = {0, 0};                                          \
        PyComplexInt16_AsComplexInt16(c, a);                            \
        return PyBool_FromLong(complex_int16_##name(c));                \
    }
UNARY_BOOL_RETURNER_CI16(nonzero)
UNARY_BOOL_RETURNER_CI16(isnan)
UNARY_BOOL_RETURNER_CI16(isinf)
UNARY_BOOL_RETURNER_CI16(isfinite)

#define BINARY_BOOL_RETURNER_CI16(name)                                      \
    static PyObject*                                                    \
    pycomplexint16_##name(PyObject* a, PyObject* b) {                   \
        complex_int16 p = {0, 0};                                          \
        complex_int16 q = {0, 0};                                          \
        PyComplexInt16_AsComplexInt16(p, a);                            \
        PyComplexInt16_AsComplexInt16(q, b);                            \
        return PyBool_FromLong(complex_int16_##name(p,q));              \
    }
BINARY_BOOL_RETURNER_CI16(equal)
BINARY_BOOL_RETURNER_CI16(not_equal)
BINARY_BOOL_RETURNER_CI16(less)
BINARY_BOOL_RETURNER_CI16(greater)
BINARY_BOOL_RETURNER_CI16(less_equal)
BINARY_BOOL_RETURNER_CI16(greater_equal)

#define UNARY_INT_RETURNER_CI16(name)                                      \
    static PyObject*                                                    \
        pycomplexint16_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int16 q = {0, 0};                                           \
        PyComplexInt16_AsComplexInt16(q, a);                              \
        return PyInt_FromLong(complex_int16_##name(q));              \
    }
UNARY_INT_RETURNER_CI16(real)
UNARY_INT_RETURNER_CI16(imag)

#define UNARY_FLOAT_RETURNER_CI16(name)                                      \
    static PyObject*                                                    \
        pycomplexint16_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int16 q = {0, 0};                                           \
        PyComplexInt16_AsComplexInt16(q, a);                              \
        return PyFloat_FromDouble(complex_int16_##name(q));              \
    }
UNARY_FLOAT_RETURNER_CI16(absolute)
UNARY_FLOAT_RETURNER_CI16(norm)
UNARY_FLOAT_RETURNER_CI16(angle)

#define UNARY_COMPLEX_INT16_RETURNER_CI16(name)                               \
    static PyObject*                                                    \
        pycomplexint16_##name(PyObject* a, PyObject* NPY_UNUSED(b)) {    \
        complex_int16 q = {0, 0};                          \
        PyComplexInt16_AsComplexInt16(q, a);                              \
        return PyComplexInt16_FromComplexInt16(complex_int16_##name(q));   \
    }
UNARY_COMPLEX_INT16_RETURNER_CI16(negative)
UNARY_COMPLEX_INT16_RETURNER_CI16(conjugate)

static PyObject* pycomplexint16_positive(PyObject* self, PyObject* NPY_UNUSED(b)) {
    Py_INCREF(self);
    return self;
}

#define CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(fake_name, name)        \
    static PyObject*                                                          \
    pycomplexint16_##fake_name##_array_operator(PyObject* a, PyObject* b) {    \
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
        complex_int16 p = {0, 0};                                                 \
        PyComplexInt16_AsComplexInt16(p, a);                                    \
        flags = NPY_ITER_EXTERNAL_LOOP;                                       \
        op[0] = (PyArrayObject *) b;                                          \
        op[1] = NULL;                                                         \
        op_flags[0] = NPY_ITER_READONLY;                                      \
        op_flags[1] = NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE;                 \
        op_dtypes[0] = PyArray_DESCR((PyArrayObject*) b);                     \
        op_dtypes[1] = complex_int16_descr;                                    \
        iter = NpyIter_MultiNew(2, op, flags, NPY_KEEPORDER, NPY_NO_CASTING, op_flags, op_dtypes); \
        if (iter == NULL) {                                                   \
            return NULL;                                                      \
        }                                                                     \
        iternext = NpyIter_GetIterNext(iter, NULL);                           \
        innerstride = NpyIter_GetInnerStrideArray(iter)[0];                   \
        itemsize = NpyIter_GetDescrArray(iter)[1]->elsize;                    \
        innersizeptr = NpyIter_GetInnerLoopSizePtr(iter);                     \
        dataptrarray = NpyIter_GetDataPtrArray(iter);                         \
        if(PyArray_EquivTypes(PyArray_DESCR((PyArrayObject*) b), complex_int16_descr)) { \
            npy_intp i;                                                       \
            do {                                                              \
                npy_intp size = *innersizeptr;                                \
                src = dataptrarray[0];                                        \
                dst = dataptrarray[1];                                        \
                for(i = 0; i < size; i++, src += innerstride, dst += itemsize) {  \
                    *((complex_int16 *) dst) = complex_int16_##name(p, *((complex_int16 *) src)); \
                }                                                             \
            } while (iternext(iter));                                         \
        } else if(PyArray_ISINTEGER((PyArrayObject*) b)) {                      \
            npy_intp i;                                                       \
            do {                                                              \
                npy_intp size = *innersizeptr;                                \
                src = dataptrarray[0];                                        \
                dst = dataptrarray[1];                                        \
                for(i = 0; i < size; i++, src += innerstride, dst += itemsize) {  \
                    *(complex_int16 *) dst = complex_int16_##name##_scalar(p, *((npy_long *) src)); \
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
    pycomplexint16_##fake_name(PyObject* a, PyObject* b) {                     \
        /* PyObject *a_type, *a_repr, *b_type, *b_repr, *a_repr2, *b_repr2;   \ */ \
        /* char* a_char, b_char, a_char2, b_char2;                            \ */ \
        npy_int64 val64;                                                      \
        npy_int32 val32;                                                      \
        complex_int16 p = {0, 0};                                                 \
        if(PyArray_Check(b)) { return pycomplexint16_##fake_name##_array_operator(a, b); } \
        if(PyInt_Check(a) && PyComplexInt16_Check(b)) {                      \
            return PyComplexInt16_FromComplexInt16(complex_int16_scalar_##name(PyInt_AsLong(a), ((PyComplexInt16*)b)->obval));                                     \
        }                                                                     \
        PyComplexInt16_AsComplexInt16(p, a);                                    \
        if(PyComplexInt16_Check(b)) {                                          \
            return PyComplexInt16_FromComplexInt16(complex_int16_##name(p,((PyComplexInt16*)b)->obval)); \
        } else if(PyInt_Check(b)) {                                         \
            return PyComplexInt16_FromComplexInt16(complex_int16_##name##_scalar(p,PyInt_AsLong(b))); \
        }                                                                     \
        PyErr_SetString(PyExc_TypeError, "Binary operation involving complex_int16 and neither integer nor complex_int16.");                                                      \
        return NULL;                                                          \
    }
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(add, add)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(subtract, subtract)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(multiply, multiply)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(divide, divide)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(true_divide, divide)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(floor_divide, divide)
/* CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_RETURNER(power, power) */

#define CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_INPLACE(fake_name, name)        \
    static PyObject*                                                         \
    pycomplexint16_inplace_##fake_name(PyObject* a, PyObject* b) {            \
        complex_int16* p = {0};                                               \
        /* fprintf (stderr, "file %s, line %d, pycomplexint16_inplace_"#fake_name"(PyObject* a, PyObject* b).\n", __FILE__, __LINE__); \ */                                       \
        if(PyFloat_Check(a) || PyInt_Check(a)) {                             \
            PyErr_SetString(PyExc_TypeError, "Cannot in-place "#fake_name" a scalar by a complex_int16; should be handled by python.");                                         \
            return NULL;                                                     \
        }                                                                    \
        PyComplexInt16_AsComplexInt16Pointer(p, a);                            \
        if(PyComplexInt16_Check(b)) {                                         \
            complex_int16_inplace_##name(p,((PyComplexInt16*)b)->obval);       \
            Py_INCREF(a);                                                    \
            return a;                                                        \
        } else if(PyInt_Check(b)) {                                        \
            complex_int16_inplace_##name##_scalar(p,PyInt_AsLong(b));     \
            Py_INCREF(a);                                                    \
            return a;                                                        \
        }                                                                    \
        PyErr_SetString(PyExc_TypeError, "Binary in-place operation involving complex_int16 and neither integer nor complex_int16.");                                                 \
        return NULL;                                                         \
    }
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_INPLACE(add ,add)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_INPLACE(subtract, subtract)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_INPLACE(multiply, multiply)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_INPLACE(divide, divide)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_INPLACE(true_divide, divide)
CI16CI16_CI16S_SCI16_BINARY_COMPLEX_INT16_INPLACE(floor_divide, divide)

static PyObject* pycomplexint16__reduce(PyComplexInt16* self) {
    return Py_BuildValue("O(OO)", Py_TYPE(self), PyInt_FromLong(self->obval.real), PyInt_FromLong(self->obval.imag));
}

static PyObject* pycomplexint16_getstate(PyComplexInt16* self, PyObject* args) {
    if( !PyArg_ParseTuple(args, ":getstate") ) {
        return NULL;
    }
    return Py_BuildValue("OO", PyInt_FromLong(self->obval.real),  PyInt_FromLong(self->obval.imag));
}

static PyObject* pycomplexint16_setstate(PyComplexInt16* self, PyObject* args) {
    complex_int16* c;
    c = &(self->obval);
    
    if( !PyArg_ParseTuple(args, "bb:setstate", &c->real, &c->imag) ) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// This is an array of methods (member functions) that will be
// available to use on the complex_int16 objects in python.  This is
// packaged up here, and will be used in the `tp_methods` field when
// definining the PyComplexInt16_Type below.
PyMethodDef pycomplexint16_methods[] = {
  // Unary bool returners
  {"nonzero", pycomplexint16_nonzero, METH_NOARGS,
   "True if the complex_int16 has all zero components"},
  {"isnan", pycomplexint16_isnan, METH_NOARGS,
   "True if the complex_int16 has any NAN components"},
  {"isinf", pycomplexint16_isinf, METH_NOARGS,
   "True if the complex_int16 has any INF components"},
  {"isfinite", pycomplexint16_isfinite, METH_NOARGS,
   "True if the complex_int16 has all finite components"},

  // Binary bool returners
  {"equal", pycomplexint16_equal, METH_O,
   "True if the complex_int16s are PRECISELY equal"},
  {"not_equal", pycomplexint16_not_equal, METH_O,
   "True if the complex_int16s are not PRECISELY equal"},
  {"less", pycomplexint16_less, METH_O,
   "Strict dictionary ordering"},
  {"greater", pycomplexint16_greater, METH_O,
   "Strict dictionary ordering"},
  {"less_equal", pycomplexint16_less_equal, METH_O,
   "Dictionary ordering"},
  {"greater_equal", pycomplexint16_greater_equal, METH_O,
   "Dictionary ordering"},

  // Unary ing returners
  {"real_part", pycomplexint16_real, METH_NOARGS,
   "real part of the complex value"},
  {"imag_part", pycomplexint16_imag, METH_NOARGS,
   "imaginary part of the complex value"},

  // Unary float returners
  {"absolute", pycomplexint16_absolute, METH_NOARGS,
   "Absolute value of complex_int16"},
  {"abs", pycomplexint16_absolute, METH_NOARGS,
   "Absolute value (Euclidean norm) of complex_int16"},
  {"norm", pycomplexint16_norm, METH_NOARGS,
   "Cayley norm (square of the absolute value) of complex_int16"},
  {"angle", pycomplexint16_angle, METH_NOARGS,
   "Angle through which rotor rotates"},

  // Unary complex_int16 returners
  // {"negative", pycomplexint16_negative, METH_NOARGS,
  //  "Return the negated complex_int16"},
  // {"positive", pycomplexint16_positive, METH_NOARGS,
  //  "Return the complex_int16 itself"},
  {"conjugate", pycomplexint16_conjugate, METH_NOARGS,
   "Return the complex conjugate of the complex_int16"},
  {"conj", pycomplexint16_conjugate, METH_NOARGS,
   "Return the complex conjugate of the complex_int16"},
  
  // complex_int16-complex_int16 or complex_int16-scalar binary complex_int16 returners
  // {"multiply", pycomplexint16_multiply, METH_O,
  //  "Standard (geometric) complex_int16 product"},
  // {"divide", pycomplexint16_divide, METH_O,
  //  "Standard (geometric) complex_int16 division"},
  // {"power", pycomplexint16_power, METH_O,
  //  "q.power(p) = (q.log() * p).exp()"},

  {"__reduce__", (PyCFunction)pycomplexint16__reduce, METH_NOARGS,
   "Return state information for pickling."},
  {"__getstate__", (PyCFunction)pycomplexint16_getstate, METH_VARARGS,
   "Return state information for pickling."},
  {"__setstate__", (PyCFunction)pycomplexint16_setstate, METH_VARARGS,
   "Reconstruct state information from pickle."},

  {NULL, NULL, 0, NULL}
};

static PyObject* pycomplexint16_num_negative(PyObject* a) { return pycomplexint16_negative(a,NULL); }
static PyObject* pycomplexint16_num_positive(PyObject* a) { return pycomplexint16_positive(a,NULL); }
static PyObject* pycomplexint16_num_absolute(PyObject* a) { return pycomplexint16_absolute(a,NULL); }
static int pycomplexint16_num_nonzero(PyObject* a) {
  complex_int16 q = ((PyComplexInt16*)a)->obval;
  return complex_int16_nonzero(q);
}
#define CANNOT_CONVERT_CI16(target)                                          \
  static PyObject* pycomplexint16_convert_##target(PyObject* a) {         \
    PyErr_SetString(PyExc_TypeError, "Cannot convert complex_int16 to " #target); \
    return NULL;                                                        \
  }
CANNOT_CONVERT_CI16(int)
CANNOT_CONVERT_CI16(float)
#if PY_MAJOR_VERSION < 3
CANNOT_CONVERT_CI16(long)
CANNOT_CONVERT_CI16(oct)
CANNOT_CONVERT_CI16(hex)
#endif

static PyNumberMethods pycomplexint16_as_number = {
  pycomplexint16_add,             // nb_add
  pycomplexint16_subtract,        // nb_subtract
  pycomplexint16_multiply,        // nb_multiply
  #if PY_MAJOR_VERSION < 3
  pycomplexint16_divide,          // nb_divide
  #endif
  0,                              // nb_remainder
  0,                              // nb_divmod
  0,                              // nb_power
  pycomplexint16_num_negative,    // nb_negative
  pycomplexint16_num_positive,    // nb_positive
  pycomplexint16_num_absolute,    // nb_absolute
  pycomplexint16_num_nonzero,     // nb_nonzero
  0,                              // nb_invert
  0,                              // nb_lshift
  0,                              // nb_rshift
  0,                              // nb_and
  0,                              // nb_xor
  0,                              // nb_or
  #if PY_MAJOR_VERSION < 3
  0,                              // nb_coerce
  #endif
  pycomplexint16_convert_int,       // nb_int
  #if PY_MAJOR_VERSION >= 3
  0,                              // nb_reserved
  #else
  pycomplexint16_convert_long,      // nb_long
  #endif
  pycomplexint16_convert_float,     // nb_float
  #if PY_MAJOR_VERSION < 3
  pycomplexint16_convert_oct,       // nb_oct
  pycomplexint16_convert_hex,       // nb_hex
  #endif
  pycomplexint16_inplace_add,       // nb_inplace_add
  pycomplexint16_inplace_subtract,  // nb_inplace_subtract
  pycomplexint16_inplace_multiply,  // nb_inplace_multiply
  #if PY_MAJOR_VERSION < 3
  pycomplexint16_inplace_divide,    // nb_inplace_divide
  #endif
  0,                              // nb_inplace_remainder
  0,                              // nb_inplace_power
  0,                              // nb_inplace_lshift
  0,                              // nb_inplace_rshift
  0,                              // nb_inplace_and
  0,                              // nb_inplace_xor
  0,                              // nb_inplace_or
  pycomplexint16_divide,          // nb_floor_divide
  pycomplexint16_divide,          // nb_true_divide
  pycomplexint16_inplace_divide,  // nb_inplace_floor_divide
  pycomplexint16_inplace_divide,  // nb_inplace_true_divide
  0,                              // nb_index
  #if PY_MAJOR_VERSION >= 3
  #if PY_MINOR_VERSION >= 5
  0,                              // nb_matrix_multiply
  0,                              // nb_inplace_matrix_multiply
  #endif
  #endif
};

// This is an array of members (member data) that will be available to
// use on the complex_int16 objects in python.  This is packaged up here,
// and will be used in the `tp_members` field when definining the
// PyComplexInt16_Type below.
PyMemberDef pycomplexint16_members[] = {
  {"real", T_BYTE, offsetof(PyComplexInt16, obval.real), 0,
   "The packed real component of the complex_int16"},
  {"imag", T_BYTE, offsetof(PyComplexInt16, obval.imag), 0,
   "The packed imaginary component of the complex_int16"},
  {NULL, 0, 0, 0, NULL}
};

static PyObject* pycomplexint16_richcompare(PyObject* a, PyObject* b, int op) {
    complex_int16 x = {0, 0};
    complex_int16 y = {0, 0};
    int result = 0;
    PyComplexInt16_AsComplexInt16(x,a);
    PyComplexInt16_AsComplexInt16(y,b);
    #define COMPARISONOP(py,op) case py: result = complex_int16_##op(x,y); break;
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

static long pycomplexint16_hash(PyObject *o) {
    complex_int16 q = ((PyComplexInt16 *)o)->obval;
    long value = 0x4567169;
    value = (10000004 * value) ^ _Py_HashDouble(q.real);
    value = (10000004 * value) ^ _Py_HashDouble(q.imag);
    if (value == -1) {
        value = -2;
    }
    return value;
}

static PyObject* pycomplexint16_repr(PyObject *o) {
    char str[1216];
    complex_int16 c = ((PyComplexInt16 *)o)->obval;
    sprintf(str, "complex_int16(%i, %i)", c.real, c.imag);
    return PyUString_FromString(str);
}

static PyObject* pycomplexint16_str(PyObject *o) {
    char str[64];
    complex_int16 c = ((PyComplexInt16 *)o)->obval;
    sprintf(str, "%i%+ij", c.real, c.imag);
    return PyUString_FromString(str);
}

// This establishes the complex_int16 as a python object (not yet a numpy
// scalar type).  The name may be a little counterintuitive; the idea
// is that this will be a type that can be used as an array dtype.
// Note that many of the slots below will be filled later, after the
// corresponding functions are defined.
static PyTypeObject PyComplexInt16_Type = {
#if PY_MAJOR_VERSION >= 3
  PyVarObject_HEAD_INIT(NULL, 0)
#else
  PyObject_HEAD_INIT(NULL)
  0,                                          // ob_size
#endif
  "complex_int16.complex_int16",                // tp_name
  sizeof(PyComplexInt16),                      // tp_basicsize
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
  pycomplexint16_repr,                         // tp_repr
  &pycomplexint16_as_number,                   // tp_as_number
  0,                                          // tp_as_sequence
  0,                                          // tp_as_mapping
  pycomplexint16_hash,                         // tp_hash
  0,                                          // tp_call
  pycomplexint16_str,                          // tp_str
  0,                                          // tp_getattro
  0,                                          // tp_setattro
  0,                                          // tp_as_buffer
#if PY_MAJOR_VERSION >= 3
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   // tp_flags
#else
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES, // tp_flags
#endif
  "Complex integer complex_int16 numbers",     // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  pycomplexint16_richcompare,                  // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  pycomplexint16_methods,                      // tp_methods
  pycomplexint16_members,                      // tp_members
  0,                                          // tp_getset
  0,                                          // tp_base; will be reset to &PyGenericArrType_Type after numpy import
  0,                                          // tp_dict
  0,                                          // tp_descr_get
  0,                                          // tp_descr_set
  0,                                          // tp_dictoffset
  pycomplexint16_init,                         // tp_init
  0,                                          // tp_alloc
  pycomplexint16_new,                          // tp_new
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
static PyArray_ArrFuncs _PyComplexInt16_ArrFuncs;
PyArray_Descr *complex_int16_descr;

static npy_bool CI16_nonzero (char *ip, PyArrayObject *ap) {
  complex_int16 c;
  complex_int16 zero = {0, 0};
  if (ap == NULL || PyArray_ISBEHAVED_RO(ap)) {
    c = *(complex_int16 *)ip;
  }
  else {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT8);
    descr->f->copyswap(&c.real, ip, !PyArray_ISNOTSWAPPED(ap), NULL);
    descr->f->copyswap(&c.imag, ip+1, !PyArray_ISNOTSWAPPED(ap), NULL);
    Py_DECREF(descr);
  }
  return (npy_bool) !complex_int16_equal(c, zero);
}

static void CI16_copyswap(complex_int16 *dst, complex_int16 *src, int swap, void *NPY_UNUSED(arr)) {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT8);
    descr->f->copyswapn(dst, sizeof(signed char), src, sizeof(signed char), 1, swap, NULL);
    Py_DECREF(descr);
}

static void CI16_copyswapn(complex_int16 *dst, npy_intp dstride,
                               complex_int16 *src, npy_intp sstride,
                               npy_intp n, int swap, void *NPY_UNUSED(arr)) {
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(NPY_INT8);
    descr->f->copyswapn(&dst->real, dstride, &src->real, sstride, n, swap, NULL);
    descr->f->copyswapn(&dst->imag, dstride, &src->imag, sstride, n, swap, NULL);
    Py_DECREF(descr);    
}

static int CI16_setitem(PyObject* item, complex_int16* c, void* NPY_UNUSED(ap)) {
    PyObject *element;
    if( PyComplexInt16_Check(item) ) {
        memcpy(c, &(((PyComplexInt16 *)item)->obval), sizeof(complex_int16));
    } else if( PyComplexInt8_Check(item) ) {
        complex_int8 i;
        PyComplexInt8_AsComplexInt8(i, item);
        const signed char* sc = fourBitLUT[i.real_imag];
        c->real = sc[0];
        c->imag = sc[1];
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
        PyErr_SetString(PyExc_TypeError, "Unknown input to CI16_setitem");
        return -1;
    }
    return 0;
}

// When a numpy array of dtype=complex_int16 is indexed, this function is
// called, returning a new complex_int16 object with a copy of the
// data... sometimes...
static PyObject* CI16_getitem(void* data, void* NPY_UNUSED(arr)) {
    complex_int16 q;
    memcpy(&q, data, sizeof(complex_int16));
    return PyComplexInt16_FromComplexInt16(q);
}

static int CI16_compare(complex_int16 *pa, complex_int16 *pb, PyArrayObject *NPY_UNUSED(ap)) {
    complex_int16 a = *pa, b = *pb;
    npy_bool anan, bnan;
    int ret;
    
    anan = complex_int16_isnan(a);
    bnan = complex_int16_isnan(b);
    
    if( anan ) {
        ret = bnan ? 0 : -1;
    } else if( bnan ) {
        ret = 1;
    } else if( complex_int16_less(a, b) ) {
        ret = -1;
    } else if( complex_int16_less(b, a) ) {
        ret = 1;
    } else {
        ret = 0;
    }
    
    return ret;
}

static int CI16_argmax(complex_int16 *ip, npy_intp n, npy_intp *max_ind, PyArrayObject *NPY_UNUSED(aip)) {
    npy_intp i;
    complex_int16 mp = *ip;
    
    *max_ind = 0;
    
    if( complex_int16_isnan(mp) ) {
        /* nan encountered; it's maximal */
        return 0;
    }

    for(i = 1; i < n; i++) {
        ip++;
        /*
         * Propagate nans, similarly as max() and min()
         */
        if( !(complex_int16_less_equal(*ip, mp)) ) {  /* negated, for correct nan handling */
            mp = *ip;
            *max_ind = i;
            if (complex_int16_isnan(mp)) {
                /* nan encountered, it's maximal */
                break;
            }
        }
    }
    return 0;
}

static void CI16_fillwithscalar(complex_int16 *buffer, npy_intp length, complex_int16 *value, void *NPY_UNUSED(ignored)) {
    npy_intp i;
    complex_int16 val = *value;

    for (i = 0; i < length; ++i) {
        buffer[i] = val;
    }
}

// This is a macro (followed by applications of the macro) that cast
// the input types to standard complex_int16 with only a nonzero scalar
// part.
#define MAKE_T_TO_CI16(TYPE, type)                                      \
static void TYPE ## _to_complex_int16(type *ip, complex_int16 *op, npy_intp n, \
                                  PyArrayObject *NPY_UNUSED(aip),      \
                                  PyArrayObject *NPY_UNUSED(aop)) {    \
    while (n--) {                                                      \
        op->real = (signed char) (*ip++);                              \
        op->imag = (signed char) (*ip++);                              \
        *op++;                                                         \
    }                                                                  \
}

MAKE_T_TO_CI16(BOOL, npy_bool);
MAKE_T_TO_CI16(BYTE, npy_byte);

static void CI8_to_complex_int16(complex_int8* ip, signed char *op, npy_intp n, PyArrayObject *NPY_UNUSED(aip), PyArrayObject *NPY_UNUSED(aop)) { 
    const signed char *sc;
    while (n--) {
        sc = fourBitLUT[ip->real_imag];
        *(op++) = (signed char) sc[0];
        *(op++) = (signed char) sc[1];
        (*ip++);
    }
}

// This is a macro (followed by applications of the macro) that cast
// the input complex types from complex_int16.
#define MAKE_CI16_TO_CT(TYPE, type)                                    \
static void complex_int16_to_## TYPE(complex_int16* ip, type *op, npy_intp n, \
                                PyArrayObject *NPY_UNUSED(aip),       \
                                PyArrayObject *NPY_UNUSED(aop)) {     \
    while (n--) {                                                     \
        *(op++) = (type) ip->real;                                    \
        *(op++) = (type) ip->imag;                                    \
        (*ip++);                                                      \
    }                                                                 \
}

MAKE_CI16_TO_CT(CFLOAT, npy_float);
MAKE_CI16_TO_CT(CDOUBLE, npy_double);
MAKE_CI16_TO_CT(CLONGDOUBLE, npy_longdouble);

static void register_cast_function_ci16(int sourceType, int destType, PyArray_VectorUnaryFunc *castfunc) {
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
// complex_int16 functions, so that they can be applied quickly to a
// numpy array of complex_int16s.
#define UNARY_GEN_UFUNC_CI16(ufunc_name, func_name, ret_type)                \
    static void complex_int16_##ufunc_name##_ufunc(char** args,              \
                                                  npy_intp* dimensions,     \
                                                  npy_intp* steps,          \
                                                  void* NPY_UNUSED(data)) { \
    char *ip1 = args[0], *op1 = args[1];                                    \
    npy_intp is1 = steps[0], os1 = steps[1];                                \
    npy_intp n = dimensions[0];                                             \
    npy_intp i;                                                             \
    for(i = 0; i < n; i++, ip1 += is1, op1 += os1){                         \
      const complex_int16 in1 = *(complex_int16 *)ip1;                        \
      *((ret_type *)op1) = complex_int16_##func_name(in1);};}
#define UNARY_UFUNC_CI16(name, ret_type) \
  UNARY_GEN_UFUNC_CI16(name, name, ret_type)
// And these all do the work mentioned above, using the macro
UNARY_UFUNC_CI16(isnan, npy_bool)
UNARY_UFUNC_CI16(isinf, npy_bool)
UNARY_UFUNC_CI16(isfinite, npy_bool)
UNARY_UFUNC_CI16(real, npy_long)
UNARY_UFUNC_CI16(imag, npy_long)
UNARY_UFUNC_CI16(norm, npy_double)
UNARY_UFUNC_CI16(absolute, npy_double)
UNARY_UFUNC_CI16(angle, npy_double)
UNARY_UFUNC_CI16(negative, complex_int16)
UNARY_UFUNC_CI16(conjugate, complex_int16)
static void complex_int16_positive_ufunc(char** args, npy_intp* dimensions, npy_intp* steps, void* NPY_UNUSED(data)) {
    char *ip1 = args[0], *op1 = args[1];
    npy_intp is1 = steps[0], os1 = steps[1];
    npy_intp n = dimensions[0];
    npy_intp i;
    for(i = 0; i < n; i++, ip1 += is1, op1 += os1) {
        const complex_int16 in1 = *(complex_int16 *)ip1;
        *((complex_int16 *)op1) = in1;
    }
}

// This is a macro that will be used to define the various basic binary
// complex_int16 functions, so that they can be applied quickly to a
// numpy array of complex_int16s.
#define BINARY_GEN_UFUNC_CI16(ufunc_name, func_name, arg_type1, arg_type2, ret_type) \
  static void complex_int16_##ufunc_name##_ufunc(char** args,                     \
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
      *((ret_type *)op1) = complex_int16_##func_name(in1, in2);                   \
    };                                                                          \
  };
// A couple special-case versions of the above
#define BINARY_UFUNC_CI16(name, ret_type)                    \
  BINARY_GEN_UFUNC_CI16(name, name, complex_int16, complex_int16, ret_type)
#define BINARY_SCALAR_UFUNC_CI16(name, ret_type)                             \
  BINARY_GEN_UFUNC_CI16(name##_scalar, name##_scalar, complex_int16, npy_long, ret_type) \
  BINARY_GEN_UFUNC_CI16(scalar_##name, scalar_##name, npy_long, complex_int16, ret_type)
// And these all do the work mentioned above, using the macros
BINARY_UFUNC_CI16(add, complex_int16)
BINARY_UFUNC_CI16(subtract, complex_int16)
BINARY_UFUNC_CI16(multiply, complex_int16)
BINARY_UFUNC_CI16(divide, complex_int16)
BINARY_GEN_UFUNC_CI16(true_divide, divide, complex_int16, complex_int16, complex_int16)
BINARY_GEN_UFUNC_CI16(floor_divide, divide, complex_int16, complex_int16, complex_int16)
BINARY_UFUNC_CI16(equal, npy_bool)
BINARY_UFUNC_CI16(not_equal, npy_bool)
BINARY_UFUNC_CI16(less, npy_bool)
BINARY_UFUNC_CI16(less_equal, npy_bool)
BINARY_SCALAR_UFUNC_CI16(add, complex_int16)
BINARY_SCALAR_UFUNC_CI16(subtract, complex_int16)
BINARY_SCALAR_UFUNC_CI16(multiply, complex_int16)
BINARY_SCALAR_UFUNC_CI16(divide, complex_int16)
BINARY_GEN_UFUNC_CI16(true_divide_scalar, divide_scalar, complex_int16, long, complex_int16)
BINARY_GEN_UFUNC_CI16(floor_divide_scalar, divide_scalar, complex_int16, long, complex_int16)
BINARY_GEN_UFUNC_CI16(scalar_true_divide, scalar_divide, long, complex_int16, complex_int16)
BINARY_GEN_UFUNC_CI16(scalar_floor_divide, scalar_divide, long, complex_int16, complex_int16)

static PyObject* complex_int16_arrtype_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    complex_int16 c;
    Py_complex cmplx;

    if( !PyArg_ParseTuple(args, "D", &cmplx) ) {
        return NULL;
    }
    
    // Ouch
    c.real = cmplx.real;
    c.imag = cmplx.imag;
    return PyArray_Scalar(&c, complex_int16_descr, NULL);
}

static PyObject* gentype_richcompare_ci16(PyObject *self, PyObject *other, int cmp_op) {
    PyObject *arr, *ret;

    arr = PyArray_FromScalar(self, NULL);
    if (arr == NULL) {
        return NULL;
    }
    ret = Py_TYPE(arr)->tp_richcompare(arr, other, cmp_op);
    Py_DECREF(arr);
    return ret;
}

static long complex_int16_arrtype_hash(PyObject *o) {
    complex_int16 c = ((PyComplexInt16 *)o)->obval;
    long value = 0x4567169;
    value = (10000004 * value) ^ _Py_HashDouble(c.real);
    value = (10000004 * value) ^ _Py_HashDouble(c.imag);
    if( value == -1 ) {
        value = -2;
    }
    return value;
}

static PyObject* complex_int16_arrtype_repr(PyObject *o) {
    char str[64];
    complex_int16 c = ((PyComplexInt16 *)o)->obval;
    sprintf(str, "complex_int16(%i, %i)", c.real, c.imag);
    return PyUString_FromString(str);
}

static PyObject* complex_int16_arrtype_str(PyObject *o) {
    char str[64];
    complex_int16 c = ((PyComplexInt16 *)o)->obval;
    sprintf(str, "%i%+ij", c.real, c.imag);
    return PyUString_FromString(str);
}

int complex_int16_elsize = sizeof(complex_int16);

typedef struct {
    char c;
    complex_int16 q;
} align_test_ci16;
int complex_int16_alignment = offsetof(align_test_ci16, q);

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//                                                             //
//  Everything above was preparation for the following set up  //
//                                                             //
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

int create_complex_int16(PyObject* m, PyObject* numpy_dict) {
    int complexi16Num;
    PyObject *tmp_ufunc;
    int arg_types[3];
    
    /* Register the complex_int16 array base type. */
    PyComplexInt16_Type.tp_base = &PyGenericArrType_Type;
    if (PyType_Ready(&PyComplexInt16_Type) < 0) {
        return -2;
    }
    
    // The array functions, to be used below.  This InitArrFuncs
    // function is a convenient way to set all the fields to zero
    // initially, so we don't get undefined behavior.
    PyArray_InitArrFuncs(&_PyComplexInt16_ArrFuncs);
    _PyComplexInt16_ArrFuncs.nonzero = (PyArray_NonzeroFunc*)CI16_nonzero;
    _PyComplexInt16_ArrFuncs.copyswap = (PyArray_CopySwapFunc*)CI16_copyswap;
    _PyComplexInt16_ArrFuncs.copyswapn = (PyArray_CopySwapNFunc*)CI16_copyswapn;
    _PyComplexInt16_ArrFuncs.setitem = (PyArray_SetItemFunc*)CI16_setitem;
    _PyComplexInt16_ArrFuncs.getitem = (PyArray_GetItemFunc*)CI16_getitem;
    _PyComplexInt16_ArrFuncs.compare = (PyArray_CompareFunc*)CI16_compare;
    _PyComplexInt16_ArrFuncs.argmax = (PyArray_ArgFunc*)CI16_argmax;
    _PyComplexInt16_ArrFuncs.fillwithscalar = (PyArray_FillWithScalarFunc*)CI16_fillwithscalar;
    
    // The complex_int16 array descr
    complex_int16_descr = PyObject_New(PyArray_Descr, &PyArrayDescr_Type);
    complex_int16_descr->typeobj = &PyComplexInt16_Type;
    complex_int16_descr->kind = 'V';
    complex_int16_descr->type = 'i';
    complex_int16_descr->byteorder = '=';
    complex_int16_descr->flags = NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM;
    complex_int16_descr->type_num = 0; // assigned at registration
    complex_int16_descr->elsize = complex_int16_elsize;
    complex_int16_descr->alignment = complex_int16_alignment;
    complex_int16_descr->subarray = NULL;
    complex_int16_descr->fields = NULL;
    complex_int16_descr->names = NULL;
    complex_int16_descr->f = &_PyComplexInt16_ArrFuncs;
    complex_int16_descr->metadata = NULL;
    complex_int16_descr->c_metadata = NULL;
    
    Py_INCREF(&PyComplexInt16_Type);
    complexi16Num = PyArray_RegisterDataType(complex_int16_descr);
    
    if( complexi16Num < 0 || complexi16Num != NPY_COMPLEX_INT16 ) {
        return -1;
    }
    
    register_cast_function_ci16(NPY_BOOL, complexi16Num, (PyArray_VectorUnaryFunc*)BOOL_to_complex_int16);
    register_cast_function_ci16(NPY_BYTE, complexi16Num, (PyArray_VectorUnaryFunc*)BYTE_to_complex_int16);
    
    register_cast_function_ci16(NPY_COMPLEX_INT8, complexi16Num, (PyArray_VectorUnaryFunc*)CI8_to_complex_int16);
    
    register_cast_function_ci16(complexi16Num, NPY_CFLOAT, (PyArray_VectorUnaryFunc*)complex_int16_to_CFLOAT);
    register_cast_function_ci16(complexi16Num, NPY_CDOUBLE, (PyArray_VectorUnaryFunc*)complex_int16_to_CDOUBLE);
    register_cast_function_ci16(complexi16Num, NPY_CLONGDOUBLE, (PyArray_VectorUnaryFunc*)complex_int16_to_CLONGDOUBLE);
    
    
    // These macros will be used below
#define REGISTER_UFUNC_CI16(name)\
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name),\
            complex_int16_descr->type_num, complex_int16_##name##_ufunc, arg_types, NULL)
    
#define REGISTER_SCALAR_UFUNC_CI16(name)\
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name),\
            complex_int16_descr->type_num, complex_int16_##name##_scalar_ufunc, arg_types, NULL)
    
#define REGISTER_UFUNC_SCALAR_CI16(name)                                   \
    PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(numpy_dict, #name), \
                                complex_int16_descr->type_num, complex_int16_##name##_scalar_ufunc, arg_types, NULL)
    
#define REGISTER_NEW_UFUNC_GENERAL_CI16(pyname, cname, nargin, nargout, doc) \
    tmp_ufunc = PyUFunc_FromFuncAndData(NULL, NULL, NULL, 0, nargin, nargout, \
                                        PyUFunc_None, #pyname, doc, 0); \
    PyUFunc_RegisterLoopForType((PyUFuncObject *)tmp_ufunc,             \
                                complex_int16_descr->type_num, complex_int16_##cname##_ufunc, arg_types, NULL); \
    PyDict_SetItemString(numpy_dict, #pyname, tmp_ufunc);               \
    Py_DECREF(tmp_ufunc)

  #define REGISTER_NEW_UFUNC_CI16(name, nargin, nargout, doc)                \
    REGISTER_NEW_UFUNC_GENERAL_CI16(name, name, nargin, nargout, doc)
    
    /* complex_int16 -> bool */
    arg_types[0] = complex_int16_descr->type_num;
    arg_types[1] = NPY_BOOL;
    
    REGISTER_UFUNC_CI16(isnan);
    REGISTER_UFUNC_CI16(isinf);
    REGISTER_UFUNC_CI16(isfinite);
    
    /* complex_int16 -> long */
    arg_types[1] = NPY_LONG;
    
    REGISTER_NEW_UFUNC_GENERAL_CI16(real_part, real, 1, 1, \
                                    "Return the real part of the complex number\n");
    REGISTER_NEW_UFUNC_GENERAL_CI16(imag_part, imag, 1, 1, \
                                    "Return the imaginary part of the complex number\n");
    
    /* complex_int16 -> double */
    arg_types[1] = NPY_DOUBLE;
    
    REGISTER_NEW_UFUNC_CI16(norm, 1, 1, \
                            "Return Cayley norm (square of the absolute value) of each complex_int16\n");
    REGISTER_UFUNC_CI16(absolute);
    REGISTER_NEW_UFUNC_GENERAL_CI16(angle, angle, 1, 1, \
                                    "Return the angle of the complex argument\n");
    
    /* complex_int16 -> complex_int16 */
    arg_types[1] = complex_int16_descr->type_num;
    
    REGISTER_UFUNC_CI16(negative);
    REGISTER_UFUNC_CI16(conjugate);

    /* complex_int16, complex_int16 -> bool */
    arg_types[2] = NPY_BOOL;

    REGISTER_UFUNC_CI16(equal);
    REGISTER_UFUNC_CI16(not_equal);
    REGISTER_UFUNC_CI16(less);
    REGISTER_UFUNC_CI16(less_equal);
    
    /* complex_int16, complex_int16 -> complex_int16 */
    arg_types[2] = complex_int16_descr->type_num;
    
    REGISTER_UFUNC_CI16(add);
    REGISTER_UFUNC_CI16(subtract);
    REGISTER_UFUNC_CI16(multiply);
    REGISTER_UFUNC_CI16(divide);
    REGISTER_UFUNC_CI16(true_divide);
    REGISTER_UFUNC_CI16(floor_divide);
    
    /* long, complex_int16 -> complex_int16 */
    arg_types[0] = NPY_LONG;
    arg_types[1] = complex_int16_descr->type_num;
    arg_types[2] = complex_int16_descr->type_num;
    
    REGISTER_SCALAR_UFUNC_CI16(add);
    REGISTER_SCALAR_UFUNC_CI16(subtract);
    REGISTER_SCALAR_UFUNC_CI16(multiply);
    REGISTER_SCALAR_UFUNC_CI16(divide);
    REGISTER_SCALAR_UFUNC_CI16(true_divide);
    REGISTER_SCALAR_UFUNC_CI16(floor_divide);
    
    /* complex_int16, long -> complex_int16 */
    arg_types[0] = complex_int16_descr->type_num;
    arg_types[1] = NPY_LONG;
    arg_types[2] = complex_int16_descr->type_num;
    
    REGISTER_UFUNC_SCALAR_CI16(add);
    REGISTER_UFUNC_SCALAR_CI16(subtract);
    REGISTER_UFUNC_SCALAR_CI16(multiply);
    REGISTER_UFUNC_SCALAR_CI16(divide);
    REGISTER_UFUNC_SCALAR_CI16(true_divide);
    REGISTER_UFUNC_SCALAR_CI16(floor_divide);
    
    PyModule_AddObject(m, "complex_int16", (PyObject *)&PyComplexInt16_Type);
    
    return complexi16Num;
}
