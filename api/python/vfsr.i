/* -*- C -*-  (not really, but good for syntax highlighting) */
#ifdef SWIGPYTHON

%module c_vfsr
%{

#include <stdio.h>
#include <numpy/arrayobject.h>
#include <rsf.h>

#define SWIG_FILE_WITH_INIT
%}

/* Get the Numeric typemaps */
%include numpy.i

%init %{
    import_array();
%}

%include typemaps.i

// Grab a Python function object as a Python object.
%typemap(in) PyObject *pyfunc {
  if (!PyCallable_Check ($input)) {
      PyErr_SetString (PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $1 = $input;
}

%include ../../filt/lib/vfsr_defs.h

%inline %{

typedef struct {
    PyObject *cost_func;
    int       number_parameters;
} VFSR_USER_DATA;

static double vfsr_cost_function_cb (double *cost_parameters,
                                     double *parameter_lower_bound,
                                     double *parameter_upper_bound,
                                     int *cost_flag,
                                     VFSR_DEFINES *USER_OPTIONS,
                                     void *user_data)
{
    PyObject *func, *arglist;
    PyObject *params, *lower_bound, *upper_bound, *flag;
    PyObject *result;
    int dims[1] = { ((VFSR_USER_DATA*)user_data)->number_parameters };
    int flagdims[1] = { 1 };
    double dres = 0;

    func = (PyObject *)((VFSR_USER_DATA*)user_data)->cost_func;
    params = (PyObject *)PyArray_FromDimsAndData (1, dims, PyArray_DOUBLE,
                                                  (char *)cost_parameters);
    lower_bound = (PyObject *)PyArray_FromDimsAndData (1, dims, PyArray_DOUBLE,
                                                       (char *)parameter_lower_bound);
    upper_bound = (PyObject *)PyArray_FromDimsAndData (1, dims, PyArray_DOUBLE,
                                                       (char *)parameter_upper_bound);
    flag = (PyObject *)PyArray_FromDimsAndData (1, flagdims, PyArray_INT,
                                                (char *)cost_flag);

    arglist = Py_BuildValue ("OOOO", params, lower_bound, upper_bound, flag);
    result = PyEval_CallObject (func, arglist);
    Py_DECREF (arglist);

    if (result) {
        dres = PyFloat_AsDouble (result);
    } else {
        PyErr_Print ();
        *cost_flag = FALSE;
    }

    Py_XDECREF (result);
    Py_XDECREF (params);
    Py_XDECREF (lower_bound);
    Py_XDECREF (upper_bound);
    Py_XDECREF (flag);

    return dres;
}
%}

%inline %{

double run (PyObject *cost_func,
            PyObject *parameter_type,
            PyObject *parameter_initial_final,
            PyObject *parameter_minimum, PyObject *parameter_maximum,
            PyObject *tangents, PyObject *curvature,
            PyObject *exit_status,
            VFSR_DEFINES *OPTIONS) {
    VFSR_USER_DATA user_data;
    user_data.cost_func = cost_func;
    PyArrayObject *param_type, *param_init, *param_min, *param_max,
                  *tans, *curv, *status;

    if (!PyArray_Check (parameter_type)) {
        PyErr_SetString (PyExc_TypeError, "parameter_type is not an array");
        return 0;
    }
    param_type = (PyArrayObject *)parameter_type;
    if (!PyArray_Check (parameter_initial_final)) {
        PyErr_SetString (PyExc_TypeError, "parameter_initial_final is not an array");
        return 0;
    }
    param_init = (PyArrayObject *)parameter_initial_final;
    if (!PyArray_Check (parameter_minimum)) {
        PyErr_SetString (PyExc_TypeError, "parameter_minimum is not an array");
        return 0;
    }
    param_min = (PyArrayObject *)parameter_minimum;
    if (!PyArray_Check (parameter_maximum)) {
        PyErr_SetString (PyExc_TypeError, "parameter_maximum is not an array");
        return 0;
    }
    param_max = (PyArrayObject *)parameter_maximum;
    if (!PyArray_Check (tangents)) {
        PyErr_SetString (PyExc_TypeError, "tangents is not an array");
        return 0;
    }
    tans = (PyArrayObject *)tangents;
    if (!PyArray_Check (curvature)) {
        PyErr_SetString (PyExc_TypeError, "curvature is not an array");
        return 0;
    }
    curv = (PyArrayObject *)curvature;
    if (!PyArray_Check (exit_status)) {
        PyErr_SetString (PyExc_TypeError, "exit_status is not an array");
        return 0;
    }
    status = (PyArrayObject *)exit_status;

    user_data.number_parameters = param_type->dimensions[0];
    Py_INCREF (parameter_type);
    Py_INCREF (parameter_initial_final);
    Py_INCREF (parameter_minimum);
    Py_INCREF (parameter_maximum);
    Py_INCREF (tangents);
    Py_INCREF (curvature);

    Py_INCREF (cost_func);

    double res= vfsr_std_rng (vfsr_cost_function_cb, param_type->dimensions[0],
                              (int*)param_type->data, (double*)param_init->data,
                              (double*)param_min->data, (double*)param_max->data,
                              (double*)tans->data, (double*)curv->data,
                              (int*)status->data, OPTIONS,
                              (void*)&user_data);

    Py_DECREF (parameter_type);
    Py_DECREF (parameter_initial_final);
    Py_DECREF (parameter_minimum);
    Py_DECREF (parameter_maximum);
    Py_DECREF (tangents);
    Py_DECREF (curvature);

    Py_DECREF (cost_func);

    return res;
}
%}

#endif /* SWIGPYTHON */
