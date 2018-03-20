/* -*- C -*-  (not really, but good for syntax highlighting) */

%module c_m8r

%{
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <setjmp.h>

#include <rsf.h>

#define SWIG_FILE_WITH_INIT

jmp_buf jump_buffer;
extern jmp_buf *python_except;
%}

/* Get the Numeric typemaps */
%include numpy.i

%init %{
    import_array();
%}

%include typemaps.i
%apply int *OUTPUT { int *par };
%apply float *OUTPUT { float *par };

%typemap(in) (int *pars, size_t n) {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int i;
    $2 = PyList_Size($input);
    $1 = sf_intalloc($2);
    for (i = 0; i < $2; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyInt_Check(o))
        $1[i] = PyInt_AsLong(PyList_GetItem($input,i));
      else {
        PyErr_SetString(PyExc_TypeError,"list must contain ints");
        free($1);
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

%typemap(freearg) (int *pars, size_t n) {
  free($1);
}

%typemap(out) struct farray* {
   int i;   
   $result = PyList_New($1->n);
   for (i=0; i < $1->n; i++) {
       PyList_SetItem($result, i, PyFloat_FromDouble($1->arr[i]));
   }
   free ($1);
}

%typemap(out) struct iarray* {
   int i;   
   $result = PyList_New($1->n);
   for (i=0; i < $1->n; i++) {
       PyList_SetItem($result, i, PyInt_FromLong($1->arr[i]));
   }
   free ($1);
}

%typemap(out) struct barray* {
   int i;   
   $result = PyList_New($1->n);
   for (i=0; i < $1->n; i++) {
       PyList_SetItem($result, i, PyInt_FromLong($1->arr[i]? 1: 0));
   }
   free ($1);
}

%inline %{
    struct farray {
	       float *arr;
	       int n;
        };	

    struct iarray {
	int *arr;
	int n;
    };

    struct barray {
	       bool *arr;
	       int n;
        };

    struct iarray* histints(sf_file file, const char* key,size_t n) {
	       struct iarray* a;
	       a = sf_alloc(1,sizeof(struct iarray));
	       a->n = n;
	       a->arr = sf_intalloc(n);
	       if (sf_histints (file,key,a->arr,n)) {
		  return a;
	       } else {
                  free (a->arr);
		  free (a);
                  return NULL;
               }
         }
%}

%typemap(in) char** {
if (PyList_Check($input)) {
   int size = PyList_Size($input);
   int i = 0;
   $1 = sf_alloc(size+1,sizeof(char*));
   for (i=0; i < size; i++) {
       PyObject *o = PyList_GetItem($input,i);
       if (PyString_Check(o)) {
	  $1[i] = PyString_AsString(PyList_GetItem($input,i));
       } else {
	  PyErr_SetString(PyExc_TypeError,"list must contain strings");
	  free($1);
	  return NULL;
       }
   }
   $1[i] = 0;
} else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
}
}

%typemap(freearg) char** {
free ((char*) $1);
}

%typemap(in,numinputs=0) bool* par (bool temp) {
    $1 = &temp;
}

%typemap(argout) bool *par {
    PyObject *o, *o2, *o3;
    if (*$1) {
       o = PyInt_FromLong(1);
    } else {
       o = PyInt_FromLong(0);
    }
    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        if (!PyTuple_Check($result)) {
            PyObject *o2 = $result;
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}

%inline %{
	struct farray* getfloats(const char* key,size_t n) {
	       struct farray* a;

	       a = sf_alloc(1,sizeof(struct farray));
	       a->n = n;
	       a->arr = sf_floatalloc(n);
	       if (sf_getfloats (key,a->arr,n)) {
		   return a;
	       } else {
		   free (a->arr);
		   free (a);
		   return NULL;
               }
         }

	 struct iarray* getints(const char* key,size_t n) {
	       struct iarray* a;
	       a = sf_alloc(1,sizeof(struct iarray));
	       a->n = n;
	       a->arr = sf_intalloc(n);
	       if (sf_getints (key,a->arr,n)) {
		  return a;
	       } else {
                  free (a->arr);
		  free (a);
                  return NULL;
               }
         }

	 struct barray* getbools(const char* key,size_t n) {
	       struct barray* a;
	       a = sf_alloc(1,sizeof(struct barray));
	       a->n = n;
	       a->arr = sf_boolalloc(n);
	       if (sf_getbools (key,a->arr,n)) {
		  return a;
	       } else {
                  free (a->arr);
		  free (a);
                  return NULL;
               }
         }
%}

/* Apply the Numeric typemaps for 1D input arrays */
%numpy_typemaps(float, NPY_FLOAT   , size_t)
%numpy_typemaps(int, NPY_INT   , size_t)
%numpy_typemaps(sf_complex, NPY_CFLOAT   , size_t)

%apply (float*  IN_ARRAY1, size_t DIM1) {(float*  arr, size_t size)};
%apply (int*  IN_ARRAY1, size_t DIM1) {(int*  arr, size_t size)};
%apply (sf_complex*  IN_ARRAY1, size_t DIM1) {(sf_complex*  arr, size_t size)};

#define SF_MAX_DIM 9

typedef struct sf_File *sf_file;

typedef enum {SF_UCHAR, SF_CHAR, SF_INT, SF_FLOAT, SF_COMPLEX} sf_datatype;
typedef enum {SF_ASCII, SF_XDR, SF_NATIVE} sf_dataform;

%exception {
    python_except = &jump_buffer;
    if (!setjmp(*python_except)) {
	$action
    } else {
	PyErr_SetString(PyExc_RuntimeError,"sf_error");
	return NULL;
    }
}

sf_file sf_input (/*@null@*/ char* tag);
sf_file sf_output (/*@null@*/ char* tag);
sf_datatype sf_gettype (sf_file file);
sf_dataform sf_getform (sf_file file);
void sf_settype (sf_file file, sf_datatype type);
void sf_setformat (sf_file file, const char* format);
void sf_setaformat (const char* format, int line);
void sf_fileclose (sf_file file);
bool sf_histint (sf_file file, const char* key,/*@out@*/ int* par);
bool sf_histints (sf_file file, const char* key,/*@out@*/ int* par, size_t n);
bool sf_histfloat (sf_file file, const char* key,/*@out@*/ float* par);
char* sf_histstring (sf_file file, const char* key);
void sf_fileflush (sf_file file, sf_file src);
void sf_putint (sf_file file, char* key, int par);
void sf_putints (sf_file file, char* key, int* pars, size_t n);
void sf_putfloat (sf_file file, char* key,float par);
void sf_putstring (sf_file file, char* key,const char* par);
void sf_putline (sf_file file, const char* line);
long sf_bytes (sf_file file);

void sf_floatwrite (float* arr, size_t size, sf_file file);
void sf_floatread (float* arr, size_t size, sf_file file);

void sf_intwrite (int* arr, size_t size, sf_file file);
void sf_intread (int* arr, size_t size, sf_file file);

int sf_try_charread2(char* arr, size_t size, sf_file file);

void sf_complexwrite (sf_complex* arr, size_t size, sf_file file);
void sf_complexread (sf_complex* arr, size_t size, sf_file file);

long sf_tell (sf_file file);
void sf_seek (sf_file file, long offset, int whence);
void sf_unpipe (sf_file file, size_t size);

void sf_init(int argc,char **argv);
void sf_parclose (void);
char* sf_getprog (void);
char* sf_getuser (void);
char* sf_gethost (void);

bool sf_getint (const char* key,/*@out@*/ int* par);
bool sf_getints (const char* key,/*@out@*/ int* par,size_t n);

bool sf_getfloat (const char* key,/*@out@*/ float* par);
bool sf_getfloats (const char* key,float* par,size_t n);

char* sf_getstring (const char* key);
bool sf_getstrings (const char* key,/*@out@*/ char** par,size_t n);

bool sf_getbool (const char* key,/*@out@*/ bool* par);
bool sf_getbools (const char* key,/*@out@*/ bool* par,size_t n);

int sf_filedims (sf_file file, /*@out@*/ int *dims);
int sf_filesize (sf_file file);
int sf_leftsize (sf_file file, int dim);
void sf_cp(sf_file in, sf_file out);
void sf_rm(const char* filename, bool force, bool verb, bool inquire);

%inline %{
    PyObject* rsf_array(sf_file file)
    {
	PyObject *array;
	char *filename, *data;
	int id, fd, nd, n[SF_MAX_DIM], type;
	npy_intp *dims;
	struct stat statbuf;
	
	/* memory map */
	if (NULL == (filename = sf_histstring(file,"in"))) return NULL;
	if (0 > (fd = open(filename,O_RDWR))) return NULL;
	if (0 > (fstat(fd,&statbuf))) return NULL;
	data = mmap(0,statbuf.st_size,PROT_READ,MAP_FILE | MAP_SHARED,fd,0);
	
	/* dimensions and type */
	nd = sf_filedims (file, n);
	dims = sf_alloc(nd,sizeof(*dims));
	for (id=0; id < nd; id++) dims[id]=n[nd-id-1];

	switch (sf_gettype(file)) {
	    case SF_FLOAT:
		type = NPY_FLOAT;
		break;
	    case SF_INT:
		type = NPY_INT;
		break;
	    case SF_COMPLEX:
		type = NPY_CFLOAT;
		break;
	    case SF_CHAR:
		type = NPY_CHAR;
		break;
	    case SF_UCHAR:
		type = NPY_UBYTE;
		break;
	    default:
		type = NPY_FLOAT;
		break;
	}

	array = PyArray_SimpleNewFromData(nd,dims,type,data);
	return array;
    }
%}
