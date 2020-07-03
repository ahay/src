/* -*- C -*-  (not really, but good for syntax highlighting) */

%module m8r

%{
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <setjmp.h>

#include <rsf.h>

#define SWIG_FILE_WITH_INIT

%}

%include typemaps.i
%apply int *OUTPUT { int *par };
%apply float *OUTPUT { float *par };

%typemap(freearg) (int *pars, size_t n) {
  free($1);
}
/*
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

%}
*/
/* this comes verbatim from the SWIG java tutorial */
%typemap(in) char ** (jint size) {
    int i = 0;
    size = (*jenv)->GetArrayLength(jenv, $input);
    $1 = (char **) malloc((size+1)*sizeof(char *));
    /* make a copy of each string */
    for (i = 0; i<size; i++) {
        jstring j_string = (jstring)(*jenv)->GetObjectArrayElement(jenv, $input, i);
        const char * c_string = (*jenv)->GetStringUTFChars(jenv, j_string, 0);
        $1[i] = malloc((strlen(c_string)+1)*sizeof(char));
        strcpy($1[i], c_string);
        (*jenv)->ReleaseStringUTFChars(jenv, j_string, c_string);
        (*jenv)->DeleteLocalRef(jenv, j_string);
    }
    $1[i] = 0;
}

/* This cleans up the memory we malloc'd before the function call */
%typemap(freearg) char ** {
    int i;
    for (i=0; i<size$argnum-1; i++)
      free($1[i]);
    free($1);
}

/* This allows a C function to return a char ** as a Java String array */
%typemap(out) char ** {
    int i;
    int len=0;
    jstring temp_string;
    const jclass clazz = (*jenv)->FindClass(jenv, "java/lang/String");

    while ($1[len]) len++;    
    jresult = (*jenv)->NewObjectArray(jenv, len, clazz, NULL);
    /* exception checking omitted */

    for (i=0; i<len; i++) {
      temp_string = (*jenv)->NewStringUTF(jenv, *result++);
      (*jenv)->SetObjectArrayElement(jenv, jresult, i, temp_string);
      (*jenv)->DeleteLocalRef(jenv, temp_string);
    }
}

/* These 3 typemaps tell SWIG what JNI and Java types to use */
%typemap(jni) char ** "jobjectArray"
%typemap(jtype) char ** "String[]"
%typemap(jstype) char ** "String[]"

/* These 2 typemaps handle the conversion of the jtype to jstype typemap type
   and vice versa */
%typemap(javain) char ** "$javainput"
%typemap(javaout) char ** {
    return $jnicall;
  }

%typemap(in,numinputs=0) bool* par (bool temp) {
    $1 = &temp;
}

#define SF_MAX_DIM 9

typedef struct sf_File *sf_file;

typedef enum {SF_UCHAR, SF_CHAR, SF_INT, SF_FLOAT, SF_COMPLEX} sf_datatype;
typedef enum {SF_ASCII, SF_XDR, SF_NATIVE} sf_dataform;

sf_file sf_input (/*@null@*/ char* tag);
sf_file sf_output (/*@null@*/ char* tag);
sf_datatype sf_gettype (sf_file file);
sf_dataform sf_getform (sf_file file);
void sf_settype (sf_file file, sf_datatype type);
void sf_setformat (sf_file file, const char* format);
void sf_setaformat (const char* format, int line, int strip);
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

%include "arrays_java.i"
%apply float[] {float *};

void sf_floatwrite (float* arr, size_t size, sf_file file);
void sf_floatread (float* arr, size_t size, sf_file file);

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

bool sf_getbool (const char* key,/*@out@*/ bool* par);
bool sf_getbools (const char* key,/*@out@*/ bool* par,size_t n);

char* sf_getstring(const char* key);

int sf_filedims (sf_file file, /*@out@*/ int *dims);
int sf_filesize (sf_file file);
int sf_leftsize (sf_file file, int dim);
void sf_cp(sf_file in, sf_file out);
void sf_rm(const char* filename, bool force, bool verb, bool inquire);
