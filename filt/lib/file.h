#ifndef _sf_file_h
#define _sf_file_h

#include <stdio.h>

#include "c99.h"

#define SF_MAX_DIM 9

typedef struct sf_File *sf_file;

typedef enum {SF_CHAR, SF_INT, SF_FLOAT, SF_COMPLEX} sf_datatype;
typedef enum {SF_ASCII, SF_XDR, SF_NATIVE} sf_dataform;

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
void sf_putints (sf_file file, char* key, int* par, size_t n);
void sf_putfloat (sf_file file, char* key,float par);
void sf_putstring (sf_file file, char* key,const char* par);
void sf_putline (sf_file file, const char* line);
long sf_bytes (sf_file file);
void sf_write (void* arr, size_t esize, size_t size, sf_file file);
void sf_read (/*@out@*/ void* arr, size_t esize, size_t size, sf_file file);
long sf_tell (sf_file file);
void sf_seek (sf_file file, long offset, int whence);
void sf_unpipe (sf_file file, size_t size);

#endif
