#ifndef _sf_alloc_h
#define _sf_alloc_h

#include <stdlib.h>

#include "c99.h"

/*@out@*/ void *sf_alloc (size_t n, size_t size);
void *sf_realloc (void* ptr, size_t n, size_t size);
/*@out@*/ char *sf_charalloc (size_t n);
/*@out@*/ int *sf_intalloc (size_t n);
/*@out@*/ float *sf_floatalloc (size_t n);
/*@out@*/ bool *sf_boolalloc (size_t n);
/*@out@*/ float complex *sf_complexalloc (size_t n);
/*@out@*/ float complex **sf_complexalloc2 (size_t n1, size_t n2);
/*@out@*/ float **sf_floatalloc2 (size_t n1, size_t n2);
/*@out@*/ float ***sf_floatalloc3 (size_t n1, size_t n2, size_t n3);
/*@out@*/ int **sf_intalloc2 (size_t n1, size_t n2);
/*@out@*/ char **sf_charalloc2 (size_t n1, size_t n2);

#endif
#ifndef _sf_c99_h
#define _sf_c99_h

/* Is it 1994 or 1999?? */
#if defined(__STDC__) && (__STDC_VERSION__ >= 199401L)

/* The following from C99 - must define for C90 */
#include <stdbool.h>       /* define bool, true, false */
#include <complex.h>

#else

typedef enum {false, true} bool;

/* What do we do with float functions? */
/* #define sqrtf sqrt 
#define logf log
#define expf exp
#define fabsf fabs
*/

/* What do we do with complex? */
#define complex  
#define I 0.0
#define csqrtf sqrtf
#define clogf logf
#define cexpf expf
#define cabsf fabsf
#define crealf sf_crealf
#define cimagf sf_cimagf

float sf_crealf(float complex c); 
float sf_cimagf(float complex c);

#endif

#endif
#ifndef _sf_decart_h
#define _sf_decart_h

/* index transform (vector to matrix) and its inverse */
void sf_line2cart( int dim, const int* nn, int i, int* ii);

int sf_cart2line( int dim, const int* nn, const int* ii);

#endif
#ifndef _sf_error_h
#define _sf_error_h

void sf_error( char *format, ... );
void sf_warning( char *format, ... );

#endif
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
bool sf_histfloat (sf_file file, const char* key,/*@out@*/ float* par);
char* sf_histstring (sf_file file, const char* key);
void sf_fileflush (sf_file file, sf_file src);
void sf_putint (sf_file file, char* key, int par);
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
#ifndef _sf_files_h
#define _sf_files_h

#include "c99.h"
#include "file.h"

int sf_filedims (sf_file file, /*@out@*/ int *n);
int sf_filesize (sf_file file);
int sf_leftsize (sf_file file, int dim);
void sf_cp(sf_file in, sf_file out);
void sf_rm(const char* filename, bool force, bool verb, bool inquire);

#endif
#ifndef _sf_getpar_h
#define _sf_getpar_h

#include <stdio.h>

#include "c99.h"

void sf_init(int argc,char *argv[]);
void sf_close (void);
char* sf_getprog (void);
char* sf_getuser (void);
char* sf_gethost (void);
bool sf_getint (const char* key,/*@out@*/ int* par);
bool sf_getints (const char* key,/*@out@*/ int* par,size_t n);
bool sf_getfloat (const char* key,/*@out@*/ float* par);
bool sf_getfloats (const char* key,/*@out@*/ float* par,size_t n);
char* sf_getstring (const char* key);
bool sf_getbool (const char* key,/*@out@*/ bool* par);
bool sf_getbools (const char* key,/*@out@*/ bool* par,size_t n);

#endif
#ifndef _sf_segy_h
#define _sf_segy_h

/* SEGY standard */

#define SF_SEGY_FORMAT  24
#define SF_SEGY_NS      20
#define SF_SEGY_DT      16

enum {
    SF_EBCBYTES=3200,	/* Bytes in the card image EBCDIC block */
    SF_BNYBYTES=400,	/* Bytes in the binary coded block	*/
    SF_HDRBYTES=240,	/* Bytes in the tape trace header	*/
    SF_NKEYS=71,	/* Number of mandated header fields	*/
    SF_BHKEYS=27	/* Number of mandated binary fields	*/
};

void sf_endian (void);
void sf_ebc2asc (int narr, char* arr);
void sf_asc2ebc (int narr, char* arr);
int sf_segyformat (const char* bhead);
int sf_segyns (const char* bhead);
float sf_segydt (const char* bhead);
void sf_segy2trace(const char* buf, float* trace, int ns, int format);
void sf_segy2head(const char* buf, int* head, int ns);
void sf_trace2segy(char* buf, const float* trace, int ns, int format);
void sf_head2segy(char* buf, const int* head, int ns);
int sf_segykey (const char* key);

#endif

#ifndef _sf_simtab_h
#define _sf_simtab_h

#include <stdio.h>

#include "c99.h"

/* Simbol Table structure */
typedef struct sf_SimTab *sf_simtab;

sf_simtab sf_simtab_init(int n);
void sf_simtab_close(sf_simtab table);
char *sf_simtab_get(sf_simtab table, const char *key);
bool sf_simtab_getint (sf_simtab table, const char* key,/*@out@*/ int* par);
bool sf_simtab_getints (sf_simtab table, const char* key,
			/*@out@*/ int* par,size_t n);
bool sf_simtab_getfloat (sf_simtab table, const char* key,/*@out@*/ float* par);
bool sf_simtab_getfloats (sf_simtab table, const char* key,
			  /*@out@*/ float* par,size_t n);
bool sf_simtab_getbool (sf_simtab table, const char* key,/*@out@*/ bool* par);
bool sf_simtab_getbools (sf_simtab table, const char* key,
			 /*@out@*/ bool* par,size_t n);
char* sf_simtab_getstring (sf_simtab table, const char* key);
void sf_simtab_enter(sf_simtab table, const char *key, const char* val);
void sf_simtab_put(sf_simtab table, const char *keyval);
void sf_simtab_input (sf_simtab table, FILE* fp);
void sf_simtab_output (sf_simtab table, FILE* fp);

#endif

#ifndef _sf_stack_h
#define _sf_stack_h

#include "c99.h"

typedef struct sf_Stack *sf_stack;

sf_stack sf_stack_init (int size);
int sf_stack_get (sf_stack s);
void sf_stack_set (sf_stack s, int pos);
void sf_push(sf_stack s, void *data, int type);
void* sf_pop(sf_stack s);
int sf_top(sf_stack s);
bool sf_full (sf_stack s);
void sf_stack_close(sf_stack s);
void sf_stack_print (sf_stack s);

#endif
