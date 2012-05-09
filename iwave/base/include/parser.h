/** @page parser IWAVE Parameter Parser 
    Author: Igor Terentyev
    <p>
    file: \ref parser.h
    <p>
    Decomposes file/string into pairs of strings: [key][value], and provides
    functions to convert value to prescribed type.
    Separator control symbol between key and value is '='
    (can be redefined in \ref utils.h file).
    Delimiters between key=value pairs are any non-graphical symbols.
    Keys CANNOT be empty: if key not found between separators,
    separator is ignored (treated as delimiter). Values can be empty.
    Quote control symbol is '"' (can be redefined in utils.h file).
    Everything between quotes is treated as graphical symbol.
    '"' symbol (non-control) is given by [""]. 
    Strings do not have to be null-terminated. String size is provided.
    For just in case strings produced by parser are by one symbol longer
    then their ireal sizes and have null-terminator in the extra space.
    <p>
    EXAMPLES: on left of "==>", input string (in square brackets). On right, [key] [value] produced by parser.
    <ul>
    <li>
    [a=b]
    ==> [a][b].
    </li>

    <li>
    [a = b]
    ==> [a][b].
    </li>

    <li>
    [a =b = c = d]
    ==> [a][b],[c][d]; second '=' ignored since no key to the left from it.
    </li>

    <li>
    ["  key w/ spaces" = "value w/ spaces and ""quotes"""]
    ==> [  key w/ spaces][value w/ spaces and "quotes"]
    </li>

    <li>
    [eq_sign="="]
    ==> [eq_sign][=].
    </li>

    <li>
    [a = 1 comment 1 b = 2 comment 2 c = 3]
    ==> [a][1],[b][2],[c][3].
    </li>

    <li>
    [ a

    <p>

    =

    <p>

    1]
    ==> [a][1].
    </li>

    <li>
    [a= = (empty parameter) b = 1 (non-empty parameter)]
    ==> [a][],[b][1].
    </li>
    </ul>
    <p>
    The package provides several types of functions which manipulate parameter
    arrays and access parameter values:
    <ol>
    <li>\ref create</li>
    <li>\ref print</li>
    <li>\ref ffaccess</li>
    <li>\ref flaccess</li>
    <li>\ref utility</li>
    </ol>
*/
/*============================================================================*/

#ifndef __PARSER_H_
#define __PARSER_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"

/*----------------------------------------------------------------------------*/
/*
Parameter structure.

int npar        :  number of parameters in the array.
SIZEDSTRING *ns :  keys array.
SIZEDSTRING *vs :  values array.
char *buffer    :  allocated buffer.
*/
typedef struct
{
  int npar;
  SIZEDSTRING *ns;
  SIZEDSTRING *vs;
  char *buffer;
	int buffersize;
} PARARRAY;
/*----------------------------------------------------------------------------*/
/** \defgroup create Creation and destruction */
/*@{*/
/** 
Set parameter array (all fields) to zeros - null initialization
@param[out] parr (PARARRAY *) - parameter array to initialize
@return (int) 0 if successful, else nonzero error code.
*/
int ps_setnull(PARARRAY *parr);
/*----------------------------------------------------------------------------*/
/** 
Creates parameter array (STORAGE ALLOCATION) from file
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] filename (char *) - name of parfile containing key=value info
@return (int) 0 if successful, else nonzero error code.
*/
int ps_createfile(PARARRAY *parr, const char *filename);
/** 
Creates parameter array (STORAGE ALLOCATION) from command-line argument list
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] argc (int) - arg count
@param[in] argv (char**) - arg array
@return (int) 0 if successful, else nonzero error code.
*/
int ps_createargs(PARARRAY *parr, int argc, char **argv);
/** 
Creates parameter array (STORAGE ALLOCATION) from string.
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] str (SIDEDSTRING) - sized string struct from which to draw 
key=value pairs
@return (int) 0 if successful, else nonzero error code.
*/
int ps_createstrz(PARARRAY *parr, SIZEDSTRING str);
/*----------------------------------------------------------------------------*/
/** 
Destroy parameter array (STORAGE DEALLOCATION). 
@param[out] parr (PARARRAY *) - param array destroyed on successful return
@return (int) 0 if successful, else nonzero error code.
*/
int ps_destroy(PARARRAY *parr);
/*----------------------------------------------------------------------------*/
/** 
Concatenate parameter arrays (STORAGE REALLOCATION).
"parr = parr + parr2".
@param[out] parr (PARARRAY *) - first summand on call, sum on return
@param[in] parr2 (PARARRAY)  - second summand (const)
@return (int) 0 if successful, else nonzero error code.
*/
int ps_concat(PARARRAY *parr, PARARRAY parr2);
/*----------------------------------------------------------------------------*/
/** 
Add "key=value" line from string input (STORAGE REALLOCATION).
@param[out] parr (PARARRAY *) - contains additional key=val pair on return
@param[in] key (const char *)  - null-terminated char array = key
@param[in] mk (int) - allocated size of key, should be >= strlen(key)+1;
@param[in] value (const char *) - null-terminated char array = value
@param[in] mv (int) - allocated size of value, should be >= strlen(value)+1;
@return (int) 0 if successful, else nonzero error code.
*/
int ps_addpairstring(PARARRAY *parr, 
		     const char * key, int mk,
		     const char * value, int mv);

/*----------------------------------------------------------------------------*/
/**
Full copy function: lhs = rhs (STORAGE (RE)ALLOCATED)
@param lhs (PARARRAY *) - target parameter array
@param rhs (PARARRAY)   - source parameter array (const)
@return (int) 0 if successful, else nonzero error code.
*/
int ps_copy(PARARRAY rhs, PARARRAY *lhs);

/*@}*/

/*----------------------------------------------------------------------------*/
/** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
*/

/*@{*/
int ps_ffcstring(PARARRAY parr, const char *key, char **p);
int ps_ffchar(PARARRAY parr, const char *key, char *p);
int ps_ffshort(PARARRAY parr, const char *key, short *p);
int ps_ffint(PARARRAY parr, const char *key, int *p);
int ps_fflong(PARARRAY parr, const char *key, long *p);
int ps_ffushort(PARARRAY parr, const char *key, unsigned short *p);
int ps_ffuint(PARARRAY parr, const char *key, unsigned int *p);
int ps_ffulong(PARARRAY parr, const char *key, unsigned long *p);
int ps_fffloat(PARARRAY parr, const char *key, float *p);
int ps_ffdouble(PARARRAY parr, const char *key, double *p);
int ps_ffreal(PARARRAY parr, const char *key, ireal *p);
/*@}*/

/*----------------------------------------------------------------------------*/
/** \defgroup flaccess Parameter access - last occurence
 
Get first occurence of data in the last occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
*/

/*@{*/
int ps_flcstring(PARARRAY parr, const char *key, char **p);
int ps_flchar(PARARRAY parr, const char *key, char *p);
int ps_flshort(PARARRAY parr, const char *key, short *p);
int ps_flint(PARARRAY parr, const char *key, int *p);
int ps_fllong(PARARRAY parr, const char *key, long *p);
int ps_flushort(PARARRAY parr, const char *key, unsigned short *p);
int ps_fluint(PARARRAY parr, const char *key, unsigned int *p);
int ps_flulong(PARARRAY parr, const char *key, unsigned long *p);
int ps_flfloat(PARARRAY parr, const char *key, float *p);
int ps_fldouble(PARARRAY parr, const char *key, double *p);
int ps_flreal(PARARRAY parr, const char *key, ireal *p);
/*@}*/

/** \defgroup print Output to stream
 */
/*@{*/
/**
Write contents of PARARRAY to stream
@param[in] parr (PARARRAY) - input param array
@param[out] stream (FILE *) - output stream
@return (int) 0 if successful, else nonzero error code.
*/

int ps_printall(PARARRAY parr, FILE *stream);
/*@}*/

/*----------------------------------------------------------------------------*/
/** \defgroup utility Utility functions - mostly internal use */
/*@{*/
/* 
Return number of parameter key occurences. 
*/
int ps_countnamez(PARARRAY parr, SIZEDSTRING key);
int ps_countname(PARARRAY parr, const char *key);
/*----------------------------------------------------------------------------*/
/*
Return number of parameter type value occurences in the parameter value. 

int occ :  key occurence number.
*/
int ps_counttypez(PARARRAY parr, SIZEDSTRING key, int occ, int type);
int ps_counttype(PARARRAY parr, const char *key, int occ, int type);
/*----------------------------------------------------------------------------*/
/*
Extract type value.

int occ :  key occurence number.
int n   :  type occurence number in the parameter value.
*/
int ps_gettypez(PARARRAY parr, SIZEDSTRING key, int occ, int type, int n, void *p);
int ps_gettype(PARARRAY parr, const char *key, int occ, int type, int n, void *p);
/*----------------------------------------------------------------------------*/
/* 
Returns sized string value of the parameter from key.

int occ :  key occurence number.
*/
int ps_getvalz(PARARRAY parr, SIZEDSTRING key, int occ, SIZEDSTRING *val);
int ps_getval(PARARRAY parr, const char *key, int occ, SIZEDSTRING *val);
/*----------------------------------------------------------------------------*/
/*
Return index of parameter occurence. 
Returns -1 if the occurence not found.

int occ :  key occurence number.
*/
int ps_nameindexz(PARARRAY parr, SIZEDSTRING key, int occ);
int ps_nameindex(PARARRAY parr, const char *key, int occ);
/*----------------------------------------------------------------------------*/
/*
Convert zero terminated string to sized string (NO ALLOCATION).
*/
SIZEDSTRING ps_name2z(const char *key);
/*----------------------------------------------------------------------------*/
/*
Seeks next value and converts it.
MAY ALLOCATE (if type is cstring).
*/
int ps_val2type(SIZEDSTRING val, int type, void *p, SIZEDSTRING *valend);
/*----------------------------------------------------------------------------*/
/*@}*/
#endif /*__PARSER_H_*/
