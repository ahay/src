#include "utils.h"

typedef struct {
  char * str;
} WORD;

typedef struct {
  WORD * key;
  WORD * val;
} KEYVAL;

typedef struct PSLINK {
  KEYVAL * pair;
  struct PSLINK * prev;
  struct PSLINK * next;
} PSLINK;

/*----------------------------------------------------------------------------*/
/** \defgroup create Creation and destruction */
/*@{*/
/** 
Initialze pair
*/
WORD * word_new();
void word_delete(WORD ** w);
void word_reset(WORD * w);
int word_assign(WORD * w, char * str, int len);
int word_whitechar(char c);
int word_copy(WORD * tgt, WORD src);
int word_read(WORD * w, char ** src);

KEYVAL * kv_new();
void kv_delete(KEYVAL ** pair);
void kv_reset(KEYVAL * pair);
int kv_check(KEYVAL src);
int kv_copy(KEYVAL * tgt, KEYVAL src);
int kv_read(KEYVAL * kv, char ** src);
void kv_print(KEYVAL kv);
void kv_fprint(KEYVAL kv, FILE * f);

PSLINK * pslink_new();
void pslink_delete(PSLINK ** p);
int pslink_front(PSLINK ** p);
int pslink_back(PSLINK ** p);
int pslink_read(PSLINK ** p, char ** str);
int pslink_findfirst(PSLINK * par, WORD skey, WORD * sval);
int pslink_findlast(PSLINK * par, WORD skey, WORD * sval);
int pslink_setfirst(PSLINK ** par, WORD skey, WORD sval);
//int pslink_setlast(PSLINK * par, WORD skey, WORD sval);

/*----------------------------------------------------------------------------*/
/** 
Set parameter array (all fields) to zeros - null initialization
@param[out] parr (PARARRAY *) - parameter array to initialize
@return (int) 0 if successful, else nonzero error code.
*/
//int ps_setnull(PARARRAY *parr);
/*----------------------------------------------------------------------------*/
/** 
Creates parameter array (STORAGE ALLOCATION) from file
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] filename (char *) - name of parfile containing key=value info
@return (int) 0 if successful, else nonzero error code.
*/
//int ps_createfile(PARARRAY *parr, const char *filename);
/** 
Creates parameter array (STORAGE ALLOCATION) from command-line argument list
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] argc (int) - arg count
@param[in] argv (char**) - arg array
@return (int) 0 if successful, else nonzero error code.
*/
//int ps_createargs(PARARRAY *parr, int argc, char **argv);
/** 
Creates parameter array (STORAGE ALLOCATION) from string.
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] str (SIDEDSTRING) - sized string struct from which to draw 
key=value pairs
@return (int) 0 if successful, else nonzero error code.
*/
//int ps_createstrz(PARARRAY *parr, SIZEDSTRING str);
/*----------------------------------------------------------------------------*/
/** 
Destroy parameter array (STORAGE DEALLOCATION). 
@param[out] parr (PARARRAY *) - param array destroyed on successful return
@return (int) 0 if successful, else nonzero error code.
*/
//int ps_destroy(PARARRAY *parr);
/*----------------------------------------------------------------------------*/
/** 
Concatenate parameter arrays (STORAGE REALLOCATION).
"parr = parr + parr2".
@param[out] parr (PARARRAY *) - first summand on call, sum on return
@param[in] parr2 (PARARRAY)  - second summand (const)
@return (int) 0 if successful, else nonzero error code.
*/
//int ps_concat(PARARRAY *parr, PARARRAY parr2);
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
//int ps_addpairstring(PARARRAY *parr, 
//		     const char * key, int mk,
//		     const char * value, int mv);

/*----------------------------------------------------------------------------*/
/**
Full copy function: lhs = rhs (STORAGE (RE)ALLOCATED)
//@param lhs (PARARRAY *) - target parameter array
//@param rhs (PARARRAY)   - source parameter array (const)
//@return (int) 0 if successful, else nonzero error code.
*/
//int ps_copy(PARARRAY rhs, PARARRAY *lhs);

/*@}*/
