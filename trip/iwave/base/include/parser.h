#ifndef __PARSER_H_
#define __PARSER_H_

#include "std_cpp_includes.hh"
#include "except.hh"
#include "utils.h"
#include "iwave_fopen.h"

typedef struct s_WORD {
  char * str;
} WORD;

typedef struct s_KEYVAL {
  WORD * key;
  WORD * val;
} KEYVAL;

typedef struct s_PSLINK {
  KEYVAL * pair;
  struct s_PSLINK * prev;
  struct s_PSLINK * next;
} PSLINK;

typedef struct s_PARARRAY {
  PSLINK * list;
} PARARRAY;

WORD * word_new();
void word_delete(WORD ** w);
void word_reset(WORD * w);
int word_assign(WORD * w, const char * str, int len);
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
void pslink_setnull(PSLINK ** p);
int pslink_front(PSLINK ** p);
int pslink_back(PSLINK ** p);
int pslink_read(PSLINK ** p, char ** str);
int pslink_findfirst(PSLINK * par, WORD skey, WORD * sval);
int pslink_findlast(PSLINK * par, WORD skey, WORD * sval);
int pslink_setfirst(PSLINK ** par, WORD skey, WORD sval);
int pslink_setlast(PSLINK ** par, WORD skey, WORD sval);

/*----------------------------------------------------------------------------*/
/** \defgroup create PARARRAY: Creation and destruction */
/*@{*/
/** 
Default - no keyval pairs, only one link
@return "bare" PARARRAY
*/
PARARRAY * ps_new();
/*----------------------------------------------------------------------------*/
/**
Destroys all links, frees all memory associated with object, returns *p = NULL
@param[out] pointer to object address
*/
void ps_delete(PARARRAY ** p);
/*----------------------------------------------------------------------------*/
/** 
Deletes, reinitializes list - provided for backward compatibility
@param[out] parr (PARARRAY *) - parameter array to initialize
@return (int) 0 if successful, else nonzero error code.
*/
int ps_setnull(PARARRAY *parr);
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

/*----------------------------------------------------------------------------*/
/** 
Creates parameter array (STORAGE ALLOCATION) from by reading from 
file pointer
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] fp (FILE *) - file pointer
@return (int) 0 if successful, else nonzero error code.
*/
int ps_createfp(PARARRAY *parr, FILE * fp);

/*----------------------------------------------------------------------------*/
/** 
Creates parameter array (STORAGE ALLOCATION) from file. Returns control of
file pointer, which may be allocated via iwave_fopen (or simply returned, if
file is already open). Includes optional specification of file prototype.
@param[out] parr - param array created on successful return
@param[in/out] stream - pointer to FILE pointer, which may be alloc
@param[in] proto - optional prototype filename for use by iwave_fopen, or NULL
@param[in] filename - name of parfile containing key=value info
@return (int) 0 if successful, else nonzero error code.
*/
int ps_createfile_fproto(PARARRAY *parr, 
			 FILE ** stream,  
			 const char * proto,
			 const char *filename);

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
NOTE: if "par" is amongst keys, create additional parameters (key=value pairs) 
from file whose name is value for "par".
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] argc (int) - arg count
@param[in] argv (char**) - arg array
@return (int) 0 if successful, else nonzero error code.
*/
int ps_createargs(PARARRAY *parr, int argc, char **argv);

/*----------------------------------------------------------------------------*/
/**
Full copy function: lhs = rhs (STORAGE (RE)ALLOCATED)
//@param lhs (PARARRAY **) - target parameter array
//@param rhs (PARARRAY)   - source parameter array (const)
//@return (int) 0 if successful, else nonzero error code.
*/
int ps_copy(PARARRAY ** tgt, PARARRAY src);

/*@}*/

/*----------------------------------------------------------------------------*/
/** \defgroup ffaccess Parameter access - first occurence 

Gets value defined in the first occurence of the key. 

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

Get value defined in the last occurence of the key. 

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

/*----------------------------------------------------------------------------*/
/** \defgroup sfassign Parameter assignment - first occurence 

Set value defined in the first occurence of the key, or add key=value pair
to list if key not already present. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.
*/

/*@{*/
int ps_sfcstring(PARARRAY parr, const char *key, const char *val);
int ps_sfchar(PARARRAY parr, const char *key, char p);
int ps_sfshort(PARARRAY parr, const char *key, short p);
int ps_sfint(PARARRAY parr, const char *key, int p);
int ps_sflong(PARARRAY parr, const char *key, long p);
int ps_sfushort(PARARRAY parr, const char *key, unsigned short p);
int ps_sfuint(PARARRAY parr, const char *key, unsigned int p);
int ps_sfulong(PARARRAY parr, const char *key, unsigned long p);
int ps_sffloat(PARARRAY parr, const char *key, float p);
int ps_sfdouble(PARARRAY parr, const char *key, double p);
int ps_sfreal(PARARRAY parr, const char *key, ireal p);
/*@}*/

/*----------------------------------------------------------------------------*/
/** \defgroup sfassign Parameter assignment - last occurence 

Set value defined in the last occurence of the key, or add key=value pair
to list if key not already present. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.
*/

/*@{*/
int ps_slcstring(PARARRAY parr, const char *key, const char *val);
int ps_slchar(PARARRAY parr, const char *key, char p);
int ps_slshort(PARARRAY parr, const char *key, short p);
int ps_slint(PARARRAY parr, const char *key, int p);
int ps_sllong(PARARRAY parr, const char *key, long p);
int ps_slushort(PARARRAY parr, const char *key, unsigned short p);
int ps_sluint(PARARRAY parr, const char *key, unsigned int p);
int ps_slulong(PARARRAY parr, const char *key, unsigned long p);
int ps_slfloat(PARARRAY parr, const char *key, float p);
int ps_sldouble(PARARRAY parr, const char *key, double p);
int ps_slreal(PARARRAY parr, const char *key, ireal p);
/*@}*/

/*----------------------------------------------------------------------------*/
/** \defgroup overloaded parse functions

 */
namespace RVL {
  bool parse(PARARRAY const & par, std::string name, string & val);
  bool parse(PARARRAY const & par, std::string name, char & val);
  bool parse(PARARRAY const & par, std::string name, short & val);
  bool parse(PARARRAY const & par, std::string name, int & val);
  bool parse(PARARRAY const & par, std::string name, long & val);
  bool parse(PARARRAY const & par, std::string name, unsigned short & val);
  bool parse(PARARRAY const & par, std::string name, unsigned int & val);
  bool parse(PARARRAY const & par, std::string name, unsigned long & val);
  bool parse(PARARRAY const & par, std::string name, float & val);
  bool parse(PARARRAY const & par, std::string name, double & val);
  bool parse(PARARRAY const & par, std::string name, bool & val);

  /** assigns default if returns false */
  template<typename T> T valparse(PARARRAY const & par, std::string name, T def) {
    T val;
    if (!parse(par,name,val)) val=def;
    return val;
  }
  /** throws exception if returns false */
  template<typename T> T valparse(PARARRAY const & par, std::string name) {
    T val;
    if (!parse(par,name,val)) {
      RVLException e;
      e<<"Error: parse_except\n";
      e<<"  failed to parse key = "<<name<<"\n";
      throw e;
    }
    return val;
  }

}

#endif
