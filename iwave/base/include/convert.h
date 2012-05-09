/*
convert.h
Igor Terentyev.
********************************************************************************
Functions to convert string to a value.

Unlike strtod (and others):
  1) Skips everything (not only white spaces), until finds number start.
  2) Does not change _p, if no conversion can be performed.
  3) In over/underflow cases conversion is performed and corresponding flag.
  
Examples:

strz2double of "abc!2.34Qqqq" will return 2.34 and new string "Qqqq"
--------------------------------------------------------------------------------

strz2int of "abc!2.34Qqqq" will return 2 and new string ".34Qqqq".
Subsequent strz2int will return 34.
*/
/*============================================================================*/

#ifndef __CONVERT_H_
#define __CONVERT_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
/*----------------------------------------------------------------------------*/
/*
SIZEDSTRING str     :  input string.
type *p             :  output number (can be NULL).
SIZEDSTRING *strend :  string starting after last charachter used for conversion.

int return          :  conversion error or over/underflow flag.
*/
int strz2char(SIZEDSTRING str, char *p, SIZEDSTRING *strend);

int strz2short(SIZEDSTRING str, short *p, SIZEDSTRING *strend);
int strz2int(SIZEDSTRING str, int *p, SIZEDSTRING *strend);
int strz2long(SIZEDSTRING str, long *p, SIZEDSTRING *strend);

int strz2ushort(SIZEDSTRING str, unsigned short *p, SIZEDSTRING *strend);
int strz2uint(SIZEDSTRING str, unsigned int *p, SIZEDSTRING *strend);
int strz2ulong(SIZEDSTRING str, unsigned long *p, SIZEDSTRING *strend);

int strz2float(SIZEDSTRING str, float *p, SIZEDSTRING *strend);
int strz2double(SIZEDSTRING str, double *p, SIZEDSTRING *strend);
/*----------------------------------------------------------------------------*/

#endif /*__CONVERT_H_*/
