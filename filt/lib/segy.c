/* Operations with SEGY standard files. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "segy.h"
#include "getpar.h"
#include "error.h"

#include "_bool.h"
/*^*/

#ifndef _sf_segy_h

#define SF_SEGY_FORMAT  24
#define SF_SEGY_NS      20
#define SF_SEGY_DT      16
/*^*/


enum {
    SF_EBCBYTES=3200,	/* Bytes in the card image EBCDIC block */
    SF_BNYBYTES=400,	/* Bytes in the binary coded block	*/
    SF_HDRBYTES=240,	/* Bytes in the tape trace header	*/
    SF_NKEYS=71,	/* Number of mandated header fields	*/
    SF_BHKEYS=27	/* Number of mandated binary fields	*/
};
/*^*/


#endif

typedef unsigned char byte;

static bool little_endian = false;

static byte EBCtoASC[256] = {
    0x00,0x01,0x02,0x03,0xCF,0x09,0xD3,0x7F,
    0xD4,0xD5,0xC3,0x0B,0x0C,0x0D,0x0E,0x0F,
    0x10,0x11,0x12,0x13,0xC7,0xB4,0x08,0xC9,
    0x18,0x19,0xCC,0xCD,0x83,0x1D,0xD2,0x1F,
    0x81,0x82,0x1C,0x84,0x86,0x0A,0x17,0x1B,
    0x89,0x91,0x92,0x95,0xA2,0x05,0x06,0x07,
    0xE0,0xEE,0x16,0xE5,0xD0,0x1E,0xEA,0x04,
    0x8A,0xF6,0xC6,0xC2,0x14,0x15,0xC1,0x1A,
    0x20,0xA6,0xE1,0x80,0xEB,0x90,0x9F,0xE2,
    0xAB,0x8B,0x9B,0x2E,0x3C,0x28,0x2B,0x7C,
    0x26,0xA9,0xAA,0x9C,0xDB,0xA5,0x99,0xE3,
    0xA8,0x9E,0x21,0x24,0x2A,0x29,0x3B,0x5E,
    0x2D,0x2F,0xDF,0xDC,0x9A,0xDD,0xDE,0x98,
    0x9D,0xAC,0xBA,0x2C,0x25,0x5F,0x3E,0x3F,
    0xD7,0x88,0x94,0xB0,0xB1,0xB2,0xFC,0xD6,
    0xFB,0x60,0x3A,0x23,0x40,0x27,0x3D,0x22,
    0xF8,0x61,0x62,0x63,0x64,0x65,0x66,0x67,
    0x68,0x69,0x96,0xA4,0xF3,0xAF,0xAE,0xC5,
    0x8C,0x6A,0x6B,0x6C,0x6D,0x6E,0x6F,0x70,
    0x71,0x72,0x97,0x87,0xCE,0x93,0xF1,0xFE,
    0xC8,0x7E,0x73,0x74,0x75,0x76,0x77,0x78,
    0x79,0x7A,0xEF,0xC0,0xDA,0x5B,0xF2,0xF9,
    0xB5,0xB6,0xFD,0xB7,0xB8,0xB9,0xE6,0xBB,
    0xBC,0xBD,0x8D,0xD9,0xBF,0x5D,0xD8,0xC4,
    0x7B,0x41,0x42,0x43,0x44,0x45,0x46,0x47,
    0x48,0x49,0xCB,0xCA,0xBE,0xE8,0xEC,0xED,
    0x7D,0x4A,0x4B,0x4C,0x4D,0x4E,0x4F,0x50,
    0x51,0x52,0xA1,0xAD,0xF5,0xF4,0xA3,0x8F,
    0x5C,0xE7,0x53,0x54,0x55,0x56,0x57,0x58,
    0x59,0x5A,0xA0,0x85,0x8E,0xE9,0xE4,0xD1,
    0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,
    0x38,0x39,0xB3,0xF7,0xF0,0xFA,0xA7,0xFF
};

static byte ASCtoEBC[256] = {
    0x00,0x01,0x02,0x03,0x37,0x2D,0x2E,0x2F,
    0x16,0x05,0x15,0x0B,0x0C,0x0D,0x0E,0x0F,
    0x10,0x11,0x12,0x13,0x3C,0x15,0x32,0x26,
    0x18,0x19,0x3F,0x27,0x1C,0x1D,0x1E,0x1F,
    0x40,0x5A,0x7F,0x7B,0x5B,0x6C,0x50,0x7D,
    0x4D,0x5D,0x5C,0x4E,0x6B,0x60,0x4B,0x61,
    0xF0,0xF1,0xF2,0xF3,0xF4,0xF5,0xF6,0xF7,
    0xF8,0xF9,0x7A,0x5E,0x4C,0x7E,0x6E,0x6F,
    0x7C,0xC1,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,
    0xC8,0xC9,0xD1,0xD2,0xD3,0xD4,0xD5,0xD6,
    0xD7,0xD8,0xD9,0xE2,0xE3,0xE4,0xE5,0xE6,
    0xE7,0xE8,0xE9,0xAD,0xE0,0xBD,0x5F,0x6D,
    0x79,0x81,0x82,0x83,0x84,0x85,0x86,0x87,
    0x88,0x89,0x91,0x92,0x93,0x94,0x95,0x96,
    0x97,0x98,0x99,0xA2,0xA3,0xA4,0xA5,0xA6,
    0xA7,0xA8,0xA9,0xC0,0x4F,0xD0,0xA1,0x07,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0xFF
};

struct segy {
    char *name;
    unsigned int size;
} segykey[] = {
    {"tracl",  4},  /* trace sequence number within line 0 */

    {"tracr",  4},  /* trace sequence number within reel 4 */

    {"fldr",   4},  /* field record number 8 */

    {"tracf",  4},  /* trace number within field record 12 */

    {"ep",     4},  /* energy source point number 16 */

    {"cdp",    4},  /* CDP ensemble number 20 */

    {"cdpt",   4},  /* trace number within CDP ensemble 24 */

    {"trid",   2},  /* trace identification code:
                        1 = seismic data
                        2 = dead
                        3 = dummy
                        4 = time break
                        5 = uphole
                        6 = sweep
                        7 = timing
                        8 = water break
                        9---, N = optional use (N = 32,767) 28 */

    {"nvs",    2},  /* number of vertically summed traces (see
		    vscode in bhed structure) 30 */

    {"nhs",    2},  /* number of horizontally summed traces (see
		    vscode in bhed structure) 32 */

    {"duse",   2},  /* data use:
		    1 = production
		    2 = test 34 */

    {"offset", 4},  /* distance from source point to receiver
		    group (negative if opposite to direction
		    in which the line was shot) 36 */

    {"gelev",  4},  /* receiver group elevation from sea level
		    (above sea level is positive) 40 */

    {"selev",  4},  /* source elevation from sea level
		     (above sea level is positive) 44 */

    {"sdepth", 4},  /* source depth (positive) 48 */

    {"gdel",   4},  /* datum elevation at receiver group 52 */

    {"sdel",   4},  /* datum elevation at source 56 */

    {"swdep",  4},  /* water depth at source 60 */

    {"gwdep",  4},  /* water depth at receiver group 64 */

    {"scalel", 2},  /* scale factor for previous 7 entries
		     with value plus or minus 10 to the
		     power 0, 1, 2, 3, or 4 (if positive,
		     multiply, if negative divide) 68 */

    {"scalco", 2},  /* scale factor for next 4 entries
		     with value plus or minus 10 to the
		     power 0, 1, 2, 3, or 4 (if positive,
		     multiply, if negative divide) 70 */

    {"sx",     4},  /* X source coordinate 72 */

    {"sy",     4},  /* Y source coordinate 76 */

    {"gx",     4},  /* X group coordinate 80 */

    {"gy",     4},  /* Y source coordinate 84 */

    {"counit", 2},  /* coordinate units code:
		     for previoius four entries
		     1 = length (meters or feet)
		     2 = seconds of arc (in this case, the
		     X values are unsigned intitude and the Y values
		     are latitude, a positive value designates
		     the number of seconds east of Greenwich
		     or north of the equator 88 */

    {"wevel",   2},  /* weathering velocity 90 */

    {"swevel",  2},  /* subweathering velocity 92 */

    {"sut",     2},  /* uphole time at source 94 */

    {"gut",     2},  /* uphole time at receiver group 96 */

    {"sstat",   2},  /* source static correction 98 */

    {"gstat",   2},  /* group static correction 100 */

    {"tstat",   2},  /* total static applied 102 */

    {"laga",    2},  /* lag time A, time in ms between end of 240-
		      byte trace identification header and time
		      break, positive if time break occurs after
		      end of header, time break is defined as
		      the initiation pulse which maybe recorded
		      on an auxiliary trace or as otherwise
		      specified by the recording system 104 */

    {"lagb",    2},  /* lag time B, time in ms between the time
		      break and the initiation time of the energy source,
		      may be positive or negative 106 */

    {"delrt",   2},  /* delay recording time, time in ms between
		      initiation time of energy source and time
		      when recording of data samples begins
		      (for deep water work if recording does not
		      start at zero time) 108 */

    {"muts",    2},  /* mute time--start 110 */

    {"mute",    2},  /* mute time--end 112 */

    {"ns",      2},  /* number of samples in this trace 114 */

    {"dt",      2},  /* sample interval, in micro-seconds 116 */

    {"gain",    2},  /* gain type of field instruments code:
		      1 = fixed
		      2 = binary
		      3 = floating point
		      4 ---- N = optional use 118 */

    {"igc",    2},   /* instrument gain constant 120 */

    {"igi",    2},   /* instrument early or initial gain 122 */

    {"corr",   2},   /* correlated:
		      1 = no
		      2 = yes 124 */    

    {"sfs",    2},   /* sweep frequency at start 126 */

    {"sfe",    2},   /* sweep frequency at end 128 */

    {"slen",   2},   /* sweep length in ms 130 */

    {"styp",   2},   /* sweep type code:
		      1 = linear
		      2 = cos-squared
		      3 = other 132 */   

    {"stas",   2},   /* sweep trace length at start in ms 134 */

    {"stae",   2},   /* sweep trace length at end in ms 136 */

    {"tatyp",  2},   /* taper type: 1=linear, 2=cos^2, 3=other 138 */

    {"afilf",  2},   /* alias filter frequency if used 140 */

    {"afils",  2},   /* alias filter slope 142 */

    {"nofilf", 2},   /* notch filter frequency if used 144 */

    {"nofils", 2},   /* notch filter slope 146 */

    {"lcf",    2},   /* low cut frequency if used 148 */

    {"hcf",    2},   /* high cut frequncy if used 150 */

    {"lcs",    2},   /* low cut slope 152 */

    {"hcs",    2},   /* high cut slope 154 */

    {"year",   2},   /* year data recorded 156 */

    {"day",    2},   /* day of year 158 */

    {"hour",   2},   /* hour of day (24 hour clock) 160 */

    {"minute", 2},   /* minute of hour 162 */

    {"sec",    2},   /* second of minute 164 */

    {"timbas", 2},   /* time basis code:
		      1 = local
		      2 = GMT
		      3 = other 166 */   

    {"trwf",   2},   /* trace weighting factor, defined as 1/2^N
		      volts for the least sigificant bit 168 */

    {"grnors", 2},   /* geophone group number of roll switch
		      position one 170 */

    {"grnofr", 2},   /* geophone group number of trace one within
		      original field record 172 */

    {"grnlof", 2},   /* geophone group number of last trace within
		      original field record 174 */

    {"gaps",   2},   /* gap size (total number of groups dropped) 176 */

    {"otrav",  2},   /* overtravel taper code:
		      1 = down (or behind)
		      2 = up (or ahead) 71/178 */
    /* 72/180 73/184 74/188 75/192 76/196 77/200 78/204 79/208 
       80/212 81/216 82/220 83/224 84/228 85/232 86/236 */
};


/* Big-endian to Little-endian conversion and back */
static int convert2(const char* buf);
static int convert4(const char* buf);
static void insert2(int y, char* buf);
static void insert4(int y, char* buf);
static void swapb(byte *x, byte *y);

/* IBM to IEEE float conversion and back */
static float ibm2float (const char* num);
static void float2ibm (float y, char* num);

static void swapb(byte *x, byte *y) 
/* swap two bytes */
{
    byte tmp; 

    tmp = *x; 
    *x = *y; 
    *y = tmp;
}

static int convert2(const char* buf)
/* convert buf to 2-byte int */
{
    union {
	byte b[2];
	short s;
    } x;

    memcpy(x.b,buf,2);

    if (little_endian) swapb(x.b,x.b+1);

    return (int) x.s;
}

static void insert2(int y, char* buf)
/* convert 2-byte int to buf */
{
    union {
	byte b[2];
	short s;
    } x;

    x.s=y;

    if (little_endian) swapb(x.b,x.b+1);
    
    memcpy(buf,x.b,2);
}

static int convert4(const char* buf)
/* convert buf to 4-byte int */
{
    union {
	byte b[4];
	long s;
    } x;

    memcpy(x.b,buf,4);

    if (little_endian) {
	swapb(x.b,x.b+3);
	swapb(x.b+1,x.b+2);
    }

    return (int) x.s;
}

static void insert4(int y, char* buf)
/* convert 4-byte int to buf */
{
    union {
	byte b[4];
	long s;
    } x;

    x.s=y;

    if (little_endian) {
	swapb(x.b,x.b+3);
	swapb(x.b+1,x.b+2);
    }
    
    memcpy(buf,x.b,4);
}

bool sf_endian (void)
/*< Endianness test, returns true for little-endian machines >*/
{
    union {
	byte c[4];
	long i;
    } test;

    test.i=0;
    test.c[0]=1;
    
    assert (2 == sizeof(short) && 4 == sizeof(long)); /* fix this later */
    little_endian = (0 != (test.i << 8));
    
    return little_endian;
}

void sf_ebc2asc (int narr, char* arr)
/*< Convert char array arrr[narr]: EBC to ASCII >*/
{
    int i, j;

    for (i=0; i < narr; i++) {
	j = (unsigned char) arr[i];
	arr[i] = (char) EBCtoASC[j];
    }
}

void sf_asc2ebc (int narr, char* arr)
/*< Convert char array arrr[narr]: ASCII to EBC >*/
{
    int i, j;

    for (i=0; i < narr; i++) {
	j = (unsigned char) arr[i];
	arr[i] = (char) ASCtoEBC[j];
    }
}


int sf_segyformat (const char* bhead)
/*< extracts SEGY format from binary header >*/
{
    return convert2(bhead+SF_SEGY_FORMAT);
}

int sf_segyns (const char* bhead)
/*< extracts ns (number of samples) from binary header >*/
{
    return convert2(bhead+SF_SEGY_NS);
}

float sf_segydt (const char* bhead)
/*< extracts dt (sampling) from binary header >*/
{
    return (float) (convert2(bhead+SF_SEGY_DT)/1000000.);
}

static void float2ibm (float y, char* num)
/* floating-point conversion to IBM format */
{
    unsigned long x, s, f;
    long e;
    const unsigned long fMAXIBM = 0x7FFFFFFF;

    memcpy (&x,&y,4);

    /* check for special case of zero */
    if ((x & 0x7fffffff) == 0) {
	insert4(x,num);
	return; 
    }

    /* fetch the sign, exponent (removing excess 127), and fraction */
    s =   x & 0x80000000;
    e = ((x & 0x7f800000) >> 23) - 127;
    f =   x & 0x007fffff;

    /* convert 23 bit fraction to 24 bit fraction */
    f <<= 1; 

    /* restore the '1' preceeded the IEEE binary point */
    f |= 0x01000000; 
    
    /* convert scale factor from base-2 to base-16 */
    if (e >= 0) {
      f <<= (e & 3); 
      e >>= 2;
    } else {
	f >>= ((-e) & 3); 
	e = -((-e) >> 2);
    }

    /* reduce fraction to 24 bits */
    if (f & 0x0f000000) {
	f >>= 4;
	e += 1;
    }

    /* convert exponent to excess 64 and store the number */
    if ((e += 64) > 127) {
	s |= fMAXIBM;
    } else if (e >= 0) {
	s |= (e << 24) | f;
    }

    insert4(s,num);
}

static float ibm2float (const char* num)
/* floating point conversion from IBM format */
{
    unsigned long x, s, f;
    const unsigned long fMAXIEEE = 0x7F7FFFFF;
    long e;         
    float y;
                                                                     
    x = convert4 (num);
    
    /* check for special case of zero */
    if ((x & 0x7fffffff) == 0) return 0.0; 

    /* fetch the sign, exponent (removing excess 64), and fraction */   
    s =   x & 0x80000000;                                               
    e = ((x & 0x7f000000) >> 24) - 64;                                   
    f =   x & 0x00ffffff;                                                
                                                                    
    /* convert scale factor from base-16 to base-2 */        
    if (e >= 0) {
	e <<= 2;  
    } else { 
	e = -((-e) << 2); 
    }
                                                                        
    /* convert exponent for 24 bit fraction to 23 bit fraction */           
    e -= 1;                                                               
                                                                            
    /* normalize the fraction */                                            
    if (0 != f) {
	while ((f & 0x00800000) == 0) {         
	    f <<= 1;
	    e -= 1;
	}       
    }                                                               
	
    /* drop the '1' preceeding the binary point */                       
    f &= 0x007fffff;                                                         
    
    /* convert exponent to excess 127 and store the number */
    if ((e += 127) >= 255) {
	s |= fMAXIEEE;
    } else if (e > 0) {
	s |= (e << 23) | f; 	    
    }    

    memcpy (&y,&s,4);
    return y;
}

void sf_segy2trace(const char* buf, float* trace, int ns, int format)
/*< Extract a floating-point trace[nt] from buffer buf.
---
format: 1: IBM, 2: int4, 3: int2
>*/
{
    int i;

    for (i=0; i < ns; i++, buf += 4) {
	switch (format) {
	    case 1: trace[i] = ibm2float (buf);       break; /* IBM float */
	    case 2: trace[i] = (float) convert4(buf); break; /* int4 */
	    case 3: trace[i] = (float) convert2(buf); break; /* int2 */
	    default: sf_error("Unknown format %d",format); break;
	}
    }
}

void sf_trace2segy(char* buf, const float* trace, int ns, int format)
/*< Convert a floating-point trace[ns] to buffer buf.
---
format: 1: IBM, 2: int4, 3: int2
>*/
{
    int i;

    for (i=0; i < ns; i++, buf += 4) {
	switch (format) {
	    case 1: float2ibm(trace[i],buf);     break; /* IBM float */
	    case 2: insert4((int) trace[i],buf); break; /* int4 */
	    case 3: insert2((int) trace[i],buf); break; /* int2 */
	    default: sf_error("Unknown format %d",format); break;
	}
    }
}

void sf_segy2head(const char* buf, int* trace, int nk)
/*< Create an integer trace header trace[nk] from buffer buf >*/
{
    int i, byte;
    const char *buf0, *bufi;

    buf0 = buf;
    for (i=0; i < nk; i++) {
	/* allow to remap header keys */
	bufi = sf_getint(segykey[i].name,&byte)? buf0+byte:buf;

	switch (segykey[i].size) {
	    case 2:
		trace[i] = convert2(bufi);
		buf += 2;
		break;
	    case 4:
		trace[i] = convert4(bufi);
		buf += 4;
		break;
	    default:
		sf_error("Unknown size %d",segykey[i].size);
		break;
	}
    }
}

int sf_segykey (const char* key) 
/*< Extract a SEGY key value >*/
{
    
    int i;

    for (i=0; i < SF_NKEYS; i++) {
	if (0==strcmp(key,segykey[i].name)) return i;
    }
    sf_error("no such key %s",key);
    return 0;
}

char* sf_segykeyword (int k) 
/*< Find a SEGY key from its number >*/
{
    return segykey[k].name;
}

void sf_head2segy(char* buf, const int* trace, int nk)
/*< Convert an integer trace[nk] to buffer buf >*/
{
    int i;

    for (i=0; i < nk; i++) {
	if (i < SF_NKEYS && 2 == segykey[i].size) {
	    insert2(trace[i],buf);
	    buf += 2;
	} else {
	    insert4(trace[i],buf);
	    buf += 4;
	}
    }
}

/* 	$Id$	 */
