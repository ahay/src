#ifndef _sf_segy_h
#define _sf_segy_h

#include "_bool.h"

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

bool sf_endian (void);
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
char* sf_segykeyword (int k);

#endif

/* 	$Id$	 */
