/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/* HDRPKGE: $Revision: 1.11 $ ; $Date: 1996/09/09 17:13:34 $	*/

/*********************** self documentation **********************/
/************************************************************************** 
HDRPKGE - routines to access the SEGY header via the hdr structure.
 
gethval		get a trace header word by index
puthval		put a trace header word by index
getbhval	get a binary header word by index
putbhval	put a binary header word by index
gethdval	get a trace header word by name
puthdval	put a trace header word by name


hdtype		get the data type of a trace header word by name
getkey		get the name of a trace header word from its index

getindex	get the index of a trace header word from the name

swaphval	swap the trace header words by index
swapbhval	swap the binary header words by index
gettapehval	get a tape trace header word by index
puttapehval	put a tape trace header word by index
gettapebhval	get a tape binary header word by index
puttapebhval	put a tape binary header word by index
printheader	display non-null header field values

*************************************************************************** 
Function Prototypes:
void gethval(const segy *tr, int index, Value *valp);
void puthval(segy *tr, int index, Value *valp);
void putbhval(bhed *bh, int index, Value *valp);
void getbhval(const bhed *bh, int index, Value *valp);
void gethdval(const segy *tr, char *key, Value *valp);
void puthdval(segy *tr, char *key, Value *valp);
char *hdtype(const char *key);
char *getkey(const int index);
int getindex(const char *key);
void swaphval(segy *tr, int index);
void swapbhval(bhed *bh, int index);
void gettapehval(tapesegy *tapetr, int index, Value *valp);
void puttapehval(tapesegy *tapetr, int index, Value *valp);
void gettapebhval(tapebhed *tapetr, int index, Value *valp);
void puttapebhval(tapebhed *tapetr, int index, Value *valp);
void printheader(const segy *tp);

*************************************************************************** 
Notes:
This package includes only those routines that directly access
the "hdr" or "bhdr" structures.  It does not include routines
such as printfval, printftype, printfhead that use the routines
in this package to indirectly access these structures.

Note that while gethdval and puthdval are more convenient to use
than gethval and puthval, they incur an inefficiency in the
common case of iterating code over a set of traces with a fixed
key or keys.  In such cases, it is advisable to set the index
or indices outside the loop using getindex.

swaphval:
Byte-swapping is needed for converting SU data from big-endian to little-
endian formats, and vice versa. The swap_.... subroutines are based
on subroutines provided by Jens Hartmann of the Institut fur Geophysik
in Hamburg. These are found in .../cwp/lib/swapbyte.c.
**************************************************************************
Authors: SEP: Einar Kjartansson	CWP: Jack Cohen, Shuki Ronen
swaphval: CWP: John Stockwell
*************************************************************************/
/**************** end self doc ********************************/

#include <string.h>

#include <rsf.h>

#include "hdrpkge.h"
#include "swapbyte.h"

#include "valpkge.h"
#include "segy.h"
/*^*/

static struct {
	const char *key;	const char *type;	int offs;
} hdr[] = {
	{   "tracl",		"i",		0},
	{   "tracr",		"i",		4},
	{    "fldr",		"i",		8},
	{   "tracf",		"i",		12},
	{      "ep",		"i",		16},
	{     "cdp",		"i",		20},
	{    "cdpt",		"i",		24},
	{    "trid",		"h",		28},
	{     "nvs",		"h",		30},
	{     "nhs",		"h",		32},
	{    "duse",		"h",		34},
	{  "offset",		"i",		36},
	{   "gelev",		"i",		40},
	{   "selev",		"i",		44},
	{  "sdepth",		"i",		48},
	{    "gdel",		"i",		52},
	{    "sdel",		"i",		56},
	{   "swdep",		"i",		60},
	{   "gwdep",		"i",		64},
	{  "scalel",		"h",		68},
	{  "scalco",		"h",		70},
	{      "sx",		"i",		72},
	{      "sy",		"i",		76},
	{      "gx",		"i",		80},
	{      "gy",		"i",		84},
	{  "counit",		"h",		88},
	{   "wevel",		"h",		90},
	{  "swevel",		"h",		92},
	{     "sut",		"h",		94},
	{     "gut",		"h",		96},
	{   "sstat",		"h",		98},
	{   "gstat",		"h",		100},
	{   "tstat",		"h",		102},
	{    "laga",		"h",		104},
	{    "lagb",		"h",		106},
	{   "delrt",		"h",		108},
	{    "muts",		"h",		110},
	{    "mute",		"h",		112},
	{      "ns",		"u",		114},
	{      "dt",		"u",		116},
	{    "gain",		"h",		118},
	{     "igc",		"h",		120},
	{     "igi",		"h",		122},
	{    "corr",		"h",		124},
	{     "sfs",		"h",		126},
	{     "sfe",		"h",		128},
	{    "slen",		"h",		130},
	{    "styp",		"h",		132},
	{    "stas",		"h",		134},
	{    "stae",		"h",		136},
	{   "tatyp",		"h",		138},
	{   "afilf",		"h",		140},
	{   "afils",		"h",		142},
	{  "nofilf",		"h",		144},
	{  "nofils",		"h",		146},
	{     "lcf",		"h",		148},
	{     "hcf",		"h",		150},
	{     "lcs",		"h",		152},
	{     "hcs",		"h",		154},
	{    "year",		"h",		156},
	{     "day",		"h",		158},
	{    "hour",		"h",		160},
	{  "minute",		"h",		162},
	{     "sec",		"h",		164},
	{  "timbas",		"h",		166},
	{    "trwf",		"h",		168},
	{  "grnors",		"h",		170},
	{  "grnofr",		"h",		172},
	{  "grnlof",		"h",		174},
	{    "gaps",		"h",		176},
	{   "otrav",		"h",		178},
	{      "d1",		"f",		180},
	{      "f1",		"f",		184},
	{      "d2",		"f",		188},
	{      "f2",		"f",		192},
	{  "ungpow",		"f",		196},
	{ "unscale",		"f",		200},
	{     "ntr",		"i",		204},
	{    "mark",		"h",		208},
	{"shortpad",		"h",		210},
};

static struct {
        const char *key;      const char *type;     int offs;
} bhdr[] = {
            {"jobid",             "i",            0},
            {"lino",              "i",            4},
            {"reno",              "i",            8},
            {"ntrpr",             "h",            12},
            {"nart",              "h",            14},
            {"hdt",               "h",            16},
            {"dto",               "h",            18},
            {"hns",               "h",            20},
            {"nso",               "h",            22},
            {"format",            "h",            24},
            {"fold",              "h",            26},
            {"tsort",             "h",            28},
            {"vscode",            "h",            30},
            {"hsfs",              "h",            32},
            {"hsfe",              "h",            34},
            {"hslen",             "h",            36},
            {"hstyp",             "h",            38},
            {"schn",              "h",            40},
            {"hstas",             "h",            42},
            {"hstae",             "h",            44},
            {"htatyp",            "h",            46},
            {"hcorr",             "h",            48},
            {"bgrcv",             "h",            50},
            {"rcvm",              "h",            52},
            {"mfeet",             "h",            54},
            {"polyt",             "h",            56},
            {"vpol",              "h",            58}
};


void gethval(const segy *tr, int index, Value *valp)
/*< get a trace header word by index >*/
{
	char *tp = (char*) tr;

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(valp->s, tp + hdr[index].offs);  break;
	case 'h': valp->h = *((short*)  (tp + hdr[index].offs)); break;
	case 'u': valp->u = *((unsigned short*) (tp + hdr[index].offs)); break;
	case 'i': valp->i = *((int*)   (tp + hdr[index].offs)); break;
	case 'p': valp->p = *((unsigned int*)   (tp + hdr[index].offs)); break;
	case 'l': valp->l = *((long*)   (tp + hdr[index].offs)); break;
	case 'v': valp->v = *((unsigned long*)  (tp + hdr[index].offs)); break;
	case 'f': valp->f = *((float*)  (tp + hdr[index].offs)); break;
	case 'd': valp->d = *((double*) (tp + hdr[index].offs)); break;
	default: sf_error("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


void puthval(segy *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(tp + hdr[index].offs, valp->s);  break;
	case 'h': *((short*)  (tp + hdr[index].offs)) = valp->h; break;
	case 'u': *((unsigned short*) (tp + hdr[index].offs)) = valp->u; break;
	case 'i': *((int*)   (tp + hdr[index].offs)) = valp->i; break;
	case 'p': *((unsigned int*)   (tp + hdr[index].offs)) = valp->p; break;
	case 'l': *((long*)   (tp + hdr[index].offs)) = valp->l; break;
	case 'v': *((unsigned long*)  (tp + hdr[index].offs)) = valp->v; break;
	case 'f': *((float*)  (tp + hdr[index].offs)) = valp->f; break;
	case 'd': *((double*) (tp + hdr[index].offs)) = valp->d; break;
	default: sf_error("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


void getbhval(const bhed *bh, int index, Value *valp)
{
	char *bhp = (char*) bh;

	switch(*(bhdr[index].type)) {
	case 'h': valp->h = *((short*) (bhp + bhdr[index].offs)); break;
	case 'i': valp->i = *((int*)   (bhp + bhdr[index].offs)); break;
	default: sf_error("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void putbhval(bhed *bh, int index, Value *valp)
{
	char *bhp = (char*) bh;

	switch(*(bhdr[index].type)) {
	case 'h': *((short*) (bhp + bhdr[index].offs)) = valp->h; break;
	case 'i': *((int*)   (bhp + bhdr[index].offs)) = valp->i; break;
	default: sf_error("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void gethdval(const segy *tr, const char *key, Value *valp)
/*< get a trace header word by name >*/
{
	int index = getindex(key);
	char *tp = (char*) tr;

	if ( -1 == (index))
		sf_error("%s: key word not in segy.h: '%s'", __FILE__, key);

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(valp->s, tp + hdr[index].offs);  break;
	case 'h': valp->h = *((short*)  (tp + hdr[index].offs)); break;
	case 'u': valp->u = *((unsigned short*) (tp + hdr[index].offs)); break;
	case 'i': valp->i = *((int*)   (tp + hdr[index].offs)); break;
	case 'p': valp->p = *((unsigned int*)   (tp + hdr[index].offs)); break;
	case 'l': valp->l = *((long*)   (tp + hdr[index].offs)); break;
	case 'v': valp->v = *((unsigned long*)  (tp + hdr[index].offs)); break;
	case 'f': valp->f = *((float*)  (tp + hdr[index].offs)); break;
	case 'd': valp->d = *((double*) (tp + hdr[index].offs)); break;
	default: sf_error("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


void puthdval(segy *tr, const char *key, Value *valp)
/*< put a trace header word by name >*/
{
	int index = getindex(key);
	char *tp = (char*) tr;

	if ( -1 == (index))
		sf_error("%s: key word not in segy.h: '%s'", __FILE__, key);

	switch(*(hdr[index].type)) {
	case 's': (void) strcpy(tp + hdr[index].offs, valp->s);  break;
	case 'h': *((short*)  (tp + hdr[index].offs)) = valp->h; break;
	case 'u': *((unsigned short*) (tp + hdr[index].offs)) = valp->u; break;
	case 'i': *((int*)   (tp + hdr[index].offs)) = valp->i; break;
	case 'p': *((unsigned int*)   (tp + hdr[index].offs)) = valp->p; break;
	case 'l': *((long*)   (tp + hdr[index].offs)) = valp->l; break;
	case 'v': *((unsigned long*)  (tp + hdr[index].offs)) = valp->v; break;
	case 'f': *((float*)  (tp + hdr[index].offs)) = valp->f; break;
	case 'd': *((double*) (tp + hdr[index].offs)) = valp->d; break;
	default: sf_error("%s: %s: mysterious data type", __FILE__,__LINE__); break;
	}

	return;
}


const char *hdtype(const char *key)
/*< get the data type of a trace header word by name >*/
{
	int index = getindex(key);

	if (-1 == (index))
		sf_error("%s: key word not in segy.h: '%s'", __FILE__, key);

	return hdr[index].type;
}

void printftype(register cwp_String key)
{
	switch(*hdtype(key)) {
	case 's':
		(void) printf("char\n");
	break;
	case 'h':
		(void) printf("short\n");
	break;
	case 'u':
		(void) printf("ushort\n");
	break;
	case 'i':
		(void) printf("int\n");
	break;
	case 'p':
		(void) printf("unsigned int\n");
	break;
	case 'l':
		(void) printf("long\n");
	break;
	case 'v':
		(void) printf("ulong\n");
	break;
	case 'f':
		(void) printf("float\n");
	break;
	case 'd':
		(void) printf("double\n");
	break;
	case 'U':
		(void) printf("unsigned int:16\n");
	break;
	case 'P':
		(void) printf("unsigned int:32\n");
	break;
	default:
		sf_error("printftype: unknown type %s", hdtype(key));
	break;
	}
	return;
}

#define SU_NKEYS	80	/* Number of key header words */
#define HDRBYTES	240	/* Bytes in the trace header */
#define	MAXSEGY		262380
/*^*/

#define   STREQ(s,t) (strcmp(s,t) == 0)

const char *getkey(const int index)
{
	return (index < SU_NKEYS && index >= 0) ? hdr[index].key : NULL;
}


int getindex(const char *key)	
/*< get index for this key >*/
{
	register int i;

	for (i = 0; i < SU_NKEYS; i++)
		if (STREQ(hdr[i].key, key))
			return i;	/* key found */

	/* not found */
	return -1;
}

void swaphval(segy *tr, int index)
{
	register char *tp= (char *) tr;

        switch(*(hdr[index].type)) {
        case 'h': swap_short_2((short*)(tp + hdr[index].offs));
	break;
        case 'u': swap_u_short_2((unsigned short*)(tp + hdr[index].offs));
	break;
        case 'i': swap_int_4((int*)(tp + hdr[index].offs));
	break;
        case 'p': swap_u_int_4((unsigned int*)(tp + hdr[index].offs));
	break;
        case 'l': swap_long_4((long*)(tp + hdr[index].offs));
	break;
        case 'v': swap_u_long_4((unsigned long*)(tp + hdr[index].offs));
	break;
        case 'f': swap_float_4((float*)(tp + hdr[index].offs));
	break;
        case 'd': swap_double_8((double*)(tp + hdr[index].offs));
	break;
        default: sf_error("%s: %s: unsupported data type", __FILE__, __LINE__);
	break;
        }

        return;
}


void swapbhval(bhed *bh, int index)
{
	register char *bhp= (char *) bh;

        switch(*(bhdr[index].type)) {
        case 'h': swap_short_2((short*)(bhp + bhdr[index].offs)); break;
        case 'i': swap_int_4((int*)(bhp + bhdr[index].offs)); break;
        default: sf_error("%s: %s: unsupported data type", __FILE__, __LINE__);
	break;
        }

        return;
}

typedef struct {        /* tapesegy - trace identification header */

        unsigned int tracl:32;  /* trace sequence number within line */

        unsigned int tracr:32;  /* trace sequence number within reel */

        unsigned int fldr:32;   /* field record number */

        unsigned int tracf:32;  /* trace number within field record */

        unsigned int ep:32;     /* energy source point number */

        unsigned int cdp:32;    /* CDP ensemble number */

        unsigned int cdpt:32;   /* trace number within CDP ensemble */

        unsigned int trid:16;   /* trace identification code:
                        1 = seismic data
                        2 = dead
                        3 = dummy
                        4 = time break
                        5 = uphole
                        6 = sweep
                        7 = timing
                        8 = water break
                        9---, N = optional use (N = 32,767) */

        unsigned int nvs:16;    /* number of vertically summed traces (see
                        vscode in bhed structure) */

        unsigned int nhs:16;    /* number of horizontally summed traces (see
                        vscode in bhed structure) */

        unsigned int duse:16;   /* data use:
                                1 = production
                                2 = test */

        unsigned int offset:32; /* distance from source point to receiver
                           group (negative if opposite to direction
                           in which the line was shot) */

        unsigned int gelev:32;  /* receiver group elevation from sea level
                           (above sea level is positive) */

        unsigned int selev:32;  /* source elevation from sea level
                           (above sea level is positive) */

        unsigned int sdepth:32; /* source depth (positive) */

        unsigned int gdel:32;   /* datum elevation at receiver group */

        unsigned int sdel:32;   /* datum elevation at source */

        unsigned int swdep:32;  /* water depth at source */

        unsigned int gwdep:32;  /* water depth at receiver group */

        unsigned int scalel:16; /* scale factor for previous 7 entries
                           with value plus or minus 10 to the
                           power 0, 1, 2, 3, or 4 (if positive,
                           multiply, if negative divide) */

        unsigned int scalco:16; /* scale factor for next 4 entries
                           with value plus or minus 10 to the
                           power 0, 1, 2, 3, or 4 (if positive,
                           multiply, if negative divide) */

        unsigned int  sx:32;    /* X source coordinate */

        unsigned int  sy:32;    /* Y source coordinate */

        unsigned int  gx:32;    /* X group coordinate */

        unsigned int  gy:32;    /* Y source coordinate */

        unsigned int counit:16; /* coordinate units code:
                                for previoius four entries
                                1 = length (meters or feet)
                                2 = seconds of arc (in this case, the
                                X values are unsigned intitude and the Y values
                                are latitude, a positive value designates
                                the number of seconds east of Greenwich
                                or north of the equator */

        unsigned int wevel:16;  /* weathering velocity */

        unsigned int swevel:16; /* subweathering velocity */

        unsigned int sut:16;    /* uphole time at source */

        unsigned int gut:16;    /* uphole time at receiver group */

        unsigned int sstat:16;  /* source static correction */

        unsigned int gstat:16;  /* group static correction */

        unsigned int tstat:16;  /* total static applied */

        unsigned int laga:16;   /* lag time A, time in ms between end of 240-
                           byte trace identification header and time
                           break, positive if time break occurs after
                           end of header, time break is defined as
                           the initiation pulse which maybe recorded
                           on an auxiliary trace or as otherwise
                           specified by the recording system */

        unsigned int lagb:16;   /* lag time B, time in ms between the time
                           break and the initiation time of the energy source,
                           may be positive or negative */

        unsigned int delrt:16;  /* delay recording time, time in ms between
                           initiation time of energy source and time
                           when recording of data samples begins
                           (for deep water work if recording does not
                           start at zero time) */

        unsigned int muts:16;   /* mute time--start */

        unsigned int mute:16;   /* mute time--end */

        unsigned int ns:16;     /* number of samples in this trace */

        unsigned int dt:16;     /* sample interval; in micro-seconds */

        unsigned int gain:16;   /* gain type of field instruments code:
                                1 = fixed
                                2 = binary
                                3 = floating point
                                4 ---- N = optional use */

        unsigned int igc:16;    /* instrument gain constant */

        unsigned int igi:16;    /* instrument early or initial gain */

        unsigned int corr:16;   /* correlated:
                                1 = no
                                2 = yes */

        unsigned int sfs:16;    /* sweep frequency at start */

        unsigned int sfe:16;    /* sweep frequency at end */

        unsigned int slen:16;   /* sweep length in ms */

        unsigned int styp:16;   /* sweep type code:
                                1 = linear
                                2 = cos-squared
                                3 = other */

        unsigned int stas:16;   /* sweep trace length at start in ms */

        unsigned int stae:16;   /* sweep trace length at end in ms */

        unsigned int tatyp:16;  /* taper type: 1=linear, 2=cos^2, 3=other */

        unsigned int afilf:16;  /* alias filter frequency if used */

        unsigned int afils:16;  /* alias filter slope */

        unsigned int nofilf:16; /* notch filter frequency if used */

        unsigned int nofils:16; /* notch filter slope */

        unsigned int lcf:16;    /* low cut frequency if used */

        unsigned int hcf:16;    /* high cut frequncy if used */

        unsigned int lcs:16;    /* low cut slope */

        unsigned int hcs:16;    /* high cut slope */

        unsigned int year:16;   /* year data recorded */

        unsigned int day:16;    /* day of year */

        unsigned int hour:16;   /* hour of day (24 hour clock) */

        unsigned int minute:16; /* minute of hour */

        unsigned int sec:16;    /* second of minute */

        unsigned int timbas:16; /* time basis code:
                                1 = local
                                2 = GMT
                                3 = other */

        unsigned int trwf:16;   /* trace weighting factor, defined as 1/2^N
                           volts for the least sigificant bit */

        unsigned int grnors:16; /* geophone group number of roll switch
                           position one */

        unsigned int grnofr:16; /* geophone group number of trace one within
                           original field record */

        unsigned int grnlof:16; /* geophone group number of last trace within
                           original field record */

        unsigned int gaps:16;   /* gap size (total number of groups dropped) */

        unsigned int otrav:16;  /* overtravel taper code:
                                1 = down (or behind)
                                2 = up (or ahead) */

        unsigned char unass[60];        /* unassigned */ 

        /* pseudo_float data[SU_NFLTS]; */
        
        unsigned char data[SU_NFLTS][4];

} tapesegy;

/* tapehdr.h - include file for SEGY traces as bytes (only for segyread,write)
 *
 * Reference:
 *      K. M. Barry, D. A. Cavers and C. W. Kneale, "Special Report:
 *              Recommended Standards for Digital Tape Formats",
 *              Geophysics, vol. 40, no. 2 (April 1975), P. 344-352.
 *      
 * $Author: john $
 * $Source: /usr/local/cwp/src/su/include/RCS/tapehdr.h,v $
 * $Revision: 1.1 $ ; $Date: 1996/09/09 16:18:41 $
 */ 

static struct {
        const char *key;      const char *type;     int offs;
} tapehdr[] = {
           {"tracl",             "P",            0},
           {"tracr",             "P",            4},
            {"fldr",             "P",            8},
           {"tracf",             "P",            12},
              {"ep",             "P",            16},
             {"cdp",             "P",            20},
            {"cdpt",             "P",            24},
            {"trid",             "U",            28},
             {"nvs",             "U",            30},
             {"nhs",             "U",            32},
            {"duse",             "U",            34},
          {"offset",             "P",            36},
           {"gelev",             "P",            40},
           {"selev",             "P",            44},
          {"sdepth",             "P",            48},
            {"gdel",             "P",            52},
            {"sdel",             "P",            56},
           {"swdep",             "P",            60},
           {"gwdep",             "P",            64},
          {"scalel",             "U",            68},
          {"scalco",             "U",            70},
              {"sx",             "P",            72},
              {"sy",             "P",            76},
              {"gx",             "P",            80},
              {"gy",             "P",            84},
          {"counit",             "U",            88},
           {"wevel",             "U",            90},
          {"swevel",             "U",            92},
             {"sut",             "U",            94},
             {"gut",             "U",            96},
           {"sstat",             "U",            98},
           {"gstat",             "U",            100},
           {"tstat",             "U",            102},
            {"laga",             "U",            104},
            {"lagb",             "U",            106},
           {"delrt",             "U",            108},
            {"muts",             "U",            110},
            {"mute",             "U",            112},
              {"ns",             "U",            114},
              {"dt",             "U",            116},
            {"gain",             "U",            118},
             {"igc",             "U",            120},
             {"igi",             "U",            122},
            {"corr",             "U",            124},
             {"sfs",             "U",            126},
             {"sfe",             "U",            128},
            {"slen",             "U",            130},
            {"styp",             "U",            132},
            {"stas",             "U",            134},
            {"stae",             "U",            136},
           {"tatyp",             "U",            138},
           {"afilf",             "U",            140},
           {"afils",             "U",            142},
          {"nofilf",             "U",            144},
          {"nofils",             "U",            146},
             {"lcf",             "U",            148},
             {"hcf",             "U",            150},
             {"lcs",             "U",            152},
             {"hcs",             "U",            154},
            {"year",             "U",            156},
             {"day",             "U",            158},
            {"hour",             "U",            160},
          {"minute",             "U",            162},
             {"sec",             "U",            164},
          {"timbas",             "U",            166},
            {"trwf",             "U",            168},
          {"grnors",             "U",            170},
          {"grnofr",             "U",            172},
          {"grnlof",             "U",            174},
            {"gaps",             "U",            176},
           {"otrav",             "U",            178},
};

static struct {
        const char *key;      const char *type;     int offs;
} tapebhdr[] = {
           {"jobid",             "P",            0},
           {"lino",              "P",            4},
           {"reno",              "P",            8},
           {"ntrpr",             "U",            12},
           {"nart",              "U",            14},
           {"hdt",               "U",            16},
           {"dto",               "U",            18},
           {"hns",               "U",            20},
           {"nso",               "U",            22},
           {"format",            "U",            24},
           {"fold",              "U",            26},
           {"tsort",             "U",            28},
           {"vscode",            "U",            30},
           {"hsfs",              "U",            32},
           {"hsfe",              "U",            34},
           {"hslen",             "U",            36},
           {"hstyp",             "U",            38},
           {"schn",              "U",            40},
           {"hstas",             "U",            42},
           {"hstae",             "U",            44},
           {"htatyp",            "U",            46},
           {"hcorr",             "U",            48},
           {"bgrcv",             "U",            50},
           {"rcvm",              "U",            52},
           {"mfeet",             "U",            54},
           {"polyt",             "U",            56},
           {"vpol",              "U",            58},
};


void gettapehval(const tapesegy *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(tapehdr[index].type)) {
	case 'U': valp->h = (short) *((short*) (tp + tapehdr[index].offs));
	break;
	case 'P': valp->i = (int) *((int*) (tp + tapehdr[index].offs));
	break;
	default: sf_error("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void puttapehval(tapesegy *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(tapehdr[index].type)) {
	case 'U': *((short*) (tp + tapehdr[index].offs)) = valp->h; break;
	case 'P': *((int*)   (tp + tapehdr[index].offs)) = valp->i; break;
	default: sf_error("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}

/* tapebhdr.h - include file for SEGY-defined binary header 
 *
 * Reference:
 *      K. M. Barry, D. A. Cavers and C. W. Kneale, "Special Report:
 *              Recommended Standards for Digital Tape Formats",
 *              Geophysics, vol. 40, no. 2 (April 1975), P. 344-352.
 *      
 * $Author: john $
 * $Source: /usr/local/cwp/src/su/include/RCS/tapebhdr.h,v $
 * $Revision: 1.2 $ ; $Date: 1996/09/09 17:11:01 $
 */ 

typedef struct {        /* bhedtape - binary header */

        unsigned int jobid:32;  /* job identification number */

        unsigned int lino:32;   /* line number (only one line per reel) */

        unsigned int reno:32;   /* reel number */

        unsigned int ntrpr:16;  /* number of data traces per record */

        unsigned int nart:16;   /* number of auxiliary traces per record */

        unsigned int hdt:16;    /* sample interval (microsecs) for this reel */

        unsigned int dto:16;    /* same for original field recording */

        unsigned int hns:16;    /* number of samples per trace for this reel */

        unsigned int nso:16;    /* same for original field recording */

        unsigned int format:16; /* data sample format code:
                                1 = floating point (4 bytes)
                                2 = fixed point (4 bytes)
                                3 = fixed point (2 bytes)
                                4 = fixed point w/gain code (4 bytes) */

        unsigned int fold:16;   /* CDP fold expected per CDP ensemble */

        unsigned int tsort:16;  /* trace sorting code: 
                                1 = as recorded (no sorting)
                                2 = CDP ensemble
                                3 = single fold continuous profile
                                4 = horizontally stacked */

        unsigned int vscode:16; /* vertical sum code:
                                1 = no sum
                                2 = two sum ...
                                N = N sum (N = 32,767) */

        unsigned int hsfs:16;   /* sweep frequency at start */

        unsigned int hsfe:16;   /* sweep frequency at end */

        unsigned int hslen:16;  /* sweep length (ms) */

        unsigned int hstyp:16;  /* sweep type code:
                                1 = linear
                                2 = parabolic
                                3 = exponential
                                4 = other */

        unsigned int schn:16;   /* trace number of sweep channel */

        unsigned int hstas:16;  /* sweep trace taper length at start if
                           tapered (the taper starts at zero time
                           and is effective for this length) */

        unsigned int hstae:16;  /* sweep trace taper length at end (the ending
                           taper starts at sweep length minus the taper
                           length at end) */

        unsigned int htatyp:16; /* sweep trace taper type code:
                                1 = linear
                                2 = cos-squared
                                3 = other */

        unsigned int hcorr:16;  /* correlated data traces code:
                                1 = no
                                2 = yes */

        unsigned int bgrcv:16;  /* binary gain recovered code:
                                1 = yes
                                2 = no */

        unsigned int rcvm:16;   /* amplitude recovery method code:
                                1 = none
                                2 = spherical divergence
                                3 = AGC
                                4 = other */

        unsigned int mfeet:16;  /* measurement system code:
                                1 = meters
                                2 = feet */

        unsigned int polyt:16;  /* impulse signal polarity code:
                                1 = increase in pressure or upward
                                    geophone case movement gives
                                    negative number on tape
                                2 = increase in pressure or upward
                                    geophone case movement gives
                                    positive number on tape */

        unsigned int vpol:16;   /* vibratory polarity code:
                                code    seismic signal lags pilot by
                                1       337.5 to  22.5 degrees
                                2        22.5 to  67.5 degrees
                                3        67.5 to 112.5 degrees
                                4       112.5 to 157.5 degrees
                                5       157.5 to 202.5 degrees
                                6       202.5 to 247.5 degrees
                                7       247.5 to 292.5 degrees
                                8       293.5 to 337.5 degrees */

        unsigned char hunass[340];      /* unassigned */

} tapebhed;

void gettapebhval(const tapebhed *tr, int index, Value *valp)
{
	char *tp = (char*) tr;

	switch(*(tapebhdr[index].type)) {
	case 'U': valp->h = (short) *((short*) (tp + tapebhdr[index].offs));
	break;
	case 'P': valp->i = (int) *((int*) (tp + tapebhdr[index].offs));
	break;
	default: sf_error("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}


void puttapebhval(tapebhed *bh, int index, Value *valp)
{
	char *bhp = (char*) bh;

	switch(*(tapebhdr[index].type)) {
	case 'U': *((short*) (bhp + tapebhdr[index].offs)) = valp->h; break;
	case 'P': *((int*)   (bhp + tapebhdr[index].offs)) = valp->i; break;
	default: sf_error("%s: %s: mysterious data type", __FILE__, __LINE__);
	break;
	}

	return;
}

/* Display non-null header field values */
void printheader(const segy *tp)
{
	int i;			/* index over header fields		*/
	int j;			/* index over non-null header fields	*/
	Value val;		/* value in header field		*/
	const char* type;	/* ... its data type			*/
	const char* key;	/* ... the name of the header field	*/
	Value zeroval;		/* zero value to compare with		*/

	zeroval.l = 0;
	j = 0;
	for (i = 0; i < SU_NKEYS; i++) {
		gethval(tp, i, &val);
		key = getkey(i);
		type = hdtype(key);
		if (valcmp(type, val, zeroval)) { /* not equal to zero */
			(void) printf(" %s=", key);
			printfval(type, val);
			if ((++j % 6) == 0) putchar('\n');
		}
	}
	putchar('\n');

	return;
}
