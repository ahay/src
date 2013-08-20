/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/* tapesegy.h - include file for SEGY traces as bytes (only for segyread,write)
 *
 * declarations for:
 *      typedef struct {} segytape - the trace identification header
 *      typedef struct {} bhedtape - binary header
 *      typedef static struct {} tapehdr -  type and offsets for trace header
 *      typedef static struct {} tapebhdr - type and offsets for binary header
 *
 * Note: last two commented off for now.
 *
 * Reference:
 *      K. M. Barry, D. A. Cavers and C. W. Kneale, "Special Report:
 *              Recommended Standards for Digital Tape Formats",
 *              Geophysics, vol. 40, no. 2 (April 1975), P. 344-352.
 *      
 * $Author: john $
 * $Source: /usr/local/cwp/src/su/include/RCS/tapesegy.h,v $
 * $Revision: 1.7 $ ; $Date: 1996/09/09 17:11:38 $
 */ 

#ifndef TAPESEGY_H
#define TAPESEGY_H

/* TYPEDEFS */
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

/* FUNCTION PROTOTYPES */
#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

void gettapehval(const tapesegy *tapetr, int index, Value *valp);
void puttapehval(tapesegy *tapetr, int index, Value *valp);

#ifdef __cplusplus /* if C++, end external linkage specification */
}
#endif

#endif
