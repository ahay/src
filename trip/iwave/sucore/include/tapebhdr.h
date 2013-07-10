/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

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

#ifndef TAPEBHDR_H
#define TAPEBHDR_H

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

/* FUNCTION PROTOTYPES */
#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

void gettapebhval(const tapebhed *tapetr, int index, Value *valp);
void puttapebhval(tapebhed *tapetr, int index, Value *valp);

#ifdef __cplusplus /* if C++, end external linkage specification */
}
#endif

#endif	/* end  TAPEBHDR_H */
