#ifndef CJBSEGY_H
#define CJBSEGY_H


/* TYPEDEFS */
typedef struct {	/* segy - trace identification header */

        int tracl;      /* 1-4, trace sequence number within line */

        int tracr;      /* 5-8, trace sequence number within reel */

        int fldr;       /* 9-12, field record number */

	int tracf;	/* 13-16, trace number within field record */
       
	int ep;	        /* 17-20, energy source point number */

	int cdp;	/* 21-24, CDP ensemble number */

	int cdpt;	/* 25-28, trace number within CDP ensemble */

	short trid;	/* 29-30, trace identification code:
			1 = seismic data
			2 = dead
			3 = dummy
			4 = time break
			5 = uphole
			6 = sweep
			7 = timing
			8 = water break
			9---, N = optional use (N = 32,767)

			Following are CWP id flags:

			 9 = autocorrelation

			10 = Fourier transformed - no packing
			     xr[0],xi[0], ..., xr[N-1],xi[N-1]

			11 = Fourier transformed - unpacked Nyquist
			     xr[0],xi[0],...,xr[N/2],xi[N/2]

			12 = Fourier transformed - packed Nyquist
	 		     even N:
			     xr[0],xr[N/2],xr[1],xi[1], ...,
				xr[N/2 -1],xi[N/2 -1]
				(note the exceptional second entry)
			     odd N:
			     xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
				xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
				(note the exceptional second & last entries)

			13 = Complex signal in the time domain
			     xr[0],xi[0], ..., xr[N-1],xi[N-1]

			14 = Fourier transformed - amplitude/phase
			     a[0],p[0], ..., a[N-1],p[N-1]

			15 = Complex time signal - amplitude/phase
			     a[0],p[0], ..., a[N-1],p[N-1]

			16 = Real part of complex trace from 0 to Nyquist

			17 = Imag part of complex trace from 0 to Nyquist

			18 = Amplitude of complex trace from 0 to Nyquist

			19 = Phase of complex trace from 0 to Nyquist

			21 = Wavenumber time domain (k-t)

			22 = Wavenumber frequency (k-omega)

			23 = Envelope of the complex time trace

			24 = Phase of the complex time trace

			25 = Frequency of the complex time trace

			30 = Depth-Range (z-x) traces

			43 = Seismic Data, Vertical Component 

			44 = Seismic Data, Horizontal Component 1 

			45 = Seismic Data, Horizontal Component 2 

			46 = Seismic Data, Radial Component

			47 = Seismic Data, Transverse Component  

			101 = Seismic data packed to bytes (by supack1)
			
			102 = Seismic data packed to 2 bytes (by supack2)
			*/

	short nvs;	/* 31-32, number of vertically summed traces (see vscode
			   in bhed structure) */

	short nhs;	/* 33-34, number of horizontally summed traces (see vscode
			   in bhed structure) */

	short duse;	/* 35-36, data use:
				1 = production
				2 = test */

	int offset;	/* 37-40, distance from source point to receiver
			   group (negative if opposite to direction
			   in which the line was shot) */

	int gelev;	/* 41-44, receiver group elevation from sea level
			   (above sea level is positive) */

	int selev;	/* 45-48, source elevation from sea level
			   (above sea level is positive) */

	int sdepth;	/* 49-52, source depth (positive) */

	int gdel;	/* 53-56, datum elevation at receiver group */

	int sdel;	/* 57-60, datum elevation at source */

	int swdep;	/* 61-64, water depth at source */

	int gwdep;	/* 65-68, water depth at receiver group */

	short scalel;	/* 69-70, scale factor for previous 7 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */

	short scalco;	/* 71-72, scale factor for next 4 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */

	int  sx;	/* 73-76, X source coordinate */

	int  sy;	/* 77-80, Y source coordinate */

	int  gx;	/* 81-84, X group coordinate */

	int  gy;	/* 85-88, Y group coordinate */

	short counit;	/* 89-90, coordinate units code:
				for previous four entries
				1 = length (meters or feet)
				2 = seconds of arc (in this case, the
				X values are longitude and the Y values
				are latitude, a positive value designates
				the number of seconds east of Greenwich
				or north of the equator */

	short wevel;	/* 91-92, weathering velocity */

	short swevel;	/* 93-94, subweathering velocity */

	short sut;	/* 95-96, uphole time at source */

	short gut;	/* 97-98, uphole time at receiver group */

	short sstat;	/* 99-100, source static correction */

	short gstat;	/* 101-102, group static correction */

	short tstat;	/* 103-104, total static applied */

	short laga;	/* 105-106, lag time A, time in ms between end of 240-
			   byte trace identification header and time
			   break, positive if time break occurs after
			   end of header, time break is defined as
			   the initiation pulse which maybe recorded
			   on an auxiliary trace or as otherwise
			   specified by the recording system */

	short lagb;	/* 107-108, lag time B, time in ms between the time break
			   and the initiation time of the energy source,
			   may be positive or negative */

	short delrt;	/* 109-110, delay recording time, time in ms between
			   initiation time of energy source and time
			   when recording of data samples begins
			   (for deep water work if recording does not
			   start at zero time) */

	short muts;	/* 111-112, mute time--start */

	short mute;	/* 113-114, mute time--end */

	unsigned short ns;	/* 115-116, number of samples in this trace */

	unsigned short dt;	/* 117-118, sample interval; in micro-seconds */

	short gain;	/* 119-120, gain type of field instruments code:
				1 = fixed
				2 = binary
				3 = floating point
				4 ---- N = optional use */

	short igc;	/* 121-122, instrument gain constant */

	short igi;	/* 123-124, instrument early or initial gain */

	short corr;	/* 125-126, correlated:
				1 = no
				2 = yes */

	short sfs;	/* 127-128, sweep frequency at start */

	short sfe;	/* 129-130, sweep frequency at end */

	short slen;	/* 131-132, sweep length in ms */

	short styp;	/* 133-134, sweep type code:
				1 = linear
				2 = cos-squared
				3 = other */

	short stas;	/* 135-136, sweep trace length at start in ms */

	short stae;	/* 137-138, sweep trace length at end in ms */

	short tatyp;	/* 139-140, taper type: 1=linear, 2=cos^2, 3=other */

	short afilf;	/* 141-142, alias filter frequency if used */

	short afils;	/* 143-144, alias filter slope */

	short nofilf;	/* 145-146, notch filter frequency if used */

	short nofils;	/* 147-148, notch filter slope */

	short lcf;	/* 149-150, low cut frequency if used */

	short hcf;	/* 151-152, high cut frequncy if used */

	short lcs;	/* 153-154, low cut slope */

	short hcs;	/* 155-156, high cut slope */

	short year;	/* 157-158, year data recorded */

	short day;	/* 159-160, day of year */

	short hour;	/* 161-162, hour of day (24 hour clock) */

	short minute;	/* 163-164, minute of hour */

	short sec;	/* 165-166, second of minute */

	short timbas;	/* 167-168, time basis code:
				1 = local
				2 = GMT
				3 = other */

	short trwf;	/* 169-170, trace weighting factor, defined as 1/2^N
			   volts for the least sigificant bit */

	short grnors;	/* 171-172, geophone group number of roll switch
			   position one */

	short grnofr;	/* 173-174, geophone group number of trace one within
			   original field record */

	short grnlof;	/* 175-176, geophone group number of last trace within
			   original field record */

	short gaps;	/* 177-178, gap size (total number of groups dropped) */

	short otrav;	/* 179-180, overtravel taper code:
				1 = down (or behind)
				2 = up (or ahead) */

	/* local assignments */
	float d1;	/* 181-184, sample spacing for non-seismic data */

	float f1;	/* 185-188, first sample location for non-seismic data */

	float d2;	/* 189-192, sample spacing between traces */

	float f2;	/* 193-196, first trace location */

	float ungpow;	/* 197-200, negative of power used for dynamic
			   range compression */

	float unscale;	/* 201-204, reciprocal of scaling factor to normalize
			   range */

	int ntr; 	/* 205-208, number of traces */

	short mark;	/* 209-210, mark selected traces */

        short shortpad; /* 211-212, alignment padding */


	int rn;         /* 213-216, the receiver number added by lyz in SurfaceVel system*/

	int nstations;  /* 217-220, the receiver number added by lyz in SurfaceVel system*/
	int nstationg;  /* 221-224, the receiver number added by lyz in SurfaceVel system*/
	short unass[8];	/* unassigned--NOTE: last entry causes 
			   a break in the word alignment, if we REALLY
			   want to maintain 240 bytes, the following
			   entry should be an odd number of short/UINT2
			   OR do the insertion above the "mark" keyword
			   entry */
} cjbsegy;


typedef struct {	/* bhed - binary header */

	int jobid;	/* job identification number */

	int lino;	/* line number (only one line per reel) */

	int reno;	/* reel number */

	short ntrpr;	/* number of data traces per record */

        short nart;	/* number of auxiliary traces per record */

	unsigned short hdt; /* sample interval in micro secs for this reel */

	unsigned short dto; /* same for original field recording */

	unsigned short hns; /* number of samples per trace for this reel */

	unsigned short nso; /* same for original field recording */

	short format;	/* data sample format code:
				1 = floating point (4 bytes)
				2 = fixed point (4 bytes)
				3 = fixed point (2 bytes)
				4 = fixed point w/gain code (4 bytes) */

	short fold;	/* CDP fold expected per CDP ensemble */

	short tsort;	/* trace sorting code: 
				1 = as recorded (no sorting)
				2 = CDP ensemble
				3 = single fold continuous profile
				4 = horizontally stacked */

	short vscode;	/* vertical sum code:
				1 = no sum
				2 = two sum ...
				N = N sum (N = 32,767) */

	short hsfs;	/* sweep frequency at start */

	short hsfe;	/* sweep frequency at end */

	short hslen;	/* sweep length (ms) */

	short hstyp;	/* sweep type code:
				1 = linear
				2 = parabolic
				3 = exponential
				4 = other */

	short schn;	/* trace number of sweep channel */

	short hstas;	/* sweep trace taper length at start if
			   tapered (the taper starts at zero time
			   and is effective for this length) */

	short hstae;	/* sweep trace taper length at end (the ending
			   taper starts at sweep length minus the taper
			   length at end) */

	short htatyp;	/* sweep trace taper type code:
				1 = linear
				2 = cos-squared
				3 = other */

	short hcorr;	/* correlated data traces code:
				1 = no
				2 = yes */

	short bgrcv;	/* binary gain recovered code:
				1 = yes
				2 = no */

	short rcvm;	/* amplitude recovery method code:
				1 = none
				2 = spherical divergence
				3 = AGC
				4 = other */

	short mfeet;	/* measurement system code:
				1 = meters
				2 = feet */

	short polyt;	/* impulse signal polarity code:
				1 = increase in pressure or upward
				    geophone case movement gives
				    negative number on tape
				2 = increase in pressure or upward
				    geophone case movement gives
				    positive number on tape */

	short vpol;	/* vibratory polarity code:
				code	seismic signal lags pilot by
				1	337.5 to  22.5 degrees
				2	 22.5 to  67.5 degrees
				3	 67.5 to 112.5 degrees
				4	112.5 to 157.5 degrees
				5	157.5 to 202.5 degrees
				6	202.5 to 247.5 degrees
				7	247.5 to 292.5 degrees
				8	293.5 to 337.5 degrees */

	short hunass[170];	/* unassigned */

} bhed;

typedef struct {	/* 3200 */
char y3200[3200];
} Y_3200;

#endif
