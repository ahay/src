/* Convert a SEG-Y or SU dataset to RSF.

Data headers and trace headers are separated from the data.

"suread" is equivalent to "segyread su=y"


SEGY key names:

tracl: trace sequence number within line 0

tracr: trace sequence number within reel 4

fldr:     field record number 8 

tracf:    trace number within field record 12 

ep:       energy source point number 16 

cdp:      CDP ensemble number 20 

cdpt:     trace number within CDP ensemble 24 

trid:     trace identification code:
1 = seismic data
2 = dead
3 = dummy
4 = time break
5 = uphole
6 = sweep
7 = timing
8 = water break
9---, N = optional use (N = 32,767) 28 

nvs:      number of vertically summed traces 30 

nhs:      number of horizontally summed traces 32 

duse:     data use:
1 = production
2 = test 34

offset:   distance from source point to receiver
group (negative if opposite to direction
in which the line was shot) 36 

gelev:    receiver group elevation from sea level
(above sea level is positive) 40 

selev:    source elevation from sea level
(above sea level is positive) 44 

sdepth:   source depth (positive) 48 

gdel:     datum elevation at receiver group 52 

sdel:     datum elevation at source 56 

swdep:    water depth at source 60 

gwdep:    water depth at receiver group 64 

scalel:   scale factor for previous 7 entries
with value plus or minus 10 to the
power 0, 1, 2, 3, or 4 (if positive,
multiply, if negative divide) 68 

scalco:   scale factor for next 4 entries
with value plus or minus 10 to the
power 0, 1, 2, 3, or 4 (if positive,
multiply, if negative divide) 70 

sx:       X source coordinate 72 

sy:       Y source coordinate 76 

gx:       X group coordinate 80 

gy:       Y group coordinate 84 

counit:   coordinate units code:
for previous four entries
1 = length (meters or feet)
2 = seconds of arc (in this case, the
X values are unsigned longitude and the Y values
are latitude, a positive value designates
the number of seconds east of Greenwich
or north of the equator 88 

wevel:     weathering velocity 90 

swevel:    subweathering velocity 92 

sut:       uphole time at source 94 

gut:       uphole time at receiver group 96 

sstat:     source static correction 98 

gstat:     group static correction 100 

tstat:     total static applied 102 

laga:      lag time A, time in ms between end of 240-
byte trace identification header and time
break, positive if time break occurs after
end of header, time break is defined as
the initiation pulse which maybe recorded
on an auxiliary trace or as otherwise
specified by the recording system 104 

lagb:      lag time B, time in ms between the time
break and the initiation time of the energy source,
may be positive or negative 106 

delrt:     delay recording time, time in ms between
initiation time of energy source and time
when recording of data samples begins
(for deep water work if recording does not
start at zero time) 108 

muts:      mute time--start 110 

mute:      mute time--end 112 

ns:        number of samples in this trace 114 

dt:        sample interval, in micro-seconds 116 

gain:      gain type of field instruments code:
1 = fixed
2 = binary
3 = floating point
4 ---- N = optional use 118 

igc:       instrument gain constant 120 

igi:       instrument early or initial gain 122 

corr:      correlated:
1 = no
2 = yes 124

sfs:       sweep frequency at start 126 

sfe:       sweep frequency at end 128 

slen:      sweep length in ms 130 

styp:      sweep type code:
1 = linear
2 = cos-squared
3 = other 132

stas:      sweep trace length at start in ms 134 

stae:      sweep trace length at end in ms 136 

tatyp:     taper type: 1=linear, 2=cos^2, 3=other 138 

afilf:     alias filter frequency if used 140 

afils:     alias filter slope 142 

nofilf:    notch filter frequency if used 144 

nofils:    notch filter slope 146 

lcf:       low cut frequency if used 148 

hcf:       high cut frequncy if used 150 

lcs:       low cut slope 152 

hcs:       high cut slope 154 

year:      year data recorded 156 

day:       day of year 158 

hour:      hour of day (24 hour clock) 160 

minute:    minute of hour 162 

sec:       second of minute 164 

timbas:    time basis code:
1 = local
2 = GMT
3 = other 166

trwf:      trace weighting factor, defined as 1/2^N
volts for the least sigificant bit 168 

grnors:    geophone group number of roll switch
position one 170

grnofr:    geophone group number of trace one within
original field record 172

grnlof:    geophone group number of last trace within
original field record 174

gaps:      gap size (total number of groups dropped) 176 

otrav:     overtravel taper code: 
1 = down (or behind)
2 = up (or ahead) 178

cdpx:   X coordinate of CDP 180

cdpy:   Y coordinate of CDP 184

iline:  in-line number 188 

xline:  cross-line number 192

shnum:  shotpoint number 196

shsca:  shotpoint scalar 200

tval:   trace value meas. 202

tconst4: transduction const 204

tconst2: transduction const 208

tunits:  transduction units 210

device:  device identifier 212

tscalar: time scalar 214

stype:   source type 216

sendir:  source energy dir. 218
 
unknown: unknown 222

smeas4:  source measurement 224

smeas2:  source measurement 228

smeasu:  source measurement unit 230 

unass1:  unassigned 232

unass2:  unassigned 236

additional keys can be created in the output trace headers
    The parameters key1, key2, ... key# are used to additional keys.   
    The keys must be unique and different from the SEGY key names above.
    The input header byte location and lengths must also be defined.  
    This capability is described in an example that defines a new keys 
    iline1 from byte 220 and xline1 from byte 224:
	  key1=iline1 iline1=220 key1_len=4 \\ 
          key2=xline1 xline1=224 key2_len=4 \\

    key#_len defaults to 4
*/
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

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>

#include <stdio.h>

#include <rsf.h>

#include "segy.h"

int main(int argc, char *argv[])
{
    bool verbose, su, xdr, suxdr;
    const char *read, *headname;
    char ahead[SF_EBCBYTES], bhead[SF_BNYBYTES];
    char *filename, *trace, *prog, key[7], *name;
    sf_file out, hdr, msk=NULL;
    int format, ns, itr, ntr, n2, itrace[SF_MAXKEYS], *mask, nkeys=SF_NKEYS, ik;
    off_t pos, start, nsegy=0;
    FILE *head, *file;
    float *ftrace, dt=0.0, t0;
    extern int fseeko(FILE *stream, off_t offset, int whence);
    extern off_t ftello (FILE *stream);

    sf_init(argc, argv);

    if (!sf_getbool("verb",&verbose)) verbose=false;
    /* Verbosity flag */

    if (!sf_getbool("su",&su)) {
	/* y if input is SU, n if input is SEGY */
	prog = sf_getprog();
	if (NULL != strstr (prog, "suread")) {
	    su = true;
	} else if (NULL != strstr (prog, "segyread")) {
	    su = false;
	} else {
	    sf_warning("%s is neither suread nor segyread, assume segyread",
		       prog);
	    su = false;
	}
    }

    if (su) {
	if (!sf_getbool("suxdr",&suxdr)) suxdr=false;
	/* y, SU has XDR support */
    } else {
	suxdr = true;
    }

    if (!sf_getbool("endian",&xdr)) xdr=true;
    /* Whether to automatically estimate endianness or not */
    if (xdr) endian();

    if (NULL == (filename = sf_getstring("tape"))) {
	/* input data */ 
	file = stdin;
    } else if (NULL == (file = fopen(filename,"rb"))) {
	sf_error("Cannot open \"%s\" for reading:",filename);
    }

    if (-1 == ftello(file)) sf_error("Cannot read from a pipe");

    if (!sf_getint("n2",&ntr)) ntr=0;
    /* number of traces to read (if 0, read all traces) */

    if (0==ntr) {
	fseeko(file,0,SEEK_END);
	pos = ftello(file); /* pos is the filesize in bytes */
	fseeko(file,0,SEEK_SET);
    } else {
	pos = 0;
    }

    /* figure out the number of trace keys */
    for (ik=0; nkeys < SF_MAXKEYS; ik++, nkeys++) {
	snprintf(key,7,"key%d",ik+1);
	if (NULL == (name = sf_getstring(key)) || 
	    NULL == sf_getstring(name)) break;
	/* \n
	  ( key# extra key for trace headers )
	  For example to define a new keys iline1 from byte 220 and xline1
	  from byte 224:
	  key1=iline1 iline1=220 key1_len=4 
          key2=xline1 xline1=224 key2_len=4

          key#_len defaults to 4
	*/
    }
    segy_init(nkeys,NULL);

    if (!su) { /* read ascii and binary headers */
	if (SF_EBCBYTES != fread(ahead, 1, SF_EBCBYTES, file)) 
	    sf_error("Error reading ebcdic header");
	
	ebc2asc (SF_EBCBYTES, ahead);
	
	if (NULL == (headname = sf_getstring("hfile"))) {
	    /* output text data header file */
	    head = NULL;
	} else {
	    if (NULL == (head = fopen(headname,"w")))
	    sf_error("Cannot open file \"%s\" for writing ascii header:",
		     headname);

	    if (SF_EBCBYTES != fwrite(ahead, 1, SF_EBCBYTES, head)) 
	    sf_error("Error writing ascii header");

	    fclose (head);
	    
	    if (verbose) 
	    sf_warning("ASCII header written to \"%s\"",headname);
	}

	if (SF_BNYBYTES != fread(bhead, 1, SF_BNYBYTES, file))
	    sf_error("Error reading binary header");

	if (NULL == (headname = sf_getstring("bfile"))) {
	    /* output binary data header file */
	    head = NULL;
	} else {
	    if (NULL == (head = fopen(headname,"wb"))) 
	    sf_error("Cannot open file \"%s\" for writing binary header:",
		     headname);

	    if (SF_BNYBYTES != fwrite(bhead, 1, SF_BNYBYTES, head)) 
	    sf_error("Error writing binary header");
	 
	    fclose (head);

	    if (verbose) 
	    sf_warning("Binary header written to \"%s\"",headname);
	}

	if (!sf_getint("format",&format)) format = segyformat (bhead);
	/* [1,2,3,5] Data format. The default is taken from binary header.
	   1 is IBM floating point
	   2 is 4-byte integer
	   3 is 2-byte integer
	   5 is IEEE floating point
       6 is native_float (same as RSF binary default)
	*/

	switch (format) {
	    case 1:
		if (verbose) sf_warning("Assuming IBM floating point format");
		break;
	    case 2:
		if (verbose) sf_warning("Assuming 4 byte integer format");
		break;
	    case 3:
		if (verbose) sf_warning("Assuming 2 byte integer format");
		break;
	    case 5:
		if (verbose) sf_warning("Assuming IEEE floating point format");
		break;
        case 6:
        if (verbose) sf_warning("Assuming native_float format");
        suxdr=false;
        break;
	    default:
		sf_error("Nonstandard format: %d",format);
		break;
	}

	if (!sf_getint("ns",&ns)) ns = segyns (bhead);
	/* Number of samples. The default is taken from binary header */
	if (0>=ns) sf_error("Number of samples is not set in binary header");

	if (verbose) sf_warning("Detected trace length of %d",ns);

	dt = segydt (bhead);
	nsegy = SF_HDRBYTES + ((3 == format)? ns*2: ns*4);    
	if (0==ntr) ntr = (pos - SF_EBCBYTES - SF_BNYBYTES)/nsegy;

        start = ftello(file);
    } else {
        start = 0;
    }

    /* read first trace header */
    trace = sf_charalloc (SF_HDRBYTES);
    if (SF_HDRBYTES != fread(trace, 1, SF_HDRBYTES, file))
	sf_error ("Error reading first trace header");
    fseeko(file,start,SEEK_SET); /* rewind */

    segy2head(trace, itrace, SF_NKEYS);
    t0 = itrace[ segykey("delrt")]/1000.;

    if (su) { /* figure out ns, dt, and ntr */
	if (!sf_getint("ns",&ns)) ns = itrace[ segykey("ns")];
	if (verbose) sf_warning("Detected trace length of %d",ns);
	dt = itrace[ segykey("dt")]/1000000.;
	free (trace);

	nsegy = SF_HDRBYTES + ns*4;
	if (0==ntr) ntr = pos/nsegy;	

	if (suxdr) format=5;
    } 

    if (verbose) sf_warning("Expect %d traces",ntr);

    if (NULL != sf_getstring("mask")) {
	/* optional header mask for reading only selected traces */
	msk = sf_input("mask");
	if (SF_INT != sf_gettype(msk)) sf_error("Need integer mask");

	mask = sf_intalloc(ntr);
	sf_intread(mask,ntr,msk);
	sf_fileclose(msk);

	for (n2=itr=0; itr < ntr; itr++) {
	    if (mask[itr]) n2++;
	}
    } else {
	mask = NULL;
	n2 = ntr;
    }
 
    if (NULL == (read = sf_getstring("read"))) read = "b";
    /* what to read: h - header, d - data, b - both (default) */

    if (read[0] != 'h') { /* not only header */
	out = sf_output("out");
	sf_putint(out,"n1",ns);
	sf_putint(out,"n2",n2);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",t0);
	sf_putstring(out,"label1","Time");
	sf_putstring(out,"unit1","s");
	sf_putfloat(out,"d2",1.);
	sf_putfloat(out,"o2",0.);
	sf_putstring(out,"label2","Trace");
	sf_setformat(out, "native_float");    
	ftrace = suxdr? sf_floatalloc (ns): NULL;
    } else {
	out = NULL;
	ftrace = NULL;
    }
    
    if (read[0] != 'd') { /* not only data */
	hdr = sf_output("tfile");

	sf_putint(hdr,"n1",nkeys);
	sf_putint(hdr,"n2",n2);
	sf_setformat(hdr,"native_int");

	segy2hist(hdr,nkeys);

	if (NULL == (headname = sf_getstring("tfile"))) headname = "tfile";
	/* output trace header file */
	if (NULL != out) sf_putstring(out,"head",headname);
    } else {
	hdr = NULL;
    }

    if (NULL != out) sf_fileflush(out,NULL);

    switch (read[0]) {
	case 'h': /* header only */
	    trace = sf_charalloc (SF_HDRBYTES);
	    nsegy -= SF_HDRBYTES;

	    for (itr=0; itr < ntr; itr++) {
		if (NULL != mask && !mask[itr]) {
		    fseeko(file,SF_HDRBYTES,SEEK_CUR);
		    continue;
		} else if (SF_HDRBYTES != fread(trace, 1, SF_HDRBYTES, file)) {
			sf_error ("Error reading trace header %d",itr+1);
		}
		fseeko(file,nsegy,SEEK_CUR);

		segy2head(trace, itrace, nkeys);
		sf_intwrite(itrace, nkeys, hdr);
	    }

	    break;
	case 'd': /* data only */
	    nsegy -= SF_HDRBYTES;		    
	    trace = sf_charalloc (nsegy);

	    for (itr=0; itr < ntr; itr++) {
		fseeko(file,SF_HDRBYTES,SEEK_CUR);
		if (NULL != mask && !mask[itr]) {
		    fseeko(file,nsegy,SEEK_CUR);
		    continue;
		} else if (nsegy != fread(trace, 1, nsegy, file)) {
			sf_error ("Error reading trace data %d",itr+1);
		}

		if (suxdr) {
		    segy2trace(trace, ftrace, ns,format);
		    sf_floatwrite (ftrace,ns,out);
		} else {
		    sf_charwrite (trace,ns*sizeof(float),out);
		} 
	    }

	    break;
	default: /* both header and data */
	    trace = sf_charalloc (nsegy);

	    for (itr=0; itr < ntr; itr++) {
		if (NULL != mask && !mask[itr]) {
		    fseeko(file,nsegy,SEEK_CUR);
		    continue;
		} else if (nsegy != fread(trace, 1, nsegy, file)) {
		    sf_error ("Error reading trace header %d",itr+1);
		}

		segy2head(trace, itrace, nkeys);
		sf_intwrite(itrace, nkeys, hdr);

		if (suxdr) {
		    segy2trace(trace + SF_HDRBYTES, ftrace, ns,format);
		    sf_floatwrite (ftrace,ns,out);
		} else {
		    sf_charwrite (trace + SF_HDRBYTES,ns*sizeof(float),out);
		} 
	    }
	    break;
    }


    exit (0);
}

/* 	$Id$	 */
