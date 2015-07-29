/* Convert an RSF dataset to SEGY or SU.

Merges trace headers with data.

"suwrite" is equivalent to "segywrite su=y"

If bfile= and/or hfile= are not provided, they will be created automatically
using information from the trace headers.

The file for tfile= can be generated with sfsegyheader.
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

#include <stdio.h>
#include <rsf.h>
#include "segy.h"

int main(int argc, char *argv[])
{
    bool verbose, su, xdr, suxdr;
    char ahead[SF_EBCBYTES], bhead[SF_BNYBYTES];
    char *headname=NULL, *filename=NULL, *trace=NULL, count[4], *prog=NULL;
    const char *myheader[] = {"      This dataset was created",
			      "     with the Madagascar package",
			      "     http://www.ahay.org/"};
    sf_file in=NULL, hdr=NULL;
    size_t nsegy;
    int format=1, i, ns, nk, itr, ntr, *itrace=NULL;
    FILE *head=NULL, *file=NULL;
    float *ftrace=NULL, dt;

    sf_init(argc, argv);
    in = sf_input("in");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_getbool("verb",&verbose)) verbose=false;
    /* Verbosity flag */

    if (!sf_getbool("endian",&xdr)) xdr=true;
    /* Whether to automatically estimate endianness or not */
    if (xdr) endian();

    if (!sf_getbool("su",&su)) {
	/* y if input is SU, n if input is SEGY */
	prog = sf_getprog();
	if (NULL != strstr (prog, "suwrite")) {
	    su = true;
	} else if (NULL != strstr (prog, "segywrite")) {
	    su = false;
	} else {
	    sf_warning("%s is neither suwrite nor segywrite, assume segywrite",
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

    if (NULL == (filename = sf_getstring("tape"))) {
	/* output data */
	file = stdout;
    } else if (NULL == (file = fopen(filename,"wb"))) {
	sf_error("Cannot open \"%s\" for writing:",filename);
    }

    if (!su) {
	if (NULL != (headname = sf_getstring("hfile")) &&
	    NULL != (head = fopen(headname,"r"))) {
	    /* input text data header file */
	    if (SF_EBCBYTES != fread(ahead, 1, SF_EBCBYTES, head)) 
		sf_error("Error reading ascii header");
	    fclose (head);

	    if (verbose) sf_warning("ASCII header read from \"%s\"",headname);
	} else {
	    for (i=0; i < SF_EBCBYTES/80; i++) {
		snprintf(count,4,"C%-2d",i+1);
		snprintf(ahead+i*80,81,"%s %-76s\n",count,
			 (i < 3)? myheader[i]:"");
	    }
	    if (verbose) sf_warning("ASCII header created on the fly");
	}

	asc2ebc (SF_EBCBYTES, ahead);

	if (SF_EBCBYTES != fwrite(ahead, 1, SF_EBCBYTES, file)) 
	    sf_error("Error writing ebcdic header");

	if (NULL != (headname = sf_getstring("bfile")) &&
	    NULL != (head = fopen(headname,"rb"))) {
	    /* input binary data header file */
	    if (SF_BNYBYTES != fread(bhead, 1, SF_BNYBYTES, head)) 
		sf_error("Error reading binary header");
	    fclose (head);

	    if (verbose) sf_warning("Binary header read from \"%s\"",headname);
	} else {
	    memset(bhead,0,SF_BNYBYTES);
	    set_segyformat(bhead,1);

	    if (verbose) sf_warning("Binary header created on the fly");
	}

	if (sf_histint(in,"n1",&ns)) set_segyns(bhead,ns);
	if (sf_histfloat(in,"d1",&dt)) set_segydt(bhead,dt);
	
	binary_head(bhead);

	if (SF_BNYBYTES != fwrite(bhead, 1, SF_BNYBYTES, file))
	    sf_error("Error writing binary header");

	format = segyformat (bhead);

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
		if (verbose) sf_warning("Assuming IEEE float format");
		break;
	    default:
		sf_error("Nonstandard format: %d",format);
		break;
	}
    }

    if (!sf_histint(in,"n1",&ns)) ns = su? 0: segyns (bhead); 
    if (0 >= ns) sf_error("Failed to determine trace length");

    if (verbose) sf_warning("Detected trace length of %d",ns);

    if (su) {
	if (!sf_getbool("suxdr",&suxdr)) suxdr=false;
	/* y, SU has XDR support */

	if (suxdr) format=5;
    }

    nsegy = SF_HDRBYTES + ((3 == format)? ns*2: ns*4);

    ntr = sf_leftsize(in,1);

    if (verbose) sf_warning("Expect %d traces",ntr);

    hdr = sf_input("tfile");
    if (!sf_histint(hdr,"n1",&nk) || (SF_NKEYS > nk))
	sf_error ("Need at least n1=%d keys in tfile",SF_NKEYS);
    if (nk*ntr != sf_filesize(hdr))
	sf_error ("Wrong number of traces in tfile");
    
    segy_init(nk,hdr);

    trace = sf_charalloc (nsegy);
    ftrace = sf_floatalloc (ns);
    itrace = sf_intalloc (nk); 

    for (itr=0; itr < ntr; itr++) {
	sf_intread(itrace,nk,hdr);	
	head2segy(trace, itrace, SF_NKEYS);
	
	if (suxdr) {
	    sf_floatread (ftrace,ns,in);
	    trace2segy(trace + SF_HDRBYTES, ftrace, ns,format);
	} else {
	    sf_charread(trace + SF_HDRBYTES,ns*sizeof(float),in);
	} 

	if (nsegy != fwrite(trace, 1, nsegy, file))
	    sf_error ("Error writing trace %d",itr+1);
    }

    exit (0);
}
