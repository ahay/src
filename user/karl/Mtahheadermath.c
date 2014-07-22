/* Mathematical operations with headers.  Output whole header.

Known functions for float data: 
cos,  sin,  tan,  acos,  asin,  atan, 
cosh, sinh, tanh, acosh, asinh, atanh,
exp,  log,  sqrt, abs, erf, erfc, sign

Known functions for int data: sign, abs

modeled on sfheademath.  This program outputs the whole header.  It was
a building block for sftahheadermath.

tah is the abbreviation of Trace And Header.  Madagascar programs 
that begin with sftah are a designed to:
1- read trace and headers from separate rsf files and write them to 
   standard output (ie sftahread)
2- filter programs that read and write standard input/output and 
   process the tah data (eg sftahnmo, sftahstack)
3- read tah data from standard input and write separate rsf files for 
   the trace and headers data (ie sftahwrite)

These programs allow Seismic Unix (su) like processing in Madagascar.  
Some programs have su like names.

Some programs in this suite are sftahread, sftahgethw, ftahhdrmath, 
and sftahwrite.

The sftahheadermath updates a trace header with a new value computed from
input trace headers.  

See also sftahmakeskey.

EXAMPLE:

sftahread \\
   verbose=1 \\
   input=npr3_gathers.rsf \\
| sftahgethw \\
   verbose=0  \\
   key=sx,sy,gx,gy,offset  \\
>/dev/null

The headers are in the file npr3_gathers_hdr.rsf, 
the headers parameter default.  The headers are merged with the trace 
amplitudes and the tah data sent down the pipe for sftahgethw.  The 
source and group coordinates and offset (sx,sy,gx,gy,offset) are 
printed to STDERR.  Traces are sent to STDOUT, which is directed to
/dev/null (the bit bucket).

PARAMETERS
   string output= no default

        An equation to compute using the header keys.  Equations should
	problable be enclosed in quotes, ie ", to the equation can include
	multiplication, *, or blanks.  
	For example, to compute the midpoint x input:
	output="(sx+gx)/2.0)"

   string outoutkey= no default
        the name of the output trace header key to put the evaluation of
	output.  For example to put the average of sx and gx into cdpx input:
	outputkey=cdpx

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

#include <string.h>
#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"

int main(int argc, char* argv[])
{
  int verbose;
  int i, i1, n1_headers, n1_traces, len, tempint, outkeyindx;
    sf_file in, out;
    char *output=NULL, *outputkey=NULL;
    float *ftra=NULL, **fbuf=NULL, **fst=NULL;
    int *itra=NULL, **ibuf=NULL, **ist=NULL;
    char* header_format;  
    sf_datatype typehead;

    float* fheader=NULL;
    int* iheader=NULL;
    float* intrace=NULL;
    

    sf_init (argc,argv);
    /*****************************/
    /* initialize verbose switch */
    /*****************************/
    if(!sf_getint("verbose",&verbose))verbose=1;
    /* \n
       flag to control amount of print
       0 terse, 1 informative, 2 chatty, 3 debug
    */
    sf_warning("verbose=%d",verbose);
 
    /******************************************/
    /* input and output data are stdin/stdout */
    /******************************************/


    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint(in,"n1_traces",&n1_traces))
      sf_error("input data not define n1_traces");
    if (!sf_histint(in,"n1_headers",&n1_headers)) 
      sf_error("input data does not define n1_headers");

    /* kls change type to header_format and read from history */
    header_format=sf_histstring(in,"header_format");
    if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
    else                                       typehead=SF_FLOAT;

    if(verbose>0)fprintf(stderr,"allocate headers. n1_headers=%d\n",n1_headers);
    fheader = sf_floatalloc(n1_headers);
    iheader = (int*)fheader;

    if(verbose>0)fprintf(stderr,"allocate intrace. n1_traces=%d\n",n1_traces);
    intrace= sf_floatalloc(n1_traces);

    segy_init(n1_headers,in);
 
    for (i=0; i < n1_headers; i++) {
      /* see if the segy keywords are in the input history file.  If they
	 are missing or different than I think they should be add them to
	 the output file history */
      if(!sf_histint(in,segykeyword(i),&tempint) || tempint!=i){
	sf_putint(out,segykeyword(i),i);
      }
    }


    if (NULL == (output = sf_getstring("output"))) sf_error("Need output=");
    /* Describes the output in a mathematical notation. */

    len = sf_math_parse (output,out,typehead);

    if (NULL==(outputkey=sf_getstring("outputkey")))sf_error("Need outputkey=");
    /* name of the header key to put the results of the output equation */
    if(!sf_histint(out,outputkey,&outkeyindx)){
      sf_error("need outkey is not a header key in the input data=");
    }
    if(verbose>0)fprintf(stderr,"outkeyindx=%d\n",outkeyindx);

    /* I do not like these 2d arrays with one of the lengths is 1.
       I have done this so I can continue to use sf_math_parse and 
       sf_math_evaluate.  Perhaps someday these can be refactored and the 
       alloc2 below can become sf_floatalloc (without the 2). Karl S */ 
    if (SF_FLOAT == typehead) { /* float typehead */
	ftra = sf_floatalloc(n1_headers);
	fbuf = sf_floatalloc2(1,n1_headers);
	fst  = sf_floatalloc2(1,len+3);
    } else {               /* int typehead */
	itra = sf_intalloc(n1_headers);
	ibuf = sf_intalloc2(1,n1_headers);
	ist  = sf_intalloc2(1,len+3);
    }

    /* put the history from the input file to the output */
    sf_fileflush(out,in);

    /***************************/
    /* start trace loop        */
    /***************************/
    if(verbose>0)fprintf(stderr,"start trace loop\n");
    while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
      if(verbose>1)fprintf(stderr,"process the tah in sftahgethw\n");

      /********************/
      /* process the tah. */
      /********************/
 
      if (SF_FLOAT == typehead) { 
	for (i1=0; i1 < n1_headers; i1++) {
	  fbuf[i1][0]=fheader[i1];
	}
	sf_math_evaluate (len, 1, fbuf, fst);
        if(verbose>2){
	  fprintf(stderr,"after math_evaluate fst[1][0]=%f\n",fst[1][0]);
	}
	fheader[outkeyindx]=fst[1][0];	
      } else {
	for (i1=0; i1 < n1_headers; i1++) {
	  /* iheader point to same place as fheader */
	  ibuf[i1][0]=iheader[i1];
	  if(verbose>2)fprintf(stderr,"iheader[i1]=%d\n",iheader[i1]);
	}
	sf_int_math_evaluate (len, 1, ibuf, ist);
        if(verbose>2){
	  fprintf(stderr,"after int_math_evaluate ist[1][0]=%d\n",ist[1][0]);
	}
	iheader[outkeyindx]=ist[1][0];
      }
      
      /***************************/
      /* write trace and headers */
      /***************************/
      put_tah(intrace, fheader, n1_traces, n1_headers, out);
      if(verbose>1)fprintf(stderr,"returned from writing the tah\n");
	
    }

    exit(0);
}
