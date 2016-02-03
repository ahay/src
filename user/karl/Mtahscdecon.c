/* Trace And Header Surface Consistant Decon.

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

The sftahsfdecon program applies surface consistant decon.

EXAMPLE:

sftahsort \\
   input=npr3_field.rsf \\
   sort="sx,sy" \\
| sftahgain  \\
   tpow=2 \\
| sftahmute \\
   tmute=0.0,3.0 \\
   xmute=0,18000  \\
| sftahscdecon \\
   key="scaon sx,sy" \\
   length=140 \\
   pmoise=.001 \\
   verbose=0  \\
|  sftahwrite output=shotdecon.rsf \\
   mode=seq \\
>/dev/null

sftahsort \\
   input=shotdecon.rsf \\
   sort="gx,gy"
| sftahscdecon \\
   length=140 \\
   pmoise=.001 \\
   verbose=0  \\
|  sftahwrite output=s-g-decon.rsf \\
   mode=seq \\
>/dev/null


Traces are in the file npr_field.rsf anf headers in the npr_field_hdr.rsf 
file.  The headers are merged with the trace amplitudes and the tah data 
sent down the pipe for spreading loss correction (sftahgain tpow=2),
a pre decon mute (sftahmute) and shot consistant decon.  Data is written
to the output file shotdecon.rsf and headers to shotdecon_hdr.rsf.  In
the next sequence the data is sorted in to common receiver domain so 
receiver consistent decon can be applied.

PARAMETERS
   length=140 
   pmoise=.001 
   verbose=0  

*/

/*
  Copyright (C) 2013 University of Texas at Austin
  
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
/*
   Program change history:
   date       Who             What
   09/22/2014 Karl Schleicher Original program based on sftah5dinterp
*/
#include <string.h>
#include <rsf.h>
#include <rsf_su.h>
#include <rsfsegy.h>

#include "tahsub.h"
#include "extraalloc.h"

int main(int argc, char* argv[])
{
  int verbose;
  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;

  char* header_format=NULL;
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  float* intrace=NULL;
  int numkeys;
  int ikey;
  char** list_of_keys;
  int *indx_of_keys;
  int indx_sx, indx_sy, indx_gx, indx_gy;
  int indx_iline, indx_xline, indx_offset;
  int indx_cdpx, indx_cdpy;
  int gather_maxfold;
  float** gather_traces;
  float** gather_headers;
  int itrace;
  int iitrace;
  bool eof_get_tah;
  bool end_of_gather;
  float sx, sy, gx, gy, offx, offy;
  int fold;
  int maxfold=0;
  float* outtrace=NULL;
  float* outheader=NULL;  
  int* ioutheader=NULL; 
  float* first_gather_header=NULL;

  /* variables from sfpef (or sftahpef or supef) */
  int ilag;		/* counters				*/
  float dt;		/* time sample interval (sec)		*/
  float *wiener=NULL;	/* Wiener error filter coefficients	*/
  float *spiker=NULL;	/* spiking decon filter			*/

  float pnoise;		/* pef additive noise level		*/

  float minlag;	/* start of error filter (sec)		*/
  int iminlag;	/* ... in samples			*/
  float maxlag;	/* end of error filter (sec)		*/
  int imaxlag;	/* ... in samples			*/
  int nlag;		/* length of error filter in samples	*/
  int ncorr;		/* length of corr window in samples	*/
  int lcorr;		/* length of autocorr in samples	*/

  float *crosscorr=NULL;	/* right hand side of Wiener eqs	*/
  float *autocorr=NULL;	/* vector of autocorrelations		*/

  float mincorr;		/* start time of correlation window	*/
  int imincorr;		/* .. in samples			*/
  float maxcorr;		/* end time of correlation window	*/
  int imaxcorr;		/* .. in samples			*/
	
  size_t lagbytes;	/* bytes in wiener and spiker filters	*/
  size_t maxlagbytes;	/* bytes in autocorrelation		*/
  float *temp=NULL;	/* temporary array			*/
  sf_file wien;
 
  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  sf_init (argc,argv);

  if(!sf_getint("verbose",&verbose))verbose=1;
      /* \n
     flag to control amount of print
     0 terse, 1 informative, 2 chatty, 3 debug
  */
  sf_warning("verbose=%d",verbose);
 
  /******************************************/
  /* input and output data are stdin/stdout */
  /******************************************/

  if(verbose>0)fprintf(stderr,"read in file name\n");  
  in = sf_input ("in");

  if(verbose>0)fprintf(stderr,"read out file name\n");
  out = sf_output ("out");

  if (!sf_histint(in,"n1_traces",&n1_traces))
    sf_error("input data not define n1_traces");
  if (!sf_histint(in,"n1_headers",&n1_headers)) 
    sf_error("input data does not define n1_headers");

  header_format=sf_histstring(in,"header_format");
  if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
  else                                       typehead=SF_FLOAT;

  if(verbose>0)fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  fheader = sf_floatalloc(n1_headers);
  first_gather_header = sf_floatalloc(n1_headers);
 
  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);

  /* maybe I should add some validation that n1== n1_traces+n1_headers+2
     and the record length read in the second word is consistent with 
     n1.  */

  /**********************************************************/
  /* end code block for standard tah Trace And Header setup */
  /* continue with any sf_puthist this tah program calls to */
  /* add to the history file                                */
  /**********************************************************/

  if(verbose>0)fprintf(stderr,"call list of keys\n");
 
  list_of_keys=sf_getnstring("key",&numkeys);
  if(list_of_keys==NULL)
    sf_error("The required parameter \"key\" was not found.");
  /* I wanted to use sf_getstrings, but it seems to want a colon seperated
     list of keys (eg key=offset:ep:fldr:cdp) and I wanted a comma seperated
     list of keys (eg key=offset:ep:fldr:cdp).
  numkeys=sf_getnumpars("key");
  if(numkeys==0)
    sf_error("The required parameter \"key\" was not found.");
  fprintf(stderr,"alloc list_of_keys numkeys=%d\n",numkeys);
  list_of_keys=(char**)sf_alloc(numkeys,sizeof(char*)); 
  sf_getstrings("key",list_of_keys,numkeys);
  */
  /* print the list of keys */
  if(verbose>1){
    fprintf(stderr,"numkeys=%d\n",numkeys);
    for(ikey=0; ikey<numkeys; ikey++){
      fprintf(stderr,"list_of_keys[%d]=%s\n",ikey,list_of_keys[ikey]);
    }
  }

  if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
  
  if (!sf_getfloat("minlag",&minlag)) {
    /* first lag of prediction filter (sec) */
    minlag=0.0f;
    iminlag = 1;
  } else {
    iminlag = SF_NINT(minlag/dt);
    if(iminlag<1)iminlag=1;
  }
  if (!sf_getfloat("maxlag",&maxlag)) {
    /* last lag of prediction filter (sec) */
    maxlag=0.0f;
    imaxlag = SF_NINT(0.05 * n1_traces); /* kls use n1_traces instead of nt */
  } else {
    imaxlag = SF_NINT(maxlag/dt);
  }
  
  /* Get parameters  */
  if (!sf_getfloat("pnoise",  &pnoise)) pnoise = 0.001;
  /* relative additive noise level */
  
  /* .. mincorr and maxcorr */
  if (sf_getfloat("mincorr", &mincorr)) imincorr = SF_NINT(mincorr/dt);
  /*(mincorr start of autocorrelation window in sec )*/
  else				  imincorr = 0;
  if (imincorr < 0) sf_error("mincorr=%g too small", mincorr);
  
  if (sf_getfloat("maxcorr", &maxcorr)) imaxcorr = SF_NINT(maxcorr/dt);
  /*(maxcorr end of autocorrelation window in sec )*/
  else				  imaxcorr = n1_traces;
  if (imaxcorr > n1_traces) sf_error("maxcorr=%g too large", maxcorr);
  
  if (imincorr >= imaxcorr)
    sf_error("mincorr=%g, maxcorr=%g", mincorr, maxcorr);

  /* compute filter sizes and correlation number */
  nlag  = imaxlag - iminlag + 1;
  ncorr = imaxcorr - imincorr + 1;
  lcorr = imaxlag + 1;
  
  /* Compute byte sizes in wiener/spiker and autocorr */
  lagbytes = sizeof(float)*nlag;
  maxlagbytes = sizeof(float)*lcorr;
  
  /* Allocate memory */
  wiener	 = sf_floatalloc(nlag);
  spiker	 = sf_floatalloc(nlag);
  autocorr = sf_floatalloc(lcorr);
  temp = sf_floatalloc(lcorr);
  /* Set pointer to "cross" correlation */
  crosscorr = autocorr + iminlag;

  if (NULL != sf_getstring("wiener")) {
    /* file to output Wiener filter */
    wien = sf_output("wiener");
    sf_putint(wien,"n1",nlag);
  } else {
    wien = NULL;
  }

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************/
  /* continue initialization specific to this tah program */
  /********************************************************/

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);
  indx_of_keys=sf_intalloc(numkeys);
  for (ikey=0; ikey<numkeys; ikey++){
    /* kls need to check each of these key names are in the segy header and
       make error message for invalid keys.  Of does segykey do this? NO, just
       segmentation fault. */
    indx_of_keys[ikey]=segykey(list_of_keys[ikey]);
  }

  /* I used Mtahstack.c and Mtahgethw to for inspiration for data flow */

  /* more initalization before starting the trace read loop */
  indx_sx=segykey("sx");
  indx_sy=segykey("sy");
  indx_gx=segykey("gx");
  indx_gy=segykey("gy");
  indx_iline=segykey("iline");
  indx_xline=segykey("xline");
  indx_offset=segykey("offset");
  indx_cdpx=segykey("cdpx");
  indx_cdpy=segykey("cdpx");

  /* allocate processing arrays */
  gather_maxfold=1000;
  gather_traces =sf_floatalloc2(n1_traces , gather_maxfold);
  gather_headers=sf_floatalloc2(n1_headers, gather_maxfold);
  
  /* allocate output trace arrays */
  outtrace  = sf_floatalloc(n1_traces);
  outheader = sf_floatalloc(n1_headers);
  ioutheader=(int*)outheader;

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");

  if(verbose>0)fprintf(stderr,"start the the trace read loop\n");
  itrace=0;
  eof_get_tah=get_tah(intrace, fheader, n1_traces, n1_headers, in);
  fold=0;
  while (!eof_get_tah){
    if(verbose>1 || (verbose==1 && itrace<5)){
      fprintf(stderr,"process the tah in sftahscdecon\n");
      for(ikey=0; ikey<numkeys; ikey++){
	fprintf(stderr," %s=",list_of_keys[ikey]);
	if(typehead == SF_INT){
	  /* just cast the header to int so the print works */
	  fprintf(stderr,"%d",((int*)fheader)[indx_of_keys[ikey]]);
	} else {
	  fprintf(stderr,"%f",fheader[indx_of_keys[ikey]]);
	}
      } /* end of the for(ikey..) loop */
      fprintf(stderr,"\n");
    }
    /********************/
    /* process the tah. */
    /********************/
    /* loop reading traces in gather */
    
    if(fold==0){
      if(verbose>1)fprintf(stderr,"start a new supergather\n");
      /* apply any initialization */
      memcpy(first_gather_header,fheader,n1_headers*sizeof(int));
    }
    if(fold>=gather_maxfold){
      /* increase buffer sizes */
      if(verbose>1){
	fprintf(stderr,"increase buffer size maxfold=%d\n",gather_maxfold);
	for(ikey=0; ikey<2;ikey++){ 
	  if(typehead == SF_INT){
	    /* just cast the header to int so the print works */
	    fprintf(stderr,"%d ",((int*)fheader)[indx_of_keys[ikey]]);
	  } else {
	    fprintf(stderr,"%f ",fheader[indx_of_keys[ikey]]);
	  }
	}
	fprintf(stderr,"\n");
      }
      gather_maxfold*=1.2;
      gather_traces =sf_floatrealloc2(gather_traces, 
					n1_traces , gather_maxfold);
      gather_headers=sf_floatrealloc2(gather_headers, 
					n1_headers, gather_maxfold);
    }
    memcpy(gather_headers[fold],fheader,n1_headers*sizeof(int));
    memcpy(gather_traces [fold],intrace,n1_traces*sizeof(int));
    fold++;

    /* get the next trace so you can test the current and next headers
       to determine if this is the end of a gather */
    eof_get_tah=get_tah(intrace, fheader, n1_traces, n1_headers, in);

    /* this is the end of gather if eof was encounterred or if at least one 
       of the header keys in indx_of_keys changed. */
    end_of_gather=eof_get_tah;

    /* if there is no next trace, do not try to test the header keys.  
       You already know it is the end of the gather. */
    if(!eof_get_tah){ 
      if(fold>1){
	for(ikey=0; ikey<numkeys; ikey++){
	    if(typehead == SF_INT){
	      if(((int*)fheader  )[indx_of_keys[ikey]]!=
		 ((int*)first_gather_header)[indx_of_keys[ikey]]){
		end_of_gather=true;
		break;
	      }
	    } else {
	      /* kls need to repeat the updates from previous code segment */
	      if(fheader            [indx_of_keys[ikey]]!=
		 first_gather_header[indx_of_keys[ikey]]){
		end_of_gather=true;
		break;
	      }
	    } /* end if(typehead == SF_INT)...else */
	  } /* end for(ikey=0; ikey<0*numkeys; ikey++) */
      }
    } /* end if(!eof_get_tah) the code segment to test hdrs for key change */


    /* if this is the last trace in the gather, process the gather and
       set fold=0.  Fold=0 will cause gather initialization on the next
       iteration through the loop */
    if(end_of_gather){
      /***********************************/
      /* process the supergather         */
      /***********************************/
      /* OK.  You have a gather in memory.  Do your 5D interpolation. */
      if(verbose>1)fprintf(stderr,"end_of_gather.  process it.\n");

      if(fold>maxfold)maxfold=fold;
      if(verbose>1)fprintf(stderr,"fold=%d maxfold=%d\n",fold,maxfold);

      /* Do surface consistant decon. */

      /*****************************************************************/
      /* Design the filter for this gather.  First compute the average */
      /* autocorrelation then compute the filter using weiner-levinson */
      /* recusion (the stoepf function). Apply filter and write trace  */
      /* in a seperate loop                                            */
      /*****************************************************************/
 
      /* zero out filter vectors */
      if(verbose>1)fprintf(stderr,"zero the filter vectors\n");
      memset((void *) wiener, 0, lagbytes);
      memset((void *) spiker, 0, lagbytes);
      memset((void *) autocorr, 0, maxlagbytes);

      /* loop to compute average autocorrelation */

      if(verbose>1)fprintf(stderr,"loop to comput average autocorrelation\n");
      for (iitrace=0; iitrace<fold; iitrace++){
	memset((void *) temp, 0, maxlagbytes);
	xcor(ncorr, 0, gather_traces[iitrace]+imincorr,
	     ncorr, 0, gather_traces[iitrace]+imincorr,
	     lcorr, 0, temp);
	/* the trace autocor is in temp.  Sum temp into the gather autocorr */
	for (ilag = 0; ilag < lcorr; ++ilag) {
	  autocorr[ilag] += temp[ilag]; 
	}

      }
      /* drop gather if autocorr[0] vanishes */
      if(autocorr[0] == 0.0) {
	fprintf(stderr,"*************************************************\n");
	fprintf(stderr,"*************************************************\n");
	fprintf(stderr,"zero autocorrelation.  gather dropped from output\n");
	fprintf(stderr,"*************************************************\n");
	fprintf(stderr,"*************************************************\n");
      } else {
	/* Whiten */
	autocorr[0] *= 1.0 + pnoise;
	
	/* Get inverse filter by Wiener-Levinson */
	if(verbose>1)fprintf(stderr,"call stoepf\n");
	stoepf(nlag, autocorr, crosscorr, wiener, spiker);
	
	/* Show pefilter on request */
	if (NULL != wien) sf_floatwrite(wiener,nlag,wien);
	
	
	/* now you need loop to apply filter, and write trace and header 
	   (tah) to the output pipe with put_tah */
	/*************************************************/
	/* apply filter and write tah (trace and header) */
	/*************************************************/
	
	for (iitrace=0; iitrace<fold; iitrace++){
	  int i;  /* index for the output trace sample */
	  /* Convolve pefilter with trace - don't do zero multiplies */
	  if(verbose>2)fprintf(stderr,"apply filter iitrace=%d\n",iitrace);
	  for (i = 0; i < n1_traces; ++i) {
	    register int j;
	    register int n = SF_MIN(i, imaxlag); 
	    register float sum = gather_traces [iitrace][i];
	    
	    for (j = iminlag; j <= n; ++j)
	      sum -= wiener[j-iminlag] * gather_traces [iitrace][i-j];
	    
	    outtrace[i] = sum;
	  }
	  if(verbose>1)fprintf(stderr,"put_tah\n");
	  put_tah(outtrace, gather_headers[iitrace], 
		  n1_traces, n1_headers, out);
	  if(verbose>1){
	    fprintf(stderr,"end apply/write loop iitrace=%d\n",iitrace);
	  }
	}
      }
      fold=0;
    }
    itrace++;
  }

  exit(0);
}

  
