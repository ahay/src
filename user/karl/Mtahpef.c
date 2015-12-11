/* Wiener predictive error filtering */
/* Copyright (c) Colorado School of Mines, 2011.
   All rights reserved.
  
   Redistribution and use in source and binary forms, with or 
   without modification, are permitted provided that the following 
   conditions are met:
  
   *  Redistributions of source code must retain the above copyright 
   notice, this list of conditions and the following disclaimer.
   *  Redistributions in binary form must reproduce the above 
   copyright notice, this list of conditions and the following 
   disclaimer in the documentation and/or other materials provided 
   with the distribution.
   *  Neither the name of the Colorado School of Mines nor the names of
   its contributors may be used to endorse or promote products 
   derived from this software without specific prior written permission.
  
   Warranty Disclaimer:
   THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
   COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
   IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
   POSSIBILITY OF SUCH DAMAGE.
  
  
   Export Restriction Disclaimer:
   We believe that CWP/SU: Seismic Un*x is a low technology product that does
   not appear on the Department of Commerce CCL list of restricted exports.
   Accordingly, we believe that our product meets the qualifications of
   an ECCN (export control classification number) of EAR99 and we believe
   it fits the qualifications of NRR (no restrictions required), and
   is thus not subject to export restrictions of any variety.
  
   Approved Reference Format:
   In publications, please refer to SU as per the following example:
   Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
   Release No. __: an open source software  package for seismic 
   research and processing, 
   Center for Wave Phenomena, Colorado School of Mines.
  
   Articles about SU in peer-reviewed journals:
   Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
   Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
   Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
   Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
  
   Acknowledgements:
   SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
   School of Mines, partially based on Stanford Exploration Project (SEP) 
   software.
*/
/* SUPEF: $Revision: 1.45 $ ; $Date: 2011/11/16 17:47:47 $		*/

/* Modified for inclusion with Madagascar. */

/* Credits:
 *	CWP: Shuki Ronen, Jack K. Cohen, Ken Larner
 *      CWP: John Stockwell, added mixing feature (April 1998)
 *      CSM: Tanya Slota (September 2005) added cdp feature
 *
 *      Technical Reference:
 *	A. Ziolkowski, "Deconvolution", for value of maxlag default:
 *		page 91: imaxlag < nt/10.  I took nt/20.
 *
 * Notes:
 *	The prediction error filter is 1,0,0...,0,-wiener[0], ...,
 *	so no point in explicitly forming it.
 *
 *	If imaxlag < 2*iminlag - 1, then we don't need to compute the
 *	autocorrelation for lags:
 *		imaxlag-iminlag+1, ..., iminlag-1
 *	It doesn't seem worth the duplicated code to implement this.
 */

#include <rsf.h>
#include <rsf_su.h>
#include <rsfsegy.h>

#include "tahsub.h"

int 
main(int argc, char **argv)
{
    int nt;			/* number of points on trace		*/
    int i,ilag;		/* counters				*/

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
	
    int imix;		/* mixing counter			*/
    int nmix;		/* number of traces to average over	*/
    size_t mixbytes;	/* number of bytes = maxlagbytes*nmix	*/ 
    float *mix=NULL;	/* array of averaging weights		*/
    float **mixacorr=NULL;	/* mixing array				*/
    float *temp=NULL;	/* temporary array			*/

    int itr, ntr;
    float *trace, *trace2;
    sf_file inp, out, wien;
        

    int n1_headers;
    char* header_format=NULL;
    sf_datatype typehead;
    float* fheader=NULL;
    int verbose;




    /* Initialize */
    sf_init(argc, argv);

    if(!sf_getint("verbose",&verbose))verbose=1;
    /* \n
       flag to control amount of print
       0 terse, 1 informative, 2 chatty, 3 debug
    */
    sf_warning("verbose=%d",verbose);

    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1_traces",&nt))   
      sf_error("input data does not define n1_traces");
        
    if(!sf_histint(inp,"n1_headers",&n1_headers)) 
      sf_error("input data does not define n1_headers");

    header_format=sf_histstring(inp,"header_format");
    if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
    else                                       typehead=SF_FLOAT;

    if(verbose>0)fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
    fheader = sf_floatalloc(n1_headers);

    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    ntr = sf_leftsize(inp,1);

    if (!sf_getfloat("minlag",&minlag)) {
	/* first lag of prediction filter (sec) */
	minlag=0.0f;
	iminlag = 1;
    } else {
	iminlag = SF_NINT(minlag/dt);
    }
    if (!sf_getfloat("maxlag",&maxlag)) {
	/* last lag of prediction filter (sec) */
	maxlag=0.0f;
	imaxlag = SF_NINT(0.05 * nt);
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
    else				  imaxcorr = nt;
    if (imaxcorr > nt) sf_error("maxcorr=%g too large", maxcorr);

    if (imincorr >= imaxcorr)
	sf_error("mincorr=%g, maxcorr=%g", mincorr, maxcorr);
	
    /*... Get mix weighting values values */
    if (!sf_getint("nmix",&nmix)) nmix=1;
    /* number of weights (floats) for moving averages */

    mix = sf_floatalloc(nmix);
    if (!sf_getfloats("mix",mix,nmix)) {
	/* weights for moving average of the autocorrelations */	
	for (imix = 0; imix < nmix; ++imix) 
	    mix[imix] = 1.0;
    } 

    /* Divide mixing weight by number of traces to mix over */
    for (imix = 0; imix < nmix; ++imix)
	mix[imix]=mix[imix]/((float) nmix);

    /* compute filter sizes and correlation number */
    nlag  = imaxlag - iminlag + 1;
    ncorr = imaxcorr - imincorr + 1;
    lcorr = imaxlag + 1;

    /* Compute byte sizes in wiener/spiker and autocorr */
    lagbytes = sizeof(float)*nlag;
    maxlagbytes = sizeof(float)*lcorr;
    mixbytes = maxlagbytes*nmix;

    /* Allocate memory */
    wiener	 = sf_floatalloc(nlag);
    spiker	 = sf_floatalloc(nlag);
    autocorr = sf_floatalloc(lcorr);
    temp = sf_floatalloc(lcorr);
    mixacorr = sf_floatalloc2(lcorr,nmix);

    if (NULL != sf_getstring("wiener")) {
	/* file to output Wiener filter */
	wien = sf_output("wiener");
	sf_putint(wien,"n1",nlag);
    } else {
	wien = NULL;
    }

    /**********************************************************/
    /* end code block for standard tah Trace And Header setup */
    /* continue with any sf_puthist this tah program calls to */
    /* add to the history file                                */
    /**********************************************************/

    /* put the history from the input file to the output */
    sf_fileflush(out,inp);

    /********************************************************/
    /* continue initialization specific to this tah program */
    /********************************************************/





    /* segy_init gets the list header keys required by segykey function  */
    segy_init(n1_headers,inp);

    /* Set pointer to "cross" correlation */
    crosscorr = autocorr + iminlag;

    /* Zero out mixing array */
    memset((void *) mixacorr[0], 0, mixbytes);

    trace  = sf_floatalloc(nt);
    trace2 = sf_floatalloc(nt);

    /* Main loop over traces */
    itr=0;
    while (!(get_tah(trace, fheader, nt, n1_headers, inp))){
	if(verbose>1 || (verbose==1 && itr<5)){
	    fprintf(stderr,"process tah %d in sftahagc\n",itr);
	}

	/********************/
	/* process the tah. */
	/********************/

	/* zero out filter vectors */
	memset((void *) wiener, 0, lagbytes);
	memset((void *) spiker, 0, lagbytes);
	memset((void *) autocorr, 0, maxlagbytes);
	memset((void *) temp, 0, maxlagbytes);

	/* fix supplied by Sanyu Ye */

	xcor(ncorr, 0, trace+imincorr,
	     ncorr, 0, trace+imincorr,
	     lcorr, 0, autocorr);

	/* Leave trace alone if autocorr[0] vanishes */
	if (autocorr[0] == 0.0) {
	    sf_floatwrite(trace,nt,out);
	    if (NULL != wien) sf_floatwrite(wiener,nlag,wien);
	    continue;
	}

	/* Whiten */
	autocorr[0] *= 1.0 + pnoise;

	/* Read autocorr into first column of mixacorr[][] */
	memcpy( (void *) mixacorr[0], 
		(const void *) autocorr, maxlagbytes);

	/* Loop over values of the autocorrelation array */
	for (ilag = 0; ilag < lcorr; ++ilag) {

	    /* Weighted moving average (mix) */
	    for(imix=0; imix<nmix; ++imix)
		temp[ilag]+=mixacorr[imix][ilag]*mix[imix];

	    /* put mixed data back in seismic trace */
	    autocorr[ilag] = temp[ilag]; 

	}

	/* Bump columns of mixacorr[][] over by 1 */
	/* to make space for autocorr from next trace */
	for (imix=nmix-1; 0<imix; --imix)
	    for (ilag=0; ilag<lcorr; ++ilag) 
		mixacorr[imix][ilag] = mixacorr[imix-1][ilag];

	/* Get inverse filter by Wiener-Levinson */
	stoepf(nlag, autocorr, crosscorr, wiener, spiker);
		
	/* Convolve pefilter with trace - don't do zero multiplies */
	for (i = 0; i < nt; ++i) {
	    register int j;
	    register int n = SF_MIN(i, imaxlag); 
	    register float sum = trace[i];

	    for (j = iminlag; j <= n; ++j)
		sum -= wiener[j-iminlag] * trace[i-j];

	    trace2[i] = sum;
	}


	/* Show pefilter on request */
	if (NULL != wien) sf_floatwrite(wiener,nlag,wien);

	put_tah(trace2, fheader, nt, n1_headers, out);
	itr++;
    } 

    exit(0);
}
