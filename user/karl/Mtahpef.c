/* Trace And Header Prediction Error Filtering 

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

The sftahpef applies prediction error filtering (often called decon)

EXAMPLE:

sftahsort                                                            //
verbose=0 input=npr3_field.rsf sort='fldr:214,254,10 tracf'          //
| sftahwindow ns=2047                                                //
| sftahgain   tpow=2                                                 //
| sftahmute  tmute=-.200,-.050,.200,3.00  xmute=0,880,1760,18000     //
ntaper=80                                                            //
| sftahpef                                                           //
verbose=1 minlag=.002 maxlag=.1  pnoise=.01 mincorr=0 maxcorr=3      //
| sftahagc  wagc=1.000 verbose=1                                     //
| sftahwrite verbose=1 mode=seq  output=pefshots.rsf                 //
>/dev/null

sfgrey <pefshots.rsf | sfpen

In this example unprocessed field traces are read by sftahsort from 
the npr3_field.rsf file.  sftahsort was used select just a 5 shotpoints 
(fldr 214 to 254 incrementing by 10) from a large dataset.  The headers
are in the file npr3_filed_hdr.rsf, the headers parameter default.  
The headers are merged with the trace amplitudes and the tah data sent 
down the pipe for further processing.

sftahwindow selects the first 2047 trace amplitudes.  The last two 
samples on this data were bad, and are elliminated from further
processing.
  
sftahgain multiplies the traces by time squared (t**2).  This 
approximately compensates for spreading loss that makes amplitude at
large time smaller than amplitude at small time.

sftahmute is applied to elliminate the data at small time/offset.  

sftahpef applies prediction error filtering (or decon).  A prediction 
gap or .002 seconds, or one sample) makes this decon "spiking" decon.
A three seconds window si selected to compute the autocorrelation and 
a .1 second filter is computed and applied.

After decon a 1 second agc was applies using the sftahagc.

Data is written to the output file, pefshots.rsf, using sftahwrite.  
Traces headers are written to pefshots_hdr.rsf.  The output file data 
order is sequential, or just in the order sftahwrite reads them from 
standard input.

 
*/

/* derived from Seismic Unix supef  */

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

    int itr; /* ntr; */
    float *trace, *trace2;
    sf_file inp, out, wien;
        

    int n1_headers;
/*    char* header_format=NULL; */
/*    sf_datatype typehead; */
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

/*    header_format=sf_histstring(inp,"header_format");
      if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
      else                                       typehead=SF_FLOAT; */

    if(verbose>0)fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
    fheader = sf_floatalloc(n1_headers);

    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    /*   ntr = sf_leftsize(inp,1); */

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
