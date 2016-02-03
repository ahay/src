/* Read Trace And Header (tah) from standard input, SPECtral BALance  */

/* tah is the abbreviation of Trace And Header.  Madagascar programs 
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

   The sftahspecbal program is designed to improve the resolution and 
   appearance of the final imaged seismic section (ie after migration 
   and stack.  When applied to noisey land data early in the processing 
   sequence (after surface consisntant decon and before migration and 
   stack) it attenuates noise on post migration and stack data.  There
   are several algorithms called spectral balancing, whitening, or 
   broading.  This program implenents a popular method.  Each input 
   trace is split into several narrow frequency bands by bandpass 
   filtering.  AGC is applied to rach frequency band, and the frequency
   bands are summed.  User parameters control the filter bands and the
   AGC length.  (reference http://www.xsgeo.com/course/spec.htm#content).

   EXAMPLE:

   sftahsort input=shots-receivers-23900_headfix.rsf  \\
   sort="xline:600,601 offset"                        \\
   | sftahspecbal                                     \\
   kls //
   | sftahwrite                                       \\
   verbose=1                                          \\
   mode=seq                                           \\
   output=specbalcmps.rsf                             \\
   >/dev/null

   sfimage <specbalcmps.rsf | sfpen

   In this example the shot organized prestack data in the file 
   shots-receivers-23900_headfix.rsf are read in xline offset order by 
   sftahsort program.  The headers are in the file 
   shots-receivers-23900_headfix_hdr.rsf, the headers parameter default.
   The headers are merged with the trace amplitudes and the tah data sent 
   down the pipe for spectral balancing.  

   kls

   Sftahwrite writes the the trace data to specbalcmp.rsf and the headers 
   are written to the file specbalcmp_hdr.rsf.  The output traces are just
   sequentially written to the output file.
   kls

   PARAMETERS

   float fmin= NULL

   Center of the first frequency band.
	
*/

/*
  Program change history:
  date       Who             What
  09/30/2015 Karl Schleicher Original program.  Derived from Mtahagc.c
*/
#include <string.h>
#include <math.h>
#include <rsf.h>
#include <rsf_su.h>
#include <rsfsegy.h>

#include "tahsub.h"

void agc_apply(float* indata, float* outdata, float* scalars, 
	       int lenagc,int nt,
               int firstlive, int lastlive){
  int half_lenagc;
  int itime;
  double thispower;
  double thisfold;
  float leftamp;
  float rightamp;

  half_lenagc=lenagc/2;
  /* fprintf(stderr,"half_lenagc=%d,lastlive=%d,firstlive=%d\n",
	  half_lenagc   ,lastlive   ,firstlive);
  */
  if(2*half_lenagc+1>lastlive-firstlive+1){
    /* trace is less than half agc. scale by 1 */
    for(itime=0; itime<nt; itime++){
        outdata[itime]=indata[itime];
        scalars[itime]=1.0;
    }
  } else {
    for(itime=firstlive+half_lenagc; itime<lastlive-half_lenagc+1; itime++){
      scalars[itime]=0.0;
    }
    /* sum the first window to get started */
    for(itime=firstlive; itime<firstlive+2*half_lenagc+1; itime++){
      scalars[firstlive+half_lenagc]+=indata[itime]*indata[itime];
    }

    /* in the middle as the window moves right, add the a point on the
       right and subtract an old point on the left. You can save a few 
       multiplications by adding an array, indata2 (the input squared)
    */
    for(itime=firstlive+half_lenagc+1; itime<lastlive-half_lenagc+1; itime++){
      scalars[itime]=scalars[itime-1]-
	indata[itime-half_lenagc-1]*indata[itime-half_lenagc-1]+
	indata[itime+half_lenagc  ]*indata[itime+half_lenagc  ];
    }
    /* scalars contains sum_square_amplitude.  Compute 1/rms */
    for(itime=firstlive+half_lenagc; itime<lastlive-half_lenagc+1; itime++){
      scalars[itime]=1.0/(sqrt(scalars[itime]/(2*half_lenagc+1)));
    }
    /* copy scalars[firstlive+half_lenagc] to beginning and 
       scalars[lastlive-half_lenagc] to the ends of the scalars array */
    for(itime=0; itime<firstlive+half_lenagc; itime++){
      scalars[itime]=scalars[firstlive+half_lenagc];
    }
    for(itime=lastlive-half_lenagc+1; itime<nt; itime++){
      scalars[itime]=scalars[lastlive-half_lenagc];
    }
    /* multiply scalars onto the trace */
    for(itime=0; itime<nt; itime++){
      outdata[itime]=indata[itime]*scalars[itime];
    }
  }
}

void compute_rms(float* indata, float* rms, 
		 int lenagc,int nt,
		 int firstlive, int lastlive){
  int half_lenagc;
  int itime;
  double thispower;
  double thisfold;
  float leftamp;
  float rightamp;

  half_lenagc=lenagc/2;
  /* fprintf(stderr,"half_lenagc=%d,lastlive=%d,firstlive=%d\n",
	  half_lenagc   ,lastlive   ,firstlive);
  */
  if(2*half_lenagc+1>lastlive-firstlive+1){
    /* trace is less than half agc. make rms 1 */
    for(itime=0; itime<nt; itime++){
        rms[itime]=1.0;
    }
    fprintf(stderr,"(2*half_lenagc+1>lastlive-firstlive+1\n");
    fprintf(stderr,"half_lenagc=%d,lastlive=%d,firstlive=%d\n",
	            half_lenagc   ,lastlive   ,firstlive);
  } else {
    for(itime=firstlive+half_lenagc; itime<lastlive-half_lenagc+1; itime++){
      rms[itime]=0.0;
    }
    if(0) {
      /* I tried this rolling summation, but it made nans.  Early trace had 
	 amp like .027. emd of trac had amp like .000005 (no spreading 
	 correction applied).  end of sum_sqr trace had nans. Maybe faster
         algorithm (double precision? leaking integration?) , but I want to 
	 push ahead to look at results. */    
      /* sum the first window to get started */
      for(itime=firstlive; itime<firstlive+2*half_lenagc+1; itime++){
	rms[firstlive+half_lenagc]+=indata[itime]*indata[itime];
      }

      /* in the middle as the window moves right, add the a point on the
	 right and subtract an old point on the left. You can save a few 
	 multiplications by adding an array, indata2 (the input squared)
      */
      for(itime=firstlive+half_lenagc+1; itime<lastlive-half_lenagc+1; itime++){
	rms[itime]=rms[itime-1]-
	  indata[itime-half_lenagc-1]*indata[itime-half_lenagc-1]+
	  indata[itime+half_lenagc  ]*indata[itime+half_lenagc  ];
      }
      if(0){
	fprintf(stderr,"itime indata sum_sqt\n");
	for(itime=1900; itime<nt; itime++){
	  fprintf(stderr,"%d, %e, %e\n",itime,indata[itime],rms[itime]);
	}
      }
    } else { /* here is the slow, dumb, short (stable) algorithm */
      for(itime=firstlive+half_lenagc; itime<=lastlive-half_lenagc; itime++){
	int iitime;
	for (iitime=itime-half_lenagc; iitime<=itime+half_lenagc; iitime++){
	  rms[itime]+=indata[iitime]*indata[iitime];
	}
      }
    }

    /* rms contains sum_square_amplitude.  Compute rms (divide by npt; srqt */
    for(itime=firstlive+half_lenagc; itime<lastlive-half_lenagc+1; itime++){
      rms[itime]=sqrt(rms[itime]/(2*half_lenagc+1));
    }
    /* copy rms[firstlive+half_lenagc] to beginning and 
       rms[lastlive-half_lenagc] to the ends of the rms array */
    for(itime=0; itime<firstlive+half_lenagc; itime++){
      rms[itime]=rms[firstlive+half_lenagc];
    }
    for(itime=lastlive-half_lenagc+1; itime<nt; itime++){
      rms[itime]=rms[lastlive-half_lenagc];
    }
    if(0) { /* kls */
      fprintf(stderr,"itime indata rms\n");
      for(itime=1900; itime<nt; itime++){
	fprintf(stderr,"%d, %e, %e\n",itime,indata[itime],rms[itime]);
      }
    }
  }
}

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
    float* outtrace=NULL;
    float* bandlimittrace=NULL;
    float* rmsall=NULL;
    float* rmsband=NULL;
    float* scalars=NULL;
    float* sum_of_scalars=NULL;
    int indx_time;
    int itrace=0;
    int ntaper;
    int numxstart;
    int numtstart;
    float* taper;
    char **list_of_floats;
    float* xstart;
    float* tstart;
    int indx_of_offset;
    float offset;
    float d1;
    float o1;
    float time_start;
    int itime_start;
    int itime_stop;
    int indx_taper;
    float wagc;
    int lenagc;

    float fmin=5;
    float fmax=95;
    float finc=5;
    int num_centerfreq;
    int firstlive;
    int lastlive;
    float pnoise;

    kiss_fftr_cfg cfg,icfg;
    sf_complex *fft1;
    sf_complex *fft2;
    int nfft;
    int nf;
    float df;
    float scale;
    int ifreq;
    float cntrfreq;
    
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
 
    /*  allocate intrace and outtrace after initialize fft and 
	make them padded */
  
    /* maybe I should add some validation that n1== n1_traces+n1_headers+2
       and the record length read in the second word is consistent with 
       n1.  */

    /**********************************************************/
    /* end code block for standard tah Trace And Header setup */
    /* continue with any sf_puthist this tah program calls to */
    /* add to the history file                                */
    /**********************************************************/

    /* put the history from the input file to the output */
    sf_fileflush(out,in);

    /********************************************************/
    /* continue initialization specific to this tah program */
    /********************************************************/





    /* segy_init gets the list header keys required by segykey function  */
    segy_init(n1_headers,in);

    if(0==1) {  
      /* infuture may want to have offset dependent agc design start time */
      /* get the mute parameters */
      if(NULL==(list_of_floats=sf_getnstring("xstart",&numxstart))){
	xstart=NULL;
	sf_error("xstart is a required parameter in sftahagc");
      } else {
	xstart=sf_floatalloc(numxstart);
	if(!sf_getfloats("xstart",xstart,numxstart))sf_error("unable to read xstart");
      }
      if(NULL==(list_of_floats=sf_getnstring("tstart",&numtstart))){
	tstart=NULL;
	sf_error("xstart is a required parameter in sftahagc");
      } else {
	tstart=sf_floatalloc(numtstart);
	if(!sf_getfloats("tstart",tstart,numtstart))sf_error("unable to read tstart");
      }
      if(numxstart!=numtstart)sf_error("bad mute parameters: numxstart!=numtstart");
      if(!sf_getint("ntaper",&ntaper))ntaper=12;
    }
    /* \n
       length of the taper on the stack mute
    */
    /* is this of use for agc?
    taper=sf_floatalloc(ntaper);
    for(indx_time=0; indx_time<ntaper; indx_time++){
	float val_sin=sin((indx_time+1)*SF_PI/(2*ntaper));
	taper[indx_time]=val_sin*val_sin;
    }
    */
    indx_of_offset=segykey("offset");
    
    if (!sf_histfloat(in,"d1",&d1))
      sf_error("input data does not define d1");
    if (!sf_histfloat(in,"o1",&o1))
      sf_error("input data does not define o1");

    /* check wagc for reasonableness.  I input 250 (not .25 and ,250) */  
    if(!sf_getfloat("wagc",&wagc)){
      sf_error("wagc is a required parameter in sftahagc");
    }
    /* \n
       length of the agc window in seconds
    */
    lenagc=wagc/d1+1.5; /* length of the agc window in samples */
 
    if (!sf_getfloat("pnoise",  &pnoise)) pnoise = 0.01;
    /* relative additive noise level */
    
    if (!sf_getfloat("fmin",&fmin))fmin=5;
    /* minimum frequency band */
    if (!sf_getfloat("fmax",&fmax))fmax=95;
    /* maximum frequency band */
    if (!sf_getfloat("finc",&finc))finc=5;
    /* frequency band increment */

    /* initialize kiss_fft kls */

    nfft  = kiss_fft_next_fast_size(n1_traces+200);
    nf    = nfft/2+1;
    df    = 1./(nfft*d1);
    scale = 1./nfft;
    cfg   = kiss_fftr_alloc(nfft,0,NULL,NULL);
    icfg  = kiss_fftr_alloc(nfft,1,NULL,NULL);
    fft1  = sf_complexalloc(nf);
    fft2  = sf_complexalloc(nf);
    if(verbose>0)fprintf(stderr,"fmax=%f\n",(nf-1)*df);

    /* allocate input and output traces with enough padding for
       ffts */
    if(verbose>0)
      fprintf(stderr,"allocate padded intrace.  n1_traces=%d nfft=%d\n",
	                                        n1_traces,   nfft);
    intrace       =sf_floatalloc(nfft); /* must be padded for input to fft */
    bandlimittrace=sf_floatalloc(nfft); /* must be padded for input to fft */
    outtrace= sf_floatalloc(n1_traces);
    rmsall        =sf_floatalloc(n1_traces);
    rmsband        =sf_floatalloc(n1_traces);
    scalars       =sf_floatalloc(n1_traces);
    sum_of_scalars=sf_floatalloc(n1_traces);
    num_centerfreq=(fmax-fmin)/finc+1.95;
    if(fabs(fmax-(fmin+(num_centerfreq-1)*finc))>.01){
      fprintf(stderr,"*************************************\n");
      fprintf(stderr,"*************************************\n");
      fprintf(stderr,"*************************************\n");
      fprintf(stderr,"fmin-fmax is not a multiple of finc\n");
      fprintf(stderr,"fmin=%f, fmax=%f, finc=%f\n",fmin,fmax,finc);
      fprintf(stderr,"*************************************\n");
      fprintf(stderr,"*************************************\n");
      sf_error("illegal combination of fmin,fmax,fminc");
    }

    /***************************/
    /* start trace loop        */
    /***************************/
    if(verbose>0)fprintf(stderr,"start trace loop\n");
 
    itrace=0;
    while (!(get_tah(intrace, fheader, n1_traces, n1_headers, in))){
	if(verbose>1 || (verbose==1 && itrace<5)){
	    fprintf(stderr,"process tah %d in sftahagc\n",itrace);
	}
	/********************/
	/* process the tah. */
	/********************/

	if(typehead == SF_INT)offset=((int  *)fheader)[indx_of_offset];
	else                  offset=((float*)fheader)[indx_of_offset];
	/* maybe latter add agc design start time
	intlin(numxstart,xstart,tstart,
	       tstart[0],tstart[numxstart-1],1,
	       &offset,&time_start);
	if(time_start<o1)time_start=o1;

	itime_start=(int)(((time_start-o1)/d1)+.5);
	*/
	/* find firstlive and lastlive */
        if(verbose>2)fprintf(stderr,"find firstlive and lastlive\n");
	for (firstlive=0; firstlive<n1_traces; firstlive++){
	  if(intrace[firstlive] != 0) break;
	}
	for (lastlive=n1_traces-1;lastlive>=0; lastlive--){
	  if(intrace[lastlive] != 0) break;
	}
	/* kls need to catch a zero trace */
        if(verbose>2)
	  fprintf(stderr,"firstlive=%d, lastlive=%d\n",firstlive,lastlive);
	/* zero the padded area on the input trace */
	for (indx_time=n1_traces; indx_time<nfft; indx_time++){
	  intrace[indx_time]=0;
	}
	
	/********************************************************/
	/* apply the SPECtral BALancing:                        */
	/* 1 compute the input trace rms in agclen gates        */
	/* 2 fft to frequency domain                            */
	/* 3 for each frequency band                            */
	/*   3.1 extract (ramp) the frequency band              */
	/*   3.2 convert freqency band to time,                 */
	/*   3.3 compute the rms of the frequency band          */
	/*   3.4 agc scalars 1/(freq_band_rms + whitenoise*rms) */ 
	/*   3.3 agc the frequency band                         */
	/*   3.4 sum to output                                  */
        /*   3.5 sum scalars to sum of scalars                  */
	/* 4 divide by the sum of scalars                       */
	/********************************************************/
	compute_rms(intrace, rmsall,
		    lenagc, n1_traces,
                    firstlive, lastlive);
	/* scale rms by whitenoise factor */
	if(verbose>2)fprintf(stderr,"put_tah rmsall\n");
	if(0) put_tah(intrace, fheader, n1_traces, n1_headers, out);
        if(0) put_tah(rmsall, fheader, n1_traces, n1_headers, out);
	for (indx_time=0; indx_time<n1_traces; indx_time++){
	  rmsall[indx_time]*=(1.0/num_centerfreq);
	}

	if(verbose>2)fprintf(stderr,"fftr\n");
	kiss_fftr(cfg, intrace, (kiss_fft_cpx*) fft1);

	/* zero the output trace and the sum_of_scalars */
	for (indx_time=0; indx_time<n1_traces; indx_time++){
	  outtrace[indx_time]=0.0;
	  sum_of_scalars[indx_time]=0.0;
	}
	for (cntrfreq=fmin; cntrfreq<=fmax; cntrfreq+=finc){
	  if(verbose>2)fprintf(stderr,"cntrfreq=%f\n",cntrfreq);
          /* zero frequencies before the band */
	  for (ifreq=0; ifreq<(int)((cntrfreq-finc)/df); ifreq++){
	    fft2[ifreq]=0.0;
	  }
	  /* triangular weight the selected  frequency band */
	  for (ifreq=(int)((cntrfreq-finc)/df); 
	       ifreq<(int)((cntrfreq+finc)/df) && ifreq<nf;
	       ifreq++){
            float weight;
	    if(ifreq>0){
	      weight=(1.0-fabs(ifreq*df-cntrfreq)/finc)/nfft;
	      fft2[ifreq]=weight*fft1[ifreq];
	    }
	  }
          /* zero frequencies larger than the band */
	  for (ifreq=(int)((cntrfreq+finc)/df); ifreq<nf; ifreq++){
	    fft2[ifreq]=0.0;
	  } 
	  /* inverse fft back to time domain */
          if(verbose>2)fprintf(stderr,"fftri\n");
	  kiss_fftri(icfg,(kiss_fft_cpx*) fft2, bandlimittrace);
          /* apply agc to the bandlimittrace and sum scalars to 
             sum_of_scalars */
	  if(verbose>2)fprintf(stderr,"agc_apply\n");
	  compute_rms(bandlimittrace, rmsband,
		    lenagc, n1_traces,
                    firstlive, lastlive);
	  if(verbose>2)fprintf(stderr,"sum to ouput\n");
	  for (indx_time=0; indx_time<n1_traces; indx_time++){
	    /* old poor 
	       scalars[indx_time]=1.0/
                   (rmsband[indx_time]+pnoise*rmsall[indx_time]);
	    */
	    if(1)scalars[indx_time]=rmsband[indx_time]/
                     (rmsband[indx_time]*rmsband[indx_time]+
		      pnoise*rmsall[indx_time]*rmsall[indx_time]);
	    else scalars[indx_time]=1.0;
	  }
	  for (indx_time=0; indx_time<n1_traces; indx_time++){
	    bandlimittrace[indx_time]*=scalars[indx_time];
	    outtrace[indx_time]+=bandlimittrace[indx_time];
	  }
	  for (indx_time=0; indx_time<n1_traces; indx_time++){
	    sum_of_scalars[indx_time]+=scalars[indx_time];
	  }
          if(0)
	    put_tah(bandlimittrace, fheader, n1_traces, n1_headers, out);
	  if(0)
	    put_tah(rmsband, fheader, n1_traces, n1_headers, out);
	  if(0)
	    put_tah(scalars, fheader, n1_traces, n1_headers, out);
	}
        if(0)put_tah(sum_of_scalars, fheader, n1_traces, n1_headers, out);
        if(1){	
	   /* divide output by sum of scalars */
	  for (indx_time=0; indx_time<n1_traces; indx_time++){
	    outtrace[indx_time]/=sum_of_scalars[indx_time];
	  }
	} 
	if (1)
	  put_tah(outtrace, fheader, n1_traces, n1_headers, out);
	itrace++;
    }

    exit(0);
}

  
