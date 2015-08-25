/* Read Trace And Header (tah) from standard input, MUTE */

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

   The sftahmute program is designed to mute data. Trace and header data 
   (tah) are read from standard input (usually a pipe).  The trace xstart 
   and tmute parameter define the mute start time.  They are interpolated 
   to determine the start time for the trace using the trace header 
   offset.  The ntaper defines the length in samples of the taper.

   EXAMPLE:

   sftahsort input=shots-receivers-23900_headfix.rsf           \\
   sort="xline:600,601 offset"                              \\
   | sftahnmo tnmo=0,2,6,10.5,16 vnmo=1500,1500,2250,3250,3700 \\
   | sftahmute                                                 \\
   xstart=0,20000 tmute=0,20 ntaper=25                        \\
   | sftahnmo                                                  \\
   tnmo=0,2,6,10.5,16                                        \\
   vnmo=1500,1500,2250,3250,3700                             \\
   inv=y                                                     \\
   | sftahmakeskey pkey=xline skey=cdpt                        \\
   | sftahwrite                                                \\
   verbose=1                                                 \\
   label2=cdpt  o2=1 n2=100 d2=1                             \\
   label3=xline o3=600 n3=1 d3=1                             \\
   output=mutecmps.rsf                                       \\
   >/dev/null

   sfgrey <mutecmps.rsf | sfpen

   In this example the shot organized prestack data in the file 
   shots-receivers-23900_headfix.rsf are read in xline offset order by 
   sftahsort program.  The headers are in the file 
   shots-receivers-23900_headfix_hdr.rsf, the headers parameter default.
   The headers are merged with the trace amplitudes and the tah data sent 
   down the pipe for nmo, mute, and inverse nmo.  This sequence was used 
   to apply the mute using times that were selected from a prestack 
   gather with moveout applied.

   The sftahnmo program uses the velocity function defined by the tnmo, 
   vnmo parameters and the offset header to apply normal moveout to 
   the traces.  

   sftahmute zeros the shallow data.  TLhe time samples above the line 
   through (time,offset) pairs (0,0)(20,20000), are set to zero. There 
   is a 25 point taper applied below the zero portion of the traces.

   A second sftahnmo execution applied inverse nmoout.  Other than inv=yes 
   the parameters are the same as in the first sftahnmo. 

   The program sftahmakeskey is used to create a secondary key used 
   in the following sftahwrite to define the location to wrte the trace 
   in the output file. Sftahmakeskey makes a secondary key (skey=cdpt) 
   the count the traces starting in the a primary key gather (pkey=xline).
   The input traces gathered by xline by sftahsort. Sftahmakeskey sets 
   cdpt to 1 when the trace has a new xline.  If the trace has the same 
   xline as the previous trace cdpt is incremented

   Sftahwrite writes the the trace data to mutecmp.rsf and the headers are 
   written to the file mutecmp_hdr.rsf.  The order of the data in the output 
   file is defined by the cdpt and xline trace headers, so the  data order
   is (time,cmpt,xline).  Finally, the output volume is displayed using
   sfgrey.

   PARAMETERS
   strings key= no default

   list of header keys to monitor to determine when to break 
   between gathers.  A gather is a sequence of traces with the 
   same value for all the header keys.  Stack summs traces in 
   the gather, divides by the fold, and outputs the stack trace.

   floats xstart= NULL

   List of floats the same length as list of floats in the tstart
   parameter.  The (xstart,tstart) pairs are interpolated using the
   trace headers offset to determine trace start time.  The mute is
   NOT moved based on the first live sample.

   floats tstart= NULL

   List of floats the same length as list of floats in the xstart
   parameter.  The (xstart,tstart) pairs are interpolated using the
   trace headers offset to determine trace start time. The mute is
   NOT moved based on the first live sample.

   float ntaper=12
   the length of the taper to use at the start of the trace.
	
*/

/*
  Program change history:
  date       Who             What
  12/1/2014 Karl Schleicher Original program.  Derived from Mtahstack.c
*/
#include <string.h>
#include <rsf.h>
#include <rsf_su.h>
#include <rsfsegy.h>

#include "tahsub.h"

typedef struct Agc_struct *agc_struct;
struct Agc_struct{
  float* agc_scalar;
  float* agc_power;
  float* agc_fold;
  int agc_nt;
};

agc_struct agc_init(int nt){
  agc_struct agc_structure;
  agc_structure=(agc_struct) sf_alloc(1,sizeof(*agc_structure));
  agc_structure->agc_scalar=(float*)sf_alloc(nt,sizeof(float));
  agc_structure->agc_power=(float*)sf_alloc(nt,sizeof(float));
  agc_structure->agc_fold=(float*)sf_alloc(nt,sizeof(float));
  agc_structure->agc_nt=nt;
  return agc_structure;
}

void agc_close(agc_struct agc_structure){
  free(agc_structure->agc_scalar);
  free(agc_structure->agc_power );
  free(agc_structure->agc_fold  );
  free(agc_structure);
}

void agc_apply(agc_struct agc_structure,float* indata,float* outdata,
	       int lenagc,int nt){
  int half_lenagc=lenagc/2;
  int itime;
  double thispower;
  double thisfold;
  float leftamp;
  float rightamp;

  if(nt != agc_structure->agc_nt){
    sf_error("agc_apply called with nt=%d, but struct->agc_nt=%d",
	                            nt,agc_structure->agc_nt   );
  }

  if(half_lenagc>nt-1){
    /* trace is less than half agc.  zero trace. */
    for(itime=0; itime<nt; itime++){
        outdata[itime]=0.0;
    }
  } else {
    thisfold=0;
    thispower=0;
    /* sum the right half window to get started */
    for(itime=0; itime<half_lenagc; itime++){
      thispower+=indata[itime]*indata[itime];
      thisfold++;
    }
    for(itime=0; itime<nt; itime++){
      /* in the middle as the window moves right, add the a point on the
	 right and subtract an old point on the left.  Test so points
	 off the vecrot are not added/subtracted.  Keep track of the fold. 
      */
      if(itime+half_lenagc<nt){
	rightamp=indata[itime+half_lenagc];
	if(rightamp!=0.0){
	  thispower+=rightamp*rightamp;
	  thisfold++;
	}
      }
      if(itime-half_lenagc-1>=0){
	leftamp=indata[itime-half_lenagc];
	if(leftamp!=0.0){
	  thispower-=leftamp*leftamp;
	  thisfold--;
	}
      }
      agc_structure->agc_power[itime]=thispower;
      agc_structure->agc_fold[itime]=thisfold;
      if (thisfold>=half_lenagc){
	agc_structure->agc_scalar[itime]=thisfold/thispower;
	outdata[itime]=sqrt(thisfold/thispower)*indata[itime];
      } else {
	/* fprintf(stderr,"itime=%d,thisfold=%f,half_lenagc=%d\n",
	   itime,thisfold,  half_lenagc); */
	agc_structure->agc_scalar[itime]=0.0;
	outdata[itime]=0.0;
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
    agc_struct agc_structure;
    int lenagc=500;

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
 
    if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
    intrace = sf_floatalloc(n1_traces);
    outtrace= sf_floatalloc(n1_traces);

  
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
    indx_of_offset=segykey("offset");
    */
    if (!sf_histfloat(in,"d1",&d1))
	sf_error("input data does not define d1");
    if (!sf_histfloat(in,"o1",&o1))
	sf_error("input data does not define o1");

    /* initialize the agc structures */
    agc_structure=agc_init(n1_traces);
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

	/* apply agc to trace */
	agc_apply(agc_structure,intrace,outtrace,
		  lenagc,n1_traces);

	put_tah(outtrace, fheader, n1_traces, n1_headers, out);
	itrace++;
    }

    exit(0);
}

  
