/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/* SUFAIRFIELDBLEND: $Revision: 1.36 $ ; $Date: 2010/02/02 22:11:31 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" sufairfieldblend 	   								",
" SUFAIRFIELDBLEND - experiment with blending the fairfield    	        ", 
" 	   								",
" sufairfieldblend <stdin >stdout ntr=1120                              ",
" 									",
" parameters:							        ",
" ntr=1120    		number of input traces                          ",
"									",
" Notes: 								",
"									",
"  									",
NULL};

/* Credits:
 *
 *	karl Schleicher : derived initial program from sumute 
 */
/**************** end self doc ***********************************/

segy tr;


int
main(int argc, char **argv)
{
  segy tr;
  float** in_traces;
  segy** in_headers;
  int nt;
  int ntr;
  int itrace;
  int minshottime, maxshottime;
  int it;
  int nt_nodetrace;
  float* nodetrace;
  float dt;

  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);
  
  
  /******************
   * Get parameters *
   ******************/
  fprintf(stderr,"start reading parameters\n");
  if(!getparint("ntr", &ntr))
    ntr=560*2; /* I got the number traces from surange */
  
  /*****************************
   * Get info from first trace *
   *****************************/
  if (!gettr(&tr)) err("can't read first trace");
 
  nt=tr.ns;
  dt=tr.dt;
  in_traces=alloc2float(nt,ntr);
  in_headers=(segy**)alloc2int(60,ntr);


  
  /********************/
  /* Loop over traces */
  /********************/
  fprintf(stderr,"start the trace loop\n");
  
  itrace=0;
  do {
    memcpy(in_traces[itrace],tr.data,nt*sizeof(float));
    memcpy(in_headers[itrace],&tr,240);
    itrace++;

  } while (gettr(&tr));

  ntr=itrace;

  /**************************************
   *blend the data into a big node_trace*
   **************************************/
  
  minshottime=maxshottime=in_headers[0]->gelev;
  for(itrace=0; itrace<ntr; itrace++){
    segy* my_hdr;
    my_hdr=in_headers[itrace];
    fprintf(stderr,"%d %d %d %d %d %d %d %d %d %d \n",
	    my_hdr->tracr,my_hdr->fldr,my_hdr->sx,my_hdr->sy,
	    my_hdr->gx,my_hdr->gy,my_hdr->ns,my_hdr->dt,
	    my_hdr->gelev,my_hdr->selev);
    /* use gelev to blend shots 754 and 763.  This was originally in header 
       bytes 41-44.  These are the time delays for synchronized blending
    */
    if(minshottime>in_headers[itrace]->gelev) 
      minshottime=in_headers[itrace]->gelev;
    if(maxshottime<in_headers[itrace]->gelev) 
      maxshottime=in_headers[itrace]->gelev;
  }
  /* the node trace must be long enough to start at t=minshottime and
     extend to t=maxshottime+trace length.
  */
  nt_nodetrace=(maxshottime-minshottime)/(dt/1000.)+nt;
  nodetrace=alloc1float(nt_nodetrace);
  for(it=0; it<nt_nodetrace; it++) nodetrace[it]=0.0;
  for(itrace=0; itrace<ntr; itrace++){
    int firstindex=(in_headers[itrace]->gelev-minshottime)/(dt/1000.);
    for(it=0; it<nt; it++){
      nodetrace[firstindex+it]+=in_traces[itrace][it];
    }

  }
  /********************************************
   * extract the data from the big node_trace *
   ********************************************/
  for(itrace=0; itrace<ntr; itrace++){
    int firstindex=(in_headers[itrace]->gelev-minshottime)/(dt/1000.);
    memcpy(tr.data,&(nodetrace[firstindex]),nt*sizeof(float));
    memcpy(&tr    ,in_headers[itrace],240);
    puttr(&tr);
  }

  return(CWP_Exit());
}
