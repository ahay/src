/* custom program to blend the fairfield data 

 */
#include <string.h>
#include <rsf.h>

int main(int argc, char* argv[])
{
  int verbose;
  int nt;
  float** in_traces;
  float* delays;
  int ntr;
  int delayfile_ntr;
  float* nodetrace;
  int nt_nodetrace;
  float dt;
  float minshottime, maxshottime;
  int it;
  char* delayfilename; 
  sf_file in=NULL, out=NULL, delayfile=NULL;
  int itrace;

  /* open the input, output, and delay files */
  sf_init (argc,argv);
  in = sf_input ("in");
  out = sf_output ("out");
  
  delayfilename=sf_getstring("delays");
  delayfile = sf_input(delayfilename);
  
  if(!sf_getint("verbose",&verbose))verbose=1;
  /* 0 terse, 1 informative, 2 chatty, 3 debug */
  
  /* get the size of the input traces and the input delayfile */
  sf_histint(in,"n1",&nt);
  sf_histfloat(in,"d1",&dt);
  ntr = sf_leftsize(in,1); /* left dimensions after the first one */
  
  delayfile_ntr = sf_leftsize(delayfile,0); /* number of input points */
  if(ntr!=delayfile_ntr){
    fprintf(stderr,"ntr=%d must equal delayfile_ntr%d",ntr,delayfile_ntr);
    exit(1);
  }
  
  /* allocate space for the input traces and the delays */ 
  in_traces = sf_floatalloc2(nt,ntr);
  delays = sf_floatalloc(ntr);
  
  /* read the input traces and the delays */ 
  sf_floatread(&(in_traces[0][0]),nt*ntr,in);
  sf_floatread(&(delays[0])      ,ntr   ,delayfile);
  if(verbose>2){
    for(itrace=0; itrace<ntr; itrace++){
      fprintf(stderr,"delays[%d]=%f\n",itrace,delays[itrace]); 
    }
  }
    
  /**************************************
   *blend the data into a big node_trace*
   **************************************/
  /* delays array is in ms */
  minshottime=maxshottime=delays[0];
  for(itrace=0; itrace<ntr; itrace++){
    if(minshottime>delays[itrace]) 
      minshottime=delays[itrace];
    if(maxshottime<delays[itrace]) 
      maxshottime=delays[itrace];
  }
  fprintf(stderr,"minshottime=%f, maxshottime=%f\n",minshottime,maxshottime);
  /* the node trace must be long enough to start at t=minshottime and
     extend to t=maxshottime+trace length.
  */
  /* dt is in seconds.  delays array is in ms. */
  nt_nodetrace=(maxshottime-minshottime)/(dt*1000.)+nt;
  nodetrace=sf_floatalloc(nt_nodetrace);
  for(it=0; it<nt_nodetrace; it++) nodetrace[it]=0.0;
  for(itrace=0; itrace<ntr; itrace++){
    int firstindex=(delays[itrace]-minshottime)/(dt*1000.);
    for(it=0; it<nt; it++){
      nodetrace[firstindex+it]+=in_traces[itrace][it];
    }
  }
  /********************************************
   * extract the data from the big node_trace *
   * back in the input trace array            *   
   ********************************************/
  for(itrace=0; itrace<ntr; itrace++){
    int firstindex=(delays[itrace]-minshottime)/(dt*1000.);
    memcpy(in_traces[itrace],&(nodetrace[firstindex]),nt*sizeof(float));
  }
  
  /* write the traces back out */
  sf_floatwrite(&(in_traces[0][0]),nt*ntr,out);
  
  exit(0);
}
