/* Seismic blending operator.
 Custom program to blend the seismic data.
*/
/*
   Copyright (C) 2013 University of Texas at Austin - Karl Schleicher
     
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

int main(int argc, char* argv[])
{
  int verbose;
  int nt;
  float** in_traces;
  float** out_traces;
  float* shot_time_in;
  float* shot_time_out;
  int ntr;
  int shot_time_in_file_ntr;
  int shot_time_out_file_ntr;
  float* nodetrace;
  int nt_nodetrace;
  float dt;
  float minshottime, maxshottime;
  int it;
  char* shot_time_in_filename; 
  char* shot_time_out_filename; 
  sf_file in=NULL, out=NULL, shot_time_in_file=NULL, shot_time_out_file=NULL;
  int itrace;

  /* open the input, output, and delay files */
  sf_init (argc,argv);
  in = sf_input ("in");
  out = sf_output ("out");
  
  shot_time_in_filename=sf_getstring("shot_time_in");
  shot_time_in_file = sf_input(shot_time_in_filename);

  shot_time_out_filename=sf_getstring("shot_time_out");
  shot_time_out_file = sf_input(shot_time_out_filename);
  
  if(!sf_getint("verbose",&verbose))verbose=1;
  /* 0 terse, 1 informative, 2 chatty, 3 debug */
  
  /* get the size of the input traces and the input shot_time_in_file */
  sf_histint(in,"n1",&nt);
  sf_histfloat(in,"d1",&dt);
  ntr = sf_leftsize(in,1); /* left dimensions after the first one */
  
  shot_time_in_file_ntr = sf_leftsize(shot_time_in_file,0); /* number of input points */
  if(ntr!=shot_time_in_file_ntr){
    fprintf(stderr,"ntr=%d must equal shot_time_in_file_ntr%d",ntr,shot_time_in_file_ntr);
    exit(1);
  }
  
  shot_time_out_file_ntr = sf_leftsize(shot_time_out_file,0); /* number of output points */
  /* allocate space for the input traces and the shot_time_in */ 
  in_traces = sf_floatalloc2(nt,ntr);
  shot_time_in = sf_floatalloc(ntr);
  out_traces = sf_floatalloc2(nt,shot_time_out_file_ntr);
  shot_time_out = sf_floatalloc(shot_time_out_file_ntr);
  
  /* read the input traces and the shot_time_in */ 
  sf_floatread(&(in_traces[0][0]),nt*ntr,in);
  sf_floatread(&(shot_time_in[0])      ,ntr   ,shot_time_in_file);
  if(verbose>2){
    for(itrace=0; itrace<ntr; itrace++){
      fprintf(stderr,"shot_time_in[%d]=%f\n",itrace,shot_time_in[itrace]); 
    }
  }
  sf_floatread(&(shot_time_out[0])      ,shot_time_out_file_ntr,
	       shot_time_out_file);
    
  /**************************************
   *blend the data into a big node_trace*
   **************************************/
  /* shot_time_in array is in ms */
  minshottime=maxshottime=shot_time_in[0];
  for(itrace=0; itrace<ntr; itrace++){
    if(minshottime>shot_time_in[itrace]) 
      minshottime=shot_time_in[itrace];
    if(maxshottime<shot_time_in[itrace]) 
      maxshottime=shot_time_in[itrace];
  }
  for(itrace=0; itrace<shot_time_out_file_ntr; itrace++){
    if(minshottime>shot_time_out[itrace]) 
      minshottime=shot_time_out[itrace];
    if(maxshottime<shot_time_out[itrace]) 
      maxshottime=shot_time_out[itrace];
  }
  fprintf(stderr,"minshottime=%f, maxshottime=%f\n",minshottime,maxshottime);
  /* the node trace must be long enough to start at t=minshottime and
     extend to t=maxshottime+trace length.
  */
  /* dt is in seconds.  shot_time_in array is in ms. */
  nt_nodetrace=(maxshottime-minshottime)/(dt*1000.)+nt;
  nodetrace=sf_floatalloc(nt_nodetrace);
  for(it=0; it<nt_nodetrace; it++) nodetrace[it]=0.0;
  for(itrace=0; itrace<ntr; itrace++){
    int firstindex=(shot_time_in[itrace]-minshottime) /(dt*1000.);
    for(it=0; it<nt; it++){
      nodetrace[firstindex+it]+=in_traces[itrace][it];
    }
  }
  /********************************************
   * extract the data from the big node_trace *
   * back in the input trace array            *   
   ********************************************/
  for(itrace=0; itrace<shot_time_out_file_ntr; itrace++){
    int firstindex=(shot_time_out[itrace]-minshottime)/(dt*1000.);
    memcpy(out_traces[itrace],&(nodetrace[firstindex]),nt*sizeof(float));
  }
  
  /* write the traces back out */
  sf_floatwrite(&(out_traces[0][0]),nt*shot_time_out_file_ntr,out);
  
  exit(0);
}
