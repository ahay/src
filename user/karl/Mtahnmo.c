/* tahnmo: Trace And Header Normal MoveOut.

   tah is the abbreviation of Trace And Header.  It identifies a group of
   programs designed to:
   1- read trace and headers from separate rsf files and write them to 
      standard output
   2- filter programs that read and write standard input/output and process 
      the tah data
   3- read tah data from standard input and write separate rsf files for the
      trace and headers data

   These programs allow Seismic Unix (su) like processing in Madagascar.  
   Some programs have su like names.

   Some programs in this suite are sf_tahread, sf_tahgethw, f_tahhdrmath, 
   and sf_tahwrite.
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
   11/14/2013 Karl Schleicher Original program based on Mtahgethw.c
*/
#include <string.h>
#include <rsf.h>
#include <rsfsegy.h>
#include <math.h>

#include "tahsub.h"
#include "mksinc.h"
#include "stoep.h"
#include "intsinc8.h"
#include "inttable8.h"
#include "intlin.h"
#include "yxtoxy.h"

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
  float o1;
  float d1;
  int indx_offset;
  float* local_sloth=NULL;
  float v0;
  int indx_time;
  float offset;
  float offset2;
  float* r_index_tx_of_it0=NULL;
  float* r_index_t0_of_itx=NULL;
  int start_indx_nmo;
  int start_indx_nmo_tx;
  float* outtrace=NULL;
  float itrace=0;
  float nmostretch;
  float lmute;
  bool inv;
  /******************************************************/
  /* code block for standard tah Trace And Header setup */
  /******************************************************/

  sf_init (argc,argv);

  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  /* verbose flag controls ammount of print */
  /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
  /* fprintf(stderr,"read verbose switch.  getint reads command line.\n"); */
  if(!sf_getint("verbose",&verbose))verbose=1;
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
  intrace= sf_floatalloc(n1_traces);

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
  if (!sf_histfloat(in,"d1",&d1))
    sf_error("input data not define d1");
  if (!sf_histfloat(in,"o1",&o1))
    sf_error("input data not define o1");
  /* kls should read label1 and verify it is time  
  if (!sf_histstring(in,"label1",&label1))
    sf_error("input data not define label1");
  */

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);
  /* get index to keys I will be using */
  indx_offset=segykey("offset");
  /* kls what other header keys do I use?  inline? xline? cdp? */

  /* get the parameter for the maximum nmo stretch. */
  if (!sf_getfloat("str",&nmostretch)) nmostretch=0.5;
  if (!sf_getfloat("lmute",&lmute)) lmute=12.*d1;
  /* length of the mute zone in seconds */
  lmute/=d1;
  fprintf(stderr,"lmute=%f\n",lmute);
  if(!sf_getbool("inv",&inv)) inv=false;
  /* if y, do inverse stretching */
  if(inv)fprintf(stderr,"inv=true\n");
  else fprintf(stderr,"inv=false\n");
    /* if y, do inverse stretching */

  /* set up velocity function ( really (1/v)**2, sloth */
  local_sloth=sf_floatalloc(n1_traces);
  /* just constant velocity today */
  if(1==0){
    sf_getfloat("v0",&v0);
    for(indx_time=0; indx_time<n1_traces; indx_time++){
      local_sloth[indx_time]=1.0/(v0*v0);
    }
  } else {
    char** list_of_floats;
    float* vnmo;
    float* tnmo;
    int numvnmo;
    int numtnmo;
    float t0;

    if(verbose>1)fprintf(stderr,"read vnmo/tnmo\n");
    /* use this fundtion to find out number of velocities and time 
       input in vnmo and tnmo */
    list_of_floats=sf_getnstring("vnmo",&numvnmo);
    if(verbose>1){
      int i;
      fprintf(stderr,"numvnmo=%d\n",numvnmo);
      for (i=0; i<numvnmo; i++){
	fprintf(stderr,"list_of_float=%s\n",list_of_floats[i]);
      }
    }
    /* should free this list of strings, but that is only a little memory */\
    list_of_floats=sf_getnstring("tnmo",&numtnmo);
    if(verbose>1){
      int i;
      for (i=0; i<numtnmo; i++){
	fprintf(stderr,"list_of_float=%s\n",list_of_floats[i]);
      }
    }
    if(numvnmo!=numtnmo){
      sf_error("number vnmo floats=%d != number tnmo floats=%d",
	       numvnmo,numtnmo);
    }
    if(numvnmo==0)sf_error("vnmo parameter is required");
    vnmo=sf_floatalloc(numvnmo);
    tnmo=sf_floatalloc(numtnmo);
    if(verbose>1)fprintf(stderr,"sf_getfloats(vnmo)");
    if (!sf_getfloats("vnmo",vnmo,numvnmo))
      sf_error("unable to read vnmo");
    if(verbose>1)fprintf(stderr,"sf_getfloats(tnmo)");
    if (!sf_getfloats("tnmo",tnmo,numtnmo))
      sf_error("unable to read tnmo");
    if(verbose>1){
      for(indx_time=0; indx_time<numvnmo; indx_time++){
	fprintf(stderr,"indx=%d, vnmo=%f, tnmo=%f\n",
		indx_time,vnmo[indx_time],tnmo[indx_time]);
      }
    }
    if(verbose>1)fprintf(stderr,"interpolate the velocity\n");
    for(indx_time=0; indx_time<n1_traces; indx_time++){
      t0=indx_time*d1+o1;
      intlin(numtnmo,tnmo,vnmo,vnmo[0],vnmo[numvnmo-1],1,&t0,&v0);
      local_sloth[indx_time]=1.0/(v0*v0);
    }    
  }
    
  if(verbose>0)fprintf(stderr,"allocate arrays for the trace loop\n");
  r_index_tx_of_it0=sf_floatalloc(n1_traces);
  r_index_t0_of_itx=sf_floatalloc(n1_traces);
  outtrace         =sf_floatalloc(n1_traces);
  

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
  while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
    if(verbose>1)fprintf(stderr,"process the tah in sftahnmo\n");
    /********************/
    /* process the tah. */
    /********************/
    /* this program applies moveout */
    /* kls this should be only be done when velocity (local_sloth) 
       or offset changes */
    if(typehead == SF_INT){   
      /* just cast the header to int so the print works */
      offset=((int*)fheader)[indx_offset];
    } else {
      offset=       fheader [indx_offset];
    }
    offset2=offset*offset;
    for(indx_time=0; indx_time<n1_traces; indx_time++){
      float tx, t0;
      t0=indx_time*d1+o1;
      tx=sqrt(t0*t0+offset*offset*local_sloth[indx_time]);
      r_index_tx_of_it0[indx_time]=(tx-o1)/d1;
      if(itrace==0 && verbose>4){
	fprintf(stderr,"indx_time=%d, tx=%f, sloth=%g, offset=%f, t0=%f\n",
		        indx_time   , tx   , local_sloth[indx_time]  ,
                                                       offset   , t0);
      }
    }
    /* kls nmo start time should depend on the nmo stretch limit.
       Find the excessive stretch closest to the bottom of the trace.  This 
       is the last time nmstrewtch is violated.  It is OK to apply nmo
       to the rest of the trace. */
    for (start_indx_nmo=n1_traces-1; start_indx_nmo>1; start_indx_nmo--){
      /* pfrintf(stderr,"r_indx[it]=%f, rindx */
      if((r_index_tx_of_it0[start_indx_nmo  ]-
	  r_index_tx_of_it0[start_indx_nmo-1])  <nmostretch) break;
    }
    if(inv){
      start_indx_nmo_tx = 1.0+r_index_tx_of_it0[start_indx_nmo];
      if(start_indx_nmo_tx>n1_traces-2)start_indx_nmo_tx = n1_traces-2;
      /* compute r_index_t0_of_itx from r_index_tx_of_it0 */
      yxtoxy(n1_traces-start_indx_nmo,1.0, start_indx_nmo,
	     &r_index_tx_of_it0[start_indx_nmo],
 
	     n1_traces-start_indx_nmo_tx,1.0, start_indx_nmo_tx,
	     -1,n1_traces,&r_index_t0_of_itx[start_indx_nmo_tx]);
    }
    /* kls inverse nmo? will need more code */
    /* do nmo via 8-point sinc interpolation */
    /* void ints8r (int nxin, float dxin, float fxin,      
                    float yin[], float yinl, float yinr, 
		    int nxout, float xout[], 
		    float yout[]) */
    if(!inv){
      ints8r(n1_traces,1.0,0,
	     intrace,0.0,0.0,
	     n1_traces-start_indx_nmo,&r_index_tx_of_it0[start_indx_nmo],
	     &outtrace[start_indx_nmo]);
      /* zero above the start time */
      for(indx_time=0; indx_time<start_indx_nmo; indx_time++){
	outtrace[indx_time]=0.0;
      }
      /* apply linear ramp kls */
      for (indx_time=start_indx_nmo; 
	   indx_time<start_indx_nmo+lmute && indx_time<n1_traces;
	   indx_time++){
	outtrace[indx_time] *= (float)(indx_time-start_indx_nmo+1)/(float)lmute;
      }
    }else{
      ints8r(n1_traces,1.0,0,
	     intrace,0.0,0.0,
	     n1_traces-start_indx_nmo_tx,&r_index_t0_of_itx[start_indx_nmo_tx],
	     &outtrace[start_indx_nmo_tx]);      
      /* zero above the start time */
      for(indx_time=0; indx_time<start_indx_nmo_tx; indx_time++){
	outtrace[indx_time]=0.0;
      }
    }
    /***************************/
    /* write trace and headers */
    /***************************/
    put_tah(outtrace, fheader, n1_traces, n1_headers, out);
    itrace++;
  }
  
  exit(0);
}

  
