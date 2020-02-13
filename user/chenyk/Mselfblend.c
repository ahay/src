/* Self blending for simple test.*/

/*
  Copyright (C) 2014 Yangkang Chen, University of Texas at Austin
     
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
               
  This program is distributed in the hope that it will be useful,
  but WitHOUT ANY WARRAn1Y; without even the implied warran1y of
  MERCHAn1ABILitY or FitNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
                    
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <string.h>
#include <time.h>
#include <rsf.h>


void blend(int n1, int n2, float dt, float **din, float **dout, float *times)
/* Blending and combing together */
{

  int i1,i2,il,nl;
  float mintime, maxtime;
  float *ltrace;
  /**************************************
   * blend the data into a long trace
   **************************************/
  /* shot_time_in array is in ms */
  mintime=maxtime=times[0];
  for(i2=0; i2<n2; i2++){
    if(mintime>times[i2]) 
      mintime=times[i2];
    if(maxtime<times[i2]) 
      maxtime=times[i2];
  }

  /* dt is in seconds.  shottime array is in s. */
  nl=(maxtime-mintime)/dt+n1;
  ltrace=sf_floatalloc(nl);
  for(il=0; il<nl; il++) ltrace[il]=0.0;
  for(i2=0; i2<n2; i2++)
    {
      int firstindex=(times[i2]-mintime) /dt;
      for(i1=0; i1<n1; i1++){
	ltrace[firstindex+i1]+=din[i2][i1];
      }
    }

  /********************************************
   * extract the data from the long trace *
   * back in the input trace array            *   
   ********************************************/
  for(i2=0; i2<n2; i2++){
    int firstindex=(times[i2]-mintime)/dt;
    memcpy(dout[i2],&(ltrace[firstindex]),n1*sizeof(float));
  }
  free(ltrace);
}


int main(int argc, char* argv[])
{
    int seed;
  int i2,i3,n1,n2,n3; /*l means long trace*/
  float dt,range,var,mean,a,b,t0;
  float *times, **tracesin, **tracesout;
  bool normal;	
  sf_file in=NULL, out=NULL;
	    
	
  /***************************************************/
  /*	Initialization				   */
  /***************************************************/
  sf_init(argc,argv);
  in = sf_input ("in");
  out = sf_output ("out");	

  /***************************************************/
  /*				   */
  /***************************************************/  	
  sf_histint(in,"n1",&n1);
  sf_histint(in,"n2",&n2); 
  n3 = sf_leftsize(in,2); /* left dimensions after the first one */  	  	 	
  sf_histfloat(in,"d1",&dt);

  /***************************************************/
  /*	Allocating memory			   */
  /***************************************************/    
  times=sf_floatalloc(n2);
  tracesin=sf_floatalloc2(n1,n2);
  tracesout=sf_floatalloc2(n1,n2);

  /***************************************************/
  /*	generate random shot scheduling value			   */
  /***************************************************/   
    
  if (!sf_getint("seed",&seed)) seed = time(NULL);
  
  /* random seed */
  init_genrand((unsigned long) seed);

    if (!sf_getbool("type",&normal)) normal=true;
    /* dithering time distribution, y: normal, n: uniform */
    if (!sf_getfloat("var",&var)) {
    
	/* dithering time variance */
	if (!sf_getfloat("range",&range)) {
	    /* dithering time range (default=1) */
	    a = 1.;
	} else {
	    a = normal? 2.*range/9. : 2.*range;
	}
    } else {
	a = normal? sqrtf(var): sqrtf(12*var);
    }

    if (!sf_getfloat("mean",&mean)) mean=0;
    /* noise mean */

    if (!sf_getfloat("t0",&t0)) t0=n1*dt;
    /*scheduled time interval between two consecutive shots in one source*/
        
    b = normal? mean: mean - 0.5*a;
      
   	    if (normal) {
		for (i2=0; i2 < n2; i2++) {
		    times[i2] = a*sf_randn_one_bm() + b + i2*t0;
		    sf_warning("time[%d]=%g,st=%g",times[i2],i2*t0);
		}
	    } else {
		for (i2=0; i2 < n2; i2++) {
		    times[i2] = a*genrand_real1() + b + i2*t0;
		    sf_warning("time[%d]=%g,st=%g",times[i2],i2*t0);
		}
	    }

	       
  for(i3=0;i3<n3;i3++)
    {
      /***************************************************/
      /*	Reading rsf file		  	    */
      /***************************************************/	
      sf_floatread(tracesin[0],n1*n2,in);

	  blend(n1,n2,dt,tracesin,tracesout,times);

      /***************************************************/
      /*	Writing output	 	   */
      /***************************************************/
      sf_floatwrite(tracesout[0],n1*n2,out);
    }
  	
  /***************************************************/
  /*	Free memory	 	   */
  /***************************************************/  	
  	
  free(times);
  free(tracesin[0]);
  free(tracesin);
  free(tracesout[0]);
  free(tracesout);
  
  exit(0);
}





