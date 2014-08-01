/* Blending, or Deblending using seislet domain thresholding.
 */
/*
  Copyright (C) 2014 University of Texas at Austin
     
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
#include <rsf.h>
#include <rsfpwd.h>

void blend(int n1, int n2, float dt, float **din, float **dout, float *times);
void seisthr(int n1, int n2, int thr, char *thrtype, float **din, float **dout);

int main(int argc, char* argv[])
{
  int i1,i2,i3,n1,n2,n3,order,iter,niter,verb; /*l means long trace*/
  float dt,thr,eps,lambda;
  float *times, **tracesin, **tracesout, **tracestmp, **dips;
  char *mode,*type,*thrtype;
  bool unit=false, inv=true, ifinit;	
  sf_file in=NULL, out=NULL, shottime=NULL, dip=NULL, init=NULL;
	    
	
	
  /***************************************************/
  /*	Initialization				   */
  /***************************************************/
  sf_init(argc,argv);
  in = sf_input ("in");
  out = sf_output ("out");	
  shottime=sf_input("shottime");


  /***************************************************/
  /*				   */
  /***************************************************/  	
  sf_histint(in,"n1",&n1);
  sf_histint(in,"n2",&n2); 
  n3 = sf_leftsize(in,2); /* left dimensions after the first one */  	  	 	
  sf_histfloat(in,"d1",&dt);

  if (NULL == (mode=sf_getstring("mode"))) mode="d";
  /* [b-blending,d-deblending] function mode, the default is d  */

  if(mode[0]=='d')
    {	
      if(!sf_getint("niter",&niter)) niter=30;
      /* number of iterations */
    				
      if(!sf_getfloat("thr",&thr)) thr=10;
      /* threshold value (coefficients preserved in percentage) */

      if(!sf_getfloat("lambda",&lambda)) lambda=0.5;
      /* update step size */
    		
      if(!sf_getfloat("eps",&eps)) eps=0.01;
      /* regularization */
    			
      if(!sf_getint("order",&order))order=1;
      /* accuracy order for seislet transform*/ 

      if(!sf_getint("verb",&verb))verb=0;
      /* output verbosity information */
       
      if (NULL == (type=sf_getstring("type"))) type="linear";
      /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

      if (NULL == (thrtype=sf_getstring("thrtype"))) type="soft";
      /* [soft,hard] thresholding type, the default is soft  */   
 
 	 dip=sf_input("dip");
 	  
      if(NULL!=sf_getstring("init")) {ifinit=true; init=sf_input("init"); }  else {ifinit=false;}	
 			
      seislet_init(n1,n2,inv,unit,eps,order,type[0]);  /* unit=false inv=true */    		 		   		
    }
	


  /***************************************************/
  /*	Allocating memory			   */
  /***************************************************/    
  times=sf_floatalloc(n2);
  tracesin=sf_floatalloc2(n1,n2);
  tracesout=sf_floatalloc2(n1,n2);
  if(mode[0]=='d') {dips=sf_floatalloc2(n1,n2); tracestmp=sf_floatalloc2(n1,n2); }


  for(i3=0;i3<n3;i3++)
    {
      /***************************************************/
      /*	Reading rsf file		  	    */
      /***************************************************/	
      sf_floatread(tracesin[0],n1*n2,in);
      sf_floatread(times,n2,shottime);

      if(mode[0]=='b')
	{    
	  blend(n1,n2,dt,tracesin,tracesout,times);
	}else
	{	  
	  sf_floatread(dips[0],n1*n2,dip);
	  seislet_set(dips);

	  if(ifinit) 
	    {
	      sf_floatread(tracesout[0],n1*n2,init);
	    }else{
	    memset(tracesout[0], 0.0, n1*n2*sizeof(float));	
	  }

	  for(iter=0;iter<niter;iter++)
	    {
	      blend(n1,n2,dt,tracesout,tracestmp,times);
		for(i1=0;i1<n1;i1++)
		  for(i2=0;i2<n2;i2++)
		  {  tracestmp[i2][i1]=tracesout[i2][i1]+lambda*(tracesin[i2][i1]-tracestmp[i2][i1]); }
	
	      seisthr(n1,n2,thr,thrtype,tracestmp,tracesout);
	    if(verb==1) sf_warning("Current status: i3=%d, iter=%d",i3,iter);		
	    }

	}
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
  if(mode[0]=='d') {free(dips[0]);free(dips);free(tracestmp[0]);free(tracestmp);} 	
  exit(0);
}

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

  /* dt is in seconds.  shottime array is in ms. */
  nl=(maxtime-mintime)/(dt*1000.)+n1;
  ltrace=sf_floatalloc(nl);
  for(il=0; il<nl; il++) ltrace[il]=0.0;
  for(i2=0; i2<n2; i2++)
    {
      int firstindex=(times[i2]-mintime) /(dt*1000.);
      for(i1=0; i1<n1; i1++){
	ltrace[firstindex+i1]+=din[i2][i1];
      }
    }

  /********************************************
   * extract the data from the long trace *
   * back in the input trace array            *   
   ********************************************/
  for(i2=0; i2<n2; i2++){
    int firstindex=(times[i2]-mintime)/(dt*1000.);
    memcpy(dout[i2],&(ltrace[firstindex]),n1*sizeof(float));
  }
  free(ltrace);
}

void seisthr(int n1, int n2, int thr, char *thrtype, float **din, float **dout)
{
  int i1,i2,nthr;
  float t;
  float **dseis,**dtmp;
  dseis=sf_floatalloc2(n1,n2);
  dtmp=sf_floatalloc2(n1,n2);
		
  seislet_lop(true,false,n1*n2,n1*n2,dseis[0],din[0]); 

  /***********************************************************/
  /* percentile thresholding */
  for(i2=0;i2<n2;i2++)
    for(i1=0;i1<n1;i1++)
      {dtmp[i2][i1]=fabs(dseis[i2][i1]); }	
	
  nthr = 0.5+n1*n2*(1.-0.01*thr);  /*round off*/
  if (nthr < 0) nthr=0;
  if (nthr >= n1*n2) nthr=n1*n2-1;
  t=sf_quantile(nthr,n1*n2,dtmp[0]);
	
  for(i2=0;i2<n2;i2++)
    for(i1=0;i1<n1;i1++)
      if(fabs(dseis[i2][i1])<t) {dseis[i2][i1]=0;   /* hard thresholding */}
      else
	{if(thrtype[0]=='s' && dseis[i2][i1]>0) dseis[i2][i1]=dseis[i2][i1]-t;   /* soft thresholding */
	  else if(thrtype[0]=='s' && dseis[i2][i1]<0) dseis[i2][i1]=dseis[i2][i1]+t;}
  /***********************************************************/
	
  seislet_lop(false,false,n1*n2,n1*n2,dseis[0],dout[0]);

	
  free(dseis[0]);
  free(dseis);
  free(dtmp[0]);
  free(dtmp);	
}
