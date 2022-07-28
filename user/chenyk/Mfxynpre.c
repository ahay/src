/* FXY non-stationary predictive filtering (with the option of local processing windows). 
with frequency-axis smoothness constraint
with frequency-dependent smoothing

Reference: Wang et al. (2021), Non-stationary predictive filtering for seismic random noise suppression - A tutorial, Geophysics, 86(3), W21-W30. 
*/

/*
  Copyright (C) 2020 Yangkang Chen
*/

#include <rsf.h>
#include "fxynpre.h"
#include "memcpy.h"
#include "win.h"

int main (int argc, char *argv[])
{
  
    bool  verb, sym, opt, fs; /*fs: if apply frequency-dependent smoothing*/
    int   n1, ix,iy,nx,ny;
    int rect[SF_MAX_DIM];    
    int niter,mode;    
    int nsx,nsy,j;
    int n1win,n2win,n3win,n1pad,n2pad,n3pad;/*window sizes,and padding sizes*/
    int nw1,nw2,nw3,iw1,iw2,iw3;/*number of windows in first and second axis,indices of nw1,nw2*/
    int s1,s2,s3; /*indices in selecting windows*/
    int nov1,nov2,nov3,ov1,ov2,ov3;/*overlapping or non-overlapping points*/
    char key[6];
    float o1,o11, d1,r1,r2,r3;
    float ***din, ***dout,***dtmp;    
    sf_file Fin=NULL, Fou=NULL;

    /*------------------------------------------------------------*/
    sf_init(argc, argv);

    Fin  = sf_input ("in");
    Fou = sf_output("out");

    if(! sf_getbool("verb",&verb))  verb=false; /* Verbosity flag */
    if (!sf_getbool("sym",&sym))    sym=false; /* y, symmetric scaling for Hermitian FFT */
    if (!sf_getbool("opt",&opt))    opt=true;  /* y, determine optimal size for efficiency */
    if (!sf_getbool("fs",&fs))    fs=false;  /* y, determine frequency-dependent smoothing */

    float Nfrac; /*starting varying from nw/Nfrac*F_nyquist*/
    float Ntimes; /*largest radius is Ntimes of the ref radius*/
    float pow;   /*sharp rate of the varying curve, the lower the sharper*/
    
	if(fs)
	{
    if (!sf_getfloat("Nfrac",&Nfrac))    Nfrac=6.0;     /*frequency-dependent smoothing starts from 1/Nfrac * (Nyquist frequency) */
    if (!sf_getfloat("Ntimes",&Ntimes))    Ntimes=5.0;  /*Maximum smoothing radius is Ntimes*(reference smoothing radius) */
    if (!sf_getfloat("pow",&pow))    pow=0.5;  /*fraction parameter*/
	}

    for (j=0; j < 3; j++) {
	snprintf(key,6,"rect%d",j+1);
	if (!sf_getint(key,rect+j)) rect[j]=1; 
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
    }
    
 	if (!sf_histint  (Fin,"n1",&n1)) n1=1;
	if (!sf_histfloat(Fin,"d1",&d1)) d1=1.;
	if (!sf_histfloat(Fin,"o1",&o1)) o1=0.;
 	if (!sf_histint  (Fin,"n2",&nx)) nx=1;
 	if (!sf_histint  (Fin,"n3",&ny)) ny=1;


    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

	/*windowing parameters*/
    if (!sf_getint("n1win",&n1win)) n1win=n1; /*first window length*/
    if (!sf_getint("n2win",&n2win)) n2win=nx; /*second window length*/	
    if (!sf_getint("n3win",&n3win)) n3win=ny; /*second window length*/	
    
	if (!sf_getfloat("r1",&r1)) r1=0.5;		  /*first overlapping ratio*/
	if (!sf_getfloat("r2",&r2)) r2=0.5;		  /*second overlapping ratio*/
	if (!sf_getfloat("r3",&r3)) r3=0.5;		  /*third overlapping ratio*/

    if (!sf_getint("mode",&mode)) mode=0; /*predictive filtering mode; default: non-stationary*/	
    
	if(! sf_getint("nsx",&nsx)) nsx=2; /* number of shifts in non-causal prediction filtering */
	if(! sf_getint("nsy",&nsy)) nsy=2; /* number of shifts in non-causal prediction filtering */

	nov1=(1-r1)*n1win;nov2=(1-r2)*n2win;nov3=(1-r3)*n3win;/*non-overlapping size 1,2,3*/
	ov1=r1*n1win;ov2=r2*n2win;ov3=r3*n3win;		/*overlapping size 1,2,3*/
	
	sf_warning("nsx=%d,nsy=%d",nsx,nsy);

	/*padding */	
	n1pad=n1win;nw1=1;while(n1pad<n1){n1pad=n1pad+nov1;nw1++;	 }   
	n2pad=n2win;nw2=1;while(n2pad<nx){n2pad=n2pad+nov2;nw2++;	 }   	
	n3pad=n3win;nw3=1;while(n3pad<ny){n3pad=n3pad+nov3;nw3++;	 }   	
			    
    din =                  sf_floatalloc3  (n1pad,n2pad,n3pad); 
    dout =                  sf_floatalloc3  (n1pad,n2pad,n3pad);           
    dtmp =                   sf_floatalloc3  (n1win,n2win,n3win);    

 	for(iy=0; iy<ny; iy++)			    
 	for(ix=0; ix<nx; ix++)
		sf_floatread (din[iy][ix],n1,Fin);
		
	for(iw3=0;iw3<nw3;iw3++)
	for(iw2=0;iw2<nw2;iw2++)
		for(iw1=0;iw1<nw1;iw1++)
			{
			s1=iw1*nov1;s2=iw2*nov2;s3=iw3*nov3;
/*				sf_warning("Hello,n1win=%d,n2win=%d",n1win,n2win);
			sf_warning("nw2=%d,nw1=%d,d1=%g,o11=%g,niter=%d,ns=%d",nw2,nw1,d1,o11,niter,ns);*/
			mcp3d(dtmp,din,0,0,0,s1,s2,s3,n1win,n2win,n3win);
			o11=o1+s1*d1;
			if(mode==0)
    		{sf_warning("Non-stationary Predictive Filtering");
    		if(fs)
    		{fxynpre(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,Nfrac,Ntimes,pow,sym,opt,verb);}
    		else
    		{fxynpre0(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,sym,opt,verb);}}
    		else{
    		if(mode==1)
    		{sf_warning("Regularized Non-stationary Autoregression");
    		fxynpre1(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,sym,opt,verb);}
    		else
    		{if(mode==2){sf_warning("Stationary Predictive Filtering");
    		fxypre(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,sym,opt,verb);}else{sf_error("wrong mode");}}
    		}    		
    		win_weight3d(dtmp,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3);
    		mcp_ad3d(dout,dtmp,s1,s2,s3,0,0,0,n1win,n2win,n3win);
						
			}

 	for(iy=0; iy<ny; iy++)			    
 	for(ix=0; ix<nx; ix++)
		sf_floatwrite (dout[iy][ix],n1,Fou); 	

	/*free memory*/
	free(**din); free(*din); free(din);
	free(**dout); free(*dout); free(dout);
	free(**dtmp); free(*dtmp); free(dtmp);
	
    exit (0);
}

