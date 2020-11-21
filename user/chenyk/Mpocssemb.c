/* POCS using semblance thresholding or amplitude thresholding.
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

#include <rsf.h>

#include <stdio.h>
#include <math.h>

#include "fft.h"
#include "taup.h"

int main(int argc, char* argv[])
{


  int nt, nx, it, ix, niter, iter, ntfft, nxfft,np, ip, ikt, ikx, iktn, ikxn, ifsnr; /* iktn, ikxn, iNyquist*/
  float dt, dx, pmin, pmax = 0.0, dp=0, p, cmax, sembpmax, num, den;
  float *sembp=NULL, *mask, *gy, *fden, *fshift, *SNR=NULL;
  float **fdata, **taup=NULL, **odata, **tdata=NULL, **odatat, **semb=NULL; /* tdata is the true data */
  kiss_fft_cpx **cdata, **cdatat;
  char *type;
  sf_file inp, outp, m, spec1=NULL, spec2=NULL, trued, snr=NULL; 


  sf_init(argc,argv);

  inp=sf_input("in");
  m=sf_input("mask");
  outp=sf_output("out");

  if(!sf_histint(inp,"n1",&nt)) sf_warning("No n1 in input");
  if(!sf_histint(inp,"n2",&nx)) sf_warning("No n2 in input");
  if(!sf_histfloat(inp,"d1",&dt)) sf_warning("No n1 in input");
  if(!sf_histfloat(inp,"d2",&dx)) sf_warning("No n2 in input");

  ntfft = 2*kiss_fft_next_fast_size((nt+1)/2);
  nxfft = 2*kiss_fft_next_fast_size((nx+1)/2);
  iktn=ntfft/2; ikxn=nxfft/2;
  float dkt = 1.0/(ntfft*dt), fkt = 0.0,kt;
  float dkx = 1.0/(nxfft*dx), fkx = 0.0,kx;

    if (NULL == (type=sf_getstring("type"))) type="amplitude";
    /* [amplitude, semblance] thresholding type, the default is amplitude thresholding  */ 

  	if(!sf_getint("niter",&niter)) niter = 10;
	/* Get the number of iterations */

  	if(!sf_getint("ifsnr",&ifsnr)) ifsnr = 0;
	/* If compute SNR during iteration */


	if(type[0]=='s')
	{

  		if(!sf_getfloat("pmin",&pmin)) pmin=-2;
        /* minimum p */		
  		if(!sf_getfloat("pmax",&pmax)) pmax=2;
        /* maximum p */			
  		if(!sf_getint("np",&np)) np=nx;
        /* number of p */
		dp=(pmax-pmin)/(np-1);		
  		sembp =sf_floatalloc(np);
  		semb  =sf_floatalloc2(nt,np);
 		taup  =sf_floatalloc2(nt,np);	
	}

	/* output files */
  	if (NULL!=sf_getstring("spec2")) 
  	{
		spec2=sf_output("spec2");
		sf_putint(spec2, "n1", ntfft);
		sf_putint(spec2, "n2", nxfft);
	}	       		
  	if (NULL!=sf_getstring("spec1")) 
  	{
		spec1=sf_output("spec1");
		sf_putint(spec1, "n1", ntfft);
		sf_putint(spec1, "n2", nxfft);
	}	  
  	if (ifsnr==1 && (NULL!=sf_getstring("true"))) 
  	{
		snr=sf_output("snr");
		trued=sf_input("true");
		tdata=sf_floatalloc2(nt,nx);
		SNR=sf_floatalloc(niter);
		sf_floatread(tdata[0],nt*nx,trued);

		sf_putint(snr,"n1",niter);
		sf_putint(snr,"d1",1);
		sf_putint(snr,"n2",1);
	}	

  /* Allocate memory */
  cdata =(kiss_fft_cpx**) sf_complexalloc2(ntfft,nxfft);
  cdatat =(kiss_fft_cpx**) sf_complexalloc2(ntfft,nxfft); /* temporary file */
  fshift= sf_floatalloc(ntfft);
  fden 	= sf_floatalloc(ntfft);
  gy 	= sf_floatalloc(nxfft);
  mask  =sf_floatalloc(nx);

  odata =sf_floatalloc2(nt,nx); 
  odatat =sf_floatalloc2(ntfft,nxfft); 
  fdata =sf_floatalloc2(ntfft,nxfft);
  memset(&odata[0][0],0,ntfft*nxfft*sizeof(float)); 

  /* Read data */
  sf_floatread(odata[0],nt*nx,inp);
  sf_floatread(mask,nx,m);

	if(type[0]=='s')
	{
   	 slant(dt,nt,nx,gy,pmin,dp,np,odata,taup,semb,sembp,&sembpmax,fshift,fden);
	}

  for (iter=niter-1; iter>=0; iter--) {
    tfft(odata, cdatat, nx, ntfft);
    xfft(cdatat, cdata, ntfft, nxfft);
    cmax = findmax(nxfft,ntfft,cdata);

    if (iter==0 || iter==niter-1) { // beginning and ending spectra
  		for (ix=0; ix<nxfft; ix++)
  		for (it=0; it<ntfft; it++)
    	fdata[ix][it] = sf_cabsf(cdata[ix][it]);

  		if (iter==0 && (NULL!=sf_getstring("spec2")))        		sf_floatwrite(fdata[0],ntfft*nxfft,spec2);
  		if (iter==niter-1 && (NULL!=sf_getstring("spec1")))   	sf_floatwrite(fdata[0],ntfft*nxfft,spec1);

    }



	if(type[0]=='a')
	{
    	for (ix=0; ix<nxfft; ix++) // Abma Kabir FT amplitude thresholding
    		for (it=0; it<ntfft; it++)
      			if (sf_cabsf(cdata[ix][it])<iter*1./niter*cmax) cdata[ix][it] = cmplx(0.,0.);
	}
	else
	{
      for (ix=0; ix<nxfft; ix++) // Abma Kabir FT amplitude thresholding
    		for (it=0; it<ntfft; it++)
      		if (sf_cabsf(cdata[ix][it])<iter*1./niter*cmax) cdata[ix][it] = cmplx(0.,0.);


      for (ikx=0,kx=fkx; ikx<=ikxn; ++ikx,kx+=dkx) {
    		for (ikt=0,kt=fkt; ikt<=iktn; ++ikt,kt+=dkt) {
      		if (kx==0) {
        		if (sf_cabsf(cdata[ikx][ikt])<iter*1./niter*cmax) cdata[ikx][ikt] = cmplx(0.,0.);
        		continue;
      		}

      p = -kx/kt; ip = round((p-pmin)/dp);
      //if (ip<0 || ip>=np) { cdata[ikx][ikt] = 0.;continue; }
      if (ip<0 || ip>=np) {  }
      else if (sembp[ip] <iter*1./niter*sembpmax) cdata[ikx][ikt] = cmplx(0.,0.);

      if (ikx>0 && ikx<(nxfft+1)/2) { // kx>=0, kx<0
      p = kx/kt; ip = round((p-pmin)/dp);
      //if (ip<0 || ip>=np) { cdata[nxfft-ikx][ikt] = 0.;continue; }
      if (ip<0 || ip>=np) {  }
      else if (sembp[ip] <iter*1./niter*sembpmax) cdata[nxfft-ikx][ikt] = cmplx(0.,0.);
      }

      if (ikt>0 && ikt<(ntfft+1)/2) { // kt<0, kx>=0
      p = kx/kt; ip = round((p-pmin)/dp);
      //if (ip<0 || ip>=np) { cdata[ikx][ntfft-ikt] = 0.;continue; }
      if (ip<0 || ip>=np) {  }
      else if (sembp[ip] <iter*1./niter*sembpmax) cdata[ikx][ntfft-ikt] = cmplx(0.,0.);
      }

      if (ikx>0 && ikx<(nxfft+1)/2 && ikt>0 && ikt<(ntfft+1)/2) { // kt<0, kx<0
      p = -kx/kt; ip = round((p-pmin)/dp);
      //if (ip<0 || ip>=np) { cdata[nxfft-ikx][ntfft-ikt] = 0.;continue; }
      if (ip<0 || ip>=np) {  }
      else if (sembp[ip] <iter*1./niter*sembpmax) cdata[nxfft-ikx][ntfft-ikt] =cmplx(0.,0.);
      }
    }}
	
	}

	ixfft(cdata, cdatat, ntfft, nxfft);
	itfft(cdatat, odatat, nxfft, ntfft);

    for (ix=0; ix<nx; ix++) { // put in ORIGINAL KNOWN data
      if (mask[ix]==1) continue;
      for (it=0; it<nt; it++) odata[ix][it]=odatat[ix][it];
    }

	num=0;den=0;
	/* If output the SNR file. */
  	if (ifsnr==1 && (NULL!=sf_getstring("true"))) 
  	{
		for(ix=0;ix<nx;ix++)
			for(it=0;it<nt;it++)
			{
				num+=tdata[ix][it]*tdata[ix][it];
				den+=(tdata[ix][it]-odata[ix][it])*(tdata[ix][it]-odata[ix][it]);
			}	
		SNR[niter-iter-1]=10*logf(num/den);
	}	
  }

  	if (ifsnr==1 && (NULL!=sf_getstring("true"))) 
	{ sf_floatwrite(SNR,niter,snr);	}

  	sf_floatwrite(odata[0],nt*nx,outp);
    exit (0);
}



