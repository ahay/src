/* FX Regularized Nonstationary Autoregression subroutine 

Reference: Wang et al. (2021), Non-stationary predictive filtering for seismic random noise suppression - A tutorial, Geophysics, 86(3), W21–W30. 
*/
/*
  Copyright (C) 2016 Yangkang Chen
*/

#include <rsf.h>
#include "fft1.h"
#include "cdivn.h"
#include "cdivn_fs.h"
#include "cdivn_rnar.h"
#include "cdivnn.h"
#include "cdivns.h" /*stationary division in traditional FXY*/
void fxynpre0(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		bool sym, 
		bool opt,
		bool verb)	
/*<Non-stationary predictive filtering in 3D (without frequency-dependent smoothing)>*/ 
{   
    int nw, nt, i2, i22;
    int m[SF_MAX_DIM], dim;    
    int nd, n12;    
    int ns,is,i1,isx,isy,ix,iy;
    float dw, mean;
    
	/* Part I:  Fourier transform */	
 	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	/*	if (nt0!=nt) sf_error("nt should be consistant");*/

	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);
    
    kiss_fft_cpx **dfft;
    
    dfft = (kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		             
    /*forward transform*/
    fft1(dtime, dfft, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, false);
	
	/* Part II:  shifting the data */
	ns=(nsx*2+1)*(nsy*2+1)-1;

	kiss_fft_cpx ***dfftfilt;		
	kiss_fft_cpx **dfftpre;
	kiss_fft_cpx ***dfftshift;

	dfftfilt=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);		
	dfftshift=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);	
	dfftpre=(kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		

	is=-1;
	for(isx=-nsx;isx<nsx+1;isx++)
		for(isy=-nsy;isy<nsy+1;isy++)
		{
			if(!(isx==0 && isy==0))
			{
			is++;
			for(iy=0;iy<ny;iy++)
				for(ix=0;ix<nx;ix++)
					for(i1=0;i1<nw;i1++)
					{
					i2=iy*nx+ix;
					if((iy+isy)>=0 && (iy+isy)<ny && (ix+isx)>=0 && (ix+isx)<nx)
					{
						i22=(iy+isy)*nx+(ix+isx);
						dfftshift[is][i2][i1]= dfft[i22][i1];
					}
					}
			}
		}	
	/*Part III: Local prediction filter for complex numbers*/
	
    nd = nw*n2;
    n12=nd*ns;
    
    m[0]=nw;m[1]=nx;m[2]=ny;
	if(ny==1)	/*2D denoising*/
		dim=2;
	else		/*3D denoising*/
		dim=3;
	
	sf_warning("nsx=%d,nsy=%d,ns=%d",nsx,nsy,ns);
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);

    cmultidivn_init(ns, dim, nd, m, rect, (sf_complex *) dfftshift[0][0], verb); 
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);

    mean = 0.;
        
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	mean += sf_crealf(sf_cmul(sf_conjf(dfftshift[is][i2][i1]),dfftshift[is][i2][i1]));
    }

    mean = sqrtf (n12/mean);
    
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],mean);
    }
    
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfft[i2][i1] = sf_crmul(dfft[i2][i1],mean);
    }
    
    cmultidivn ((sf_complex *) dfft[0], (sf_complex *) dfftfilt[0][0],niter);
	

    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
	{
	    dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],1./mean);
	}
	
	cweight2_lop(false,false,n12,nd,(sf_complex *) dfftfilt[0][0], (sf_complex *) dfftpre[0]);
		
    /*inverse transform*/
    fft1(dtime, dfftpre, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, true);   
    
	cmultidivn_close();  
    free(**dfftfilt);free(*dfftfilt);free(dfftfilt);
    free(**dfftshift);free(*dfftshift);free(dfftshift);
    free(*dfftpre);free(dfftpre);
    free(*dfft);free(dfft);
    
}

void fxynpre(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		float Nfrac, /*starting varying from 1/Nfrac*F_nyquist*/
		float Ntimes, /*largest radius is Ntimes of the ref radius*/
		float pow,   /*sharp rate of the varying curve, the lower the sharper*/
		bool sym, 
		bool opt,
		bool verb)	
/*<Non-stationary predictive filtering in 3D with frequency-dependent smoothing>*/ 
{   
    int nw, nt, i2, i22;
    int m[SF_MAX_DIM], dim;    
    int nd, n12;    
    int ns,is,i1,isx,isy,ix,iy;
    float dw, mean;
    
	/* Part I:  Fourier transform */	
 	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	/*	if (nt0!=nt) sf_error("nt should be consistant");*/

	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);
    
    kiss_fft_cpx **dfft;
    
    dfft = (kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		             
    /*forward transform*/
    fft1(dtime, dfft, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, false);
	
	/* Part II:  shifting the data */
	ns=(nsx*2+1)*(nsy*2+1)-1;

	kiss_fft_cpx ***dfftfilt;		
	kiss_fft_cpx **dfftpre;
	kiss_fft_cpx ***dfftshift;

	dfftfilt=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);		
	dfftshift=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);	
	dfftpre=(kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		

	is=-1;
	for(isx=-nsx;isx<nsx+1;isx++)
		for(isy=-nsy;isy<nsy+1;isy++)
		{
			if(!(isx==0 && isy==0))
			{
			is++;
			for(iy=0;iy<ny;iy++)
				for(ix=0;ix<nx;ix++)
					for(i1=0;i1<nw;i1++)
					{
					i2=iy*nx+ix;
					if((iy+isy)>=0 && (iy+isy)<ny && (ix+isx)>=0 && (ix+isx)<nx)
					{
						i22=(iy+isy)*nx+(ix+isx);
						dfftshift[is][i2][i1]= dfft[i22][i1];
					}
					}
			}
		}	
	/*Part III: Local prediction filter for complex numbers*/
	
    nd = nw*n2;
    n12=nd*ns;
    
    m[0]=nw;m[1]=nx;m[2]=ny;
	if(ny==1)	/*2D denoising*/
		dim=2;
	else		/*3D denoising*/
		dim=3;

	int **rct;float quad;
	rct=sf_intalloc2(2,nw);
	for(i1=0;i1<nw;i1++)
	{
		if(i1<=nw/6)
		{
		rct[0][i1]=rect[1];
		rct[1][i1]=rect[2];
		}else{
		quad = ((i1-nw/Nfrac)*(i1-nw/Nfrac)/(nw-1-nw/Nfrac)/(nw-1-nw/Nfrac));/*default: Nfrac=6.0*/
		quad = powf(quad,0.5*pow);/*default: pow=0.5*/
		rct[0][i1]=rect[1]+quad*(Ntimes*rect[1]-rect[1]); /*from reference radius to Ntimes of ref gradually*//*default: Ntimes=5*/
		rct[1][i1]=rect[2]+quad*(Ntimes*rect[2]-rect[2]); /*from reference radius to Ntimes of ref gradually*/
		}
/*		sf_warning("R[%d]=%d",i1,rct[0][i1]);*/
	}
	sf_warning("Rf = %d is unchanged",rect[0]);
	sf_warning("Rx changes from %d to %d",rect[1],rct[0][nw-1]);
	sf_warning("Ry changes from %d to %d",rect[2],rct[1][nw-1]);

    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);

    cmultidivn_fs_init(ns, dim, nd, m, rect[0], rct, (sf_complex *) dfftshift[0][0], verb); 
    /*frequency-dependent smoothing*/
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);

    mean = 0.;
        
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	mean += sf_crealf(sf_cmul(sf_conjf(dfftshift[is][i2][i1]),dfftshift[is][i2][i1]));
    }

    mean = sqrtf (n12/mean);
    
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],mean);
    }
    
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfft[i2][i1] = sf_crmul(dfft[i2][i1],mean);
    }
    cmultidivn_fs ((sf_complex *) dfft[0], (sf_complex *) dfftfilt[0][0],niter);
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
	{
	    dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],1./mean);
	}
	cweight2_lop(false,false,n12,nd,(sf_complex *) dfftfilt[0][0], (sf_complex *) dfftpre[0]);
	
    /*inverse transform*/
    fft1(dtime, dfftpre, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, true);   
    
	cmultidivn_fs_close();    
    free(**dfftfilt);free(*dfftfilt);free(dfftfilt);
    free(**dfftshift);free(*dfftshift);free(dfftshift);
    free(*dfftpre);free(dfftpre);
    free(*dfft);free(dfft);
    
    
}

void fxynpre1(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		bool sym, 
		bool opt,
		bool verb)	
/*<Non-stationary predictive filtering in 3D>*/ 
{   
    int nw, nt, i2, i22;
    int m[SF_MAX_DIM], dim;    
    int nd, n12;    
    int ns,is,i1,isx,isy,ix,iy;
    float dw, mean;
    
	/* Part I:  Fourier transform */	
 	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	/*	if (nt0!=nt) sf_error("nt should be consistant");*/

	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);
    
    kiss_fft_cpx **dfft;
    
    dfft = (kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		             
    /*forward transform*/
    fft1(dtime, dfft, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, false);
	
	/* Part II:  shifting the data */
	ns=(nsx*2+1)*(nsy*2+1)-1;

	kiss_fft_cpx ***dfftfilt;		
	kiss_fft_cpx **dfftpre;
	kiss_fft_cpx ***dfftshift;

	dfftfilt=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);		
	dfftshift=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);	
	dfftpre=(kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		

	is=-1;
	for(isx=-nsx;isx<nsx+1;isx++)
		for(isy=-nsy;isy<nsy+1;isy++)
		{
			if(!(isx==0 && isy==0))
			{
			is++;
			for(iy=0;iy<ny;iy++)
				for(ix=0;ix<nx;ix++)
					for(i1=0;i1<nw;i1++)
					{
					i2=iy*nx+ix;
					if((iy+isy)>=0 && (iy+isy)<ny && (ix+isx)>=0 && (ix+isx)<nx)
					{
						i22=(iy+isy)*nx+(ix+isx);
						dfftshift[is][i2][i1]= dfft[i22][i1];
					}
					}
			}
		}	
	/*Part III: Local prediction filter for complex numbers*/
	
    nd = nw*n2;
    n12=nd*ns;
    
    m[0]=nw;m[1]=nx;m[2]=ny;
	if(ny==1)	/*2D denoising*/
		dim=2;
	else		/*3D denoising*/
		dim=3;
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);
    cmultidivn_rnar_init(ns, dim, nd, m, rect, (sf_complex *) dfftshift[0][0], verb); 
    if(verb) sf_warning("rect[0]=%d,rect[1]=%d,rect[2]=%d",rect[0],rect[1],rect[2]);
    

    mean = 0.;
        
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	mean += sf_crealf(sf_cmul(sf_conjf(dfftshift[is][i2][i1]),dfftshift[is][i2][i1]));
    }

    mean = sqrtf (n12/mean);
    
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],mean);
    }
    
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfft[i2][i1] = sf_crmul(dfft[i2][i1],mean);
    }

    cmultidivn_rnar ((sf_complex *) dfft[0], (sf_complex *) dfftfilt[0][0],niter);

    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
	{
	    dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],1./mean);
	}
	
	cweight2_lop(false,false,n12,nd,(sf_complex *) dfftfilt[0][0], (sf_complex *) dfftpre[0]);
		
    /*inverse transform*/
    fft1(dtime, dfftpre, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, true);   
    
	cmultidivn_rnar_close();    
    free(**dfftfilt);free(*dfftfilt);free(dfftfilt);
    free(**dfftshift);free(*dfftshift);free(dfftshift);
    free(*dfftpre);free(dfftpre);
    free(*dfft);free(dfft);

}

void fxynpre3(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		float **rct /* triangle lengths [ndim][nd] */,
        int **sft 	/* triangle shifts [ndim][nd] */,
		bool sym, 
		bool opt,
		bool verb)	
/*<Non-stationary predictive filtering in 3D with arbitrary smoothing strategies>*/ 
{   
    int nw, nt, i2, i22;
    int m[SF_MAX_DIM], dim;    
    int nd, n12;    
    int ns,is,i1,isx,isy,ix,iy;
    float dw, mean;
    
	/* Part I:  Fourier transform */	
 	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	/*	if (nt0!=nt) sf_error("nt should be consistant");*/

	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);
    
    kiss_fft_cpx **dfft;
    
    dfft = (kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		             
    /*forward transform*/
    fft1(dtime, dfft, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, false);
	
	/* Part II:  shifting the data */
	ns=(nsx*2+1)*(nsy*2+1)-1;

	kiss_fft_cpx ***dfftfilt;		
	kiss_fft_cpx **dfftpre;
	kiss_fft_cpx ***dfftshift;

	dfftfilt=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);		
	dfftshift=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);	
	dfftpre=(kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		

	is=-1;
	for(isx=-nsx;isx<nsx+1;isx++)
		for(isy=-nsy;isy<nsy+1;isy++)
		{
			if(!(isx==0 && isy==0))
			{
			is++;
			for(iy=0;iy<ny;iy++)
				for(ix=0;ix<nx;ix++)
					for(i1=0;i1<nw;i1++)
					{
					i2=iy*nx+ix;
					if((iy+isy)>=0 && (iy+isy)<ny && (ix+isx)>=0 && (ix+isx)<nx)
					{
						i22=(iy+isy)*nx+(ix+isx);
						dfftshift[is][i2][i1]= dfft[i22][i1];
					}
					}
			}
		}	
	/*Part III: Local prediction filter for complex numbers*/
	
    nd = nw*n2;
    n12=nd*ns;
    
    m[0]=nw;m[1]=nx;m[2]=ny;
	if(ny==1)	/*2D denoising*/
		dim=2;
	else		/*3D denoising*/
		dim=3;

// 	int **rct;
// 	rct=sf_intalloc2(2,nw);
// 	for(i1=0;i1<nw;i1++)
// 	{
// 		if(i1<=nw/2)
// 		{
// 		rct[1][i1]=7;
// 		rct[2][i1]=7;
// 		}else{
// 		rct[1][i1]=20;
// 		rct[2][i1]=20;	
// 		}
// 	}
	
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);

    cmultidivnn_init(ns, dim, nd, m, rect, rct, sft, (sf_complex *) dfftshift[0][0], verb); 
    /*frequency-dependent smoothing*/
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);

    mean = 0.;
        
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	mean += sf_crealf(sf_cmul(sf_conjf(dfftshift[is][i2][i1]),dfftshift[is][i2][i1]));
    }

    mean = sqrtf (n12/mean);
    
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],mean);
    }
    
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfft[i2][i1] = sf_crmul(dfft[i2][i1],mean);
    }
    sf_warning("cmultidivnn starts!");
    cmultidivnn ((sf_complex *) dfft[0], (sf_complex *) dfftfilt[0][0],niter);
	sf_warning("cmultidivnn finishes!");
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
	{
	    dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],1./mean);
	}
	
	sf_warning("cweight starts!");
	cweight2_lop(false,false,n12,nd,(sf_complex *) dfftfilt[0][0], (sf_complex *) dfftpre[0]);
	sf_warning("cweight finishes!");
	
    /*inverse transform*/
    fft1(dtime, dfftpre, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, true);   
    
	cmultidivnn_close();    
    free(**dfftfilt);free(*dfftfilt);free(dfftfilt);
    free(**dfftshift);free(*dfftshift);free(dfftshift);
    free(*dfftpre);free(dfftpre);
    free(*dfft);free(dfft);
    
}

void fxypre(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		bool sym, 
		bool opt,
		bool verb)	
/*<Stationary predictive filtering in 3D>*/ 
{   
    int nw, nt, i2, i22;
    int m[SF_MAX_DIM], dim;    
    int nd, n12;    
    int ns,is,i1,isx,isy,ix,iy;
    float dw, mean;
    
	/* Part I:  Fourier transform */	
 	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	/*	if (nt0!=nt) sf_error("nt should be consistant");*/

	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);
    
    kiss_fft_cpx **dfft;
    
    dfft = (kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		             
    /*forward transform*/
    fft1(dtime, dfft, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, false);
	
	/* Part II:  shifting the data */
	ns=(nsx*2+1)*(nsy*2+1)-1;

	kiss_fft_cpx **dfftfilt;		
	kiss_fft_cpx **dfftpre;
	kiss_fft_cpx ***dfftshift;

	dfftfilt=(kiss_fft_cpx**)  sf_complexalloc2(nw,ns);		
	dfftshift=(kiss_fft_cpx***)  sf_complexalloc3(nw,n2,ns);	
	dfftpre=(kiss_fft_cpx**)  sf_complexalloc2(nw,n2);
		

	is=-1;
	for(isx=-nsx;isx<nsx+1;isx++)
		for(isy=-nsy;isy<nsy+1;isy++)
		{
			if(!(isx==0 && isy==0))
			{
			is++;
			for(iy=0;iy<ny;iy++)
				for(ix=0;ix<nx;ix++)
					for(i1=0;i1<nw;i1++)
					{
					i2=iy*nx+ix;
					if((iy+isy)>=0 && (iy+isy)<ny && (ix+isx)>=0 && (ix+isx)<nx)
					{
						i22=(iy+isy)*nx+(ix+isx);
						dfftshift[is][i2][i1]= dfft[i22][i1];
					}
					}
			}
		}	
	/*Part III: Local prediction filter for complex numbers*/
	
    nd = nw*n2;
    n12=nd*ns;
    
    m[0]=nw;m[1]=1;m[2]=1;

	if(ny==1)	/*2D denoising*/
		dim=1;
	else		/*3D denoising*/
		dim=2;
		
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);
    
    rect[1]=1;rect[2]=1;
    if(verb) sf_warning("rect[0]=%d,rect[1]=%d,rect[2]=%d",rect[0],rect[1],rect[2]);
    cmultidivns_init(ns, 1, nd, m, rect, (sf_complex *) dfftshift[0][0], verb); 
    if(verb) sf_warning("n12=%d,nd=%d,ns=%d",n12,nd,ns);

    mean = 0.;
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	mean += sf_crealf(sf_cmul(sf_conjf(dfftshift[is][i2][i1]),dfftshift[is][i2][i1]));
    }
    mean = sqrtf (n12/mean);
    
    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],mean);
    }
    
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
    {
	dfft[i2][i1] = sf_crmul(dfft[i2][i1],mean);
    }
    cmultidivns ((sf_complex *) dfft[0], (sf_complex *) dfftfilt[0],niter);
	

    for(is=0;is<ns;is++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<nw;i1++)
	{
	    dfftshift[is][i2][i1] = sf_crmul(dfftshift[is][i2][i1],1./mean);
	}
	cweight2s_lop(false,false,nw*ns,nd,(sf_complex *) dfftfilt[0], (sf_complex *) dfftpre[0]);
	
    /*inverse transform*/
    fft1(dtime, dfftpre, n2,n1,d1,o1,nt,nw,dw, sym, opt, verb, true);   
    
	cmultidivns_close();    
    free(*dfftfilt);free(dfftfilt);
    free(**dfftshift);free(*dfftshift);free(dfftshift);
    free(*dfftpre);free(dfftpre);
    free(*dfft);free(dfft);
    
}


