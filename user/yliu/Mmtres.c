/* Calculate apparent resistivity and phase of MT data. */
/*
  Copyright (C) 2013 Jilin University
  
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

int main (int argc, char* argv[])
{
    bool opt, comp, phase, verb;
    double n1;
    float b, a, pi;
    int nt, nw, i, j, N, m1, M, k, f;                    
    float *TDEx,*TDEy,*TDHx,*TDHy,*outp,*ExHyres,*EyHxres,*ExHypha,*EyHxpha;
    float *bb, *pp, *qq, *dd;
    sf_file in, out, Ey, Hx, Hy;

    kiss_fft_cpx *FDEx=NULL,*FDEy=NULL,*FDHx=NULL,*FDHy=NULL;
    kiss_fft_cpx *A=NULL,*B=NULL,*ExAs1=NULL,*HyBs1=NULL,*ExBs1=NULL;
    kiss_fft_cpx *HyAs1=NULL,*HxAs1=NULL,*HxBs1=NULL;
    kiss_fft_cpx *EyAs1=NULL,*EyBs1=NULL,*ExAs2=NULL;
    kiss_fft_cpx *HyBs2=NULL,*ExBs2=NULL,*HyAs2=NULL,*HxAs2=NULL;
    kiss_fft_cpx *HxBs2=NULL,*EyAs2=NULL,*EyBs2=NULL,*ExAs3=NULL;
    kiss_fft_cpx *HyBs3=NULL,*ExBs3=NULL,*HyAs3=NULL,*HxAs3=NULL;
    kiss_fft_cpx *HxBs3=NULL,*EyAs3=NULL,*EyBs3=NULL,*ExAs4=NULL;
    kiss_fft_cpx *HyBs4=NULL,*ExBs4=NULL,*HyAs4=NULL,*HxAs4=NULL;
    kiss_fft_cpx *HxBs4=NULL,*EyAs4=NULL,*EyBs4=NULL,*ExA1=NULL;
    kiss_fft_cpx *HyB1=NULL,*ExB1=NULL,*HyA1=NULL,*HxA1=NULL,*HxB1=NULL;
    kiss_fft_cpx *EyA1=NULL,*EyB1=NULL,*ExA2=NULL,*HyB2=NULL,*ExB2=NULL;
    kiss_fft_cpx *HyA2=NULL,*HxA2=NULL,*HxB2=NULL,*EyA2=NULL,*EyB2=NULL;
    kiss_fft_cpx *ExA3=NULL,*HyB3=NULL,*ExB3=NULL,*HyA3=NULL,*HxA3=NULL;
    kiss_fft_cpx *HxB3=NULL,*EyA3=NULL,*EyB3=NULL,*ExA4=NULL,*HyB4=NULL;
    kiss_fft_cpx *ExB4=NULL,*HyA4=NULL,*HxA4=NULL,*HxB4=NULL,*EyA4=NULL;
    kiss_fft_cpx *EyB4=NULL,*Zxys1=NULL,*Zyxs1=NULL,*Zxys2=NULL,*Zyxs2=NULL;
    kiss_fft_cpx *Zxys3=NULL,*Zyxs3=NULL,*Zxys4=NULL,*Zyxs4=NULL;
    kiss_fft_cpx *Zxy=NULL,*Zyx=NULL;
    kiss_fftr_cfg cfg;
    sf_init (argc, argv);
    in= sf_input("in");
    Ey=sf_input("Ey");
    Hx=sf_input("Hx");
    Hy=sf_input("Hy");
    out = sf_output("out");
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getbool("comp",&comp)) comp=true;
    /* component selection */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

   if (!sf_getbool("phase",&phase)) phase=false;
    /* if y, calculate apparent resistivity, otherwise calculate phase */

    if (!sf_histdouble(in,"n1",&n1)) sf_error("No n1= in input");
  
    bb= sf_floatalloc(n1);
    pp= sf_floatalloc(n1);
    qq= sf_floatalloc(n1);
    dd= sf_floatalloc(n1);	
   
    sf_floatread(bb,n1,in);
    sf_floatread(pp,n1,Ey);
    sf_floatread(qq,n1,Hx);
    sf_floatread(dd,n1,Hy);
    M=49;
    outp =sf_floatalloc(M);
    pi=3.14;
    f=1; 
    b=-3;

    for(k=0;k<M;k++) {
	a=pow(10,b);
	nw=(int)(128000./(a*1000));
	nt = opt? 2*kiss_fft_next_fast_size((nw+1)/2): nw;
	
	if (nt%2) nt++;
	m1=n1;
	if(m1%nw) N=m1/nw+1;
	else N=m1/nw;
	if(verb) sf_warning("slice %d of %d, freq=%f, stack=%d;",k+1,M,a,N);
	FDEx= (kiss_fft_cpx*)sf_complexalloc(nt);
	FDEy= (kiss_fft_cpx*)sf_complexalloc(nt);
	FDHx= (kiss_fft_cpx*)sf_complexalloc(nt);
	FDHy= (kiss_fft_cpx*)sf_complexalloc(nt);
	TDEx= sf_floatalloc(nt);
	TDEy= sf_floatalloc(nt);
	TDHx= sf_floatalloc(nt);
	TDHy= sf_floatalloc(nt);
	sf_putint(out,"n1",M);
	sf_putfloat(out,"d1",0.1);
	sf_putfloat(out,"o1",-3);
	ExAs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyBs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExBs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyAs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxAs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxBs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyAs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyBs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExAs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyBs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExBs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyAs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxAs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxBs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyAs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyBs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExAs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyBs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExBs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyAs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxAs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxBs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyAs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyBs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExAs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyBs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExBs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyAs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxAs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxBs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyAs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyBs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExA1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyB1=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExB1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyA1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxA1=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxB1=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyA1=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyB1=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExA2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyB2=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExB2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyA2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxA2=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxB2=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyA2=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyB2=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExA3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyB3=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExB3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyA3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxA3=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxB3=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyA3=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyB3=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExA4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyB4=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExB4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HyA4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxA4=(kiss_fft_cpx*)sf_complexalloc(nt);
	HxB4=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyA4=(kiss_fft_cpx*)sf_complexalloc(nt);
	EyB4=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zxys1=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zyxs1=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zxys2=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zyxs2=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zxys3=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zyxs3=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zxys4=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zyxs4=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zxy=(kiss_fft_cpx*)sf_complexalloc(nt);
	Zyx=(kiss_fft_cpx*)sf_complexalloc(nt);
	ExHyres=sf_floatalloc(nt);
	EyHxres=sf_floatalloc(nt);
	ExHypha=sf_floatalloc(nt);
	EyHxpha=sf_floatalloc(nt);
	cfg = kiss_fftr_alloc(nt,0,NULL,NULL);
	
	for(j=0;j<N;j++) {
	    for(i=0;i<nt;i++){
		if((i<nw)&&((j*nw+i)<n1)) {
		    TDEx[i]=bb[j*nw+i]*0.5*(1-cos(2*pi*i/nt));
		}
		else {
		    TDEx[i]=0.;
		}
	    }
	    kiss_fftr(cfg,TDEx,FDEx);
	    
	    
	    for(i=0;i<nt;i++){
		if((i<nw)&&((j*nw+i)<n1)) {
		    TDEy[i]=pp[j*nw+i]*0.5*(1-cos(2*pi*i/nt));
		}
		else{
		    TDEy[i]=0.;
		} 
	    }
	    kiss_fftr(cfg,TDEy,FDEy);  
	    
	    
	    for(i=0;i<nt;i++){
		if((i<nw)&&((j*nw+i)<n1)) {
		    TDHx[i]=qq[j*nw+i]*0.5*(1-cos(2*pi*i/nt));
		}
		else{TDHx[i]=0.;} 
	    }
	    kiss_fftr(cfg,TDHx,FDHx); 
	    
	    
	    for(i=0;i<nt;i++){
		if((i<nw)&&((j*nw+i)<n1)) {
		    TDHy[i]=dd[j*nw+i]*0.5*(1-cos(2*pi*i/nt));
		}
		else{
		    TDHy[i]=0.;
		} 
	    }
	    kiss_fftr(cfg,TDHy,FDHy); 
	    A=FDEx;B=FDEy;
	    ExAs1[f]=sf_cadd(ExAs1[f],sf_cmul(FDEx[f],sf_conjf(A[f])));
	    HyBs1[f]=sf_cadd(HyBs1[f],sf_cmul(FDHy[f],sf_conjf(B[f])));
	    ExBs1[f]=sf_cadd(ExBs1[f],sf_cmul(FDEx[f],sf_conjf(B[f])));
	    HyAs1[f]=sf_cadd(HyAs1[f],sf_cmul(FDHy[f],sf_conjf(A[f])));
	    HxAs1[f]=sf_cadd(HxAs1[f],sf_cmul(FDHx[f],sf_conjf(A[f])));
	    HxBs1[f]=sf_cadd(HxBs1[f],sf_cmul(FDHx[f],sf_conjf(B[f])));
	    EyAs1[f]=sf_cadd(EyAs1[f],sf_cmul(FDEy[f],sf_conjf(A[f])));
	    EyBs1[f]=sf_cadd(EyBs1[f],sf_cmul(FDEy[f],sf_conjf(B[f])));
	  
	    A=FDEx;B=FDHx;
	    ExAs2[f]= sf_cadd(ExAs2[f],sf_cmul(FDEx[f],sf_conjf(A[f])));
	    HyBs2[f]= sf_cadd(HyBs2[f],sf_cmul(FDHy[f],sf_conjf(B[f])));
	    ExBs2[f]= sf_cadd(ExBs2[f],sf_cmul(FDEx[f],sf_conjf(B[f])));
	    HyAs2[f]= sf_cadd(HyAs2[f],sf_cmul(FDHy[f],sf_conjf(A[f])));
	    HxAs2[f]= sf_cadd(HxAs2[f],sf_cmul(FDHx[f],sf_conjf(A[f])));
	    HxBs2[f]= sf_cadd(HxBs2[f],sf_cmul(FDHx[f],sf_conjf(B[f])));
	    EyAs2[f]= sf_cadd(EyAs2[f],sf_cmul(FDEy[f],sf_conjf(A[f])));
	    EyBs2[f]= sf_cadd(EyBs2[f],sf_cmul(FDEy[f],sf_conjf(B[f])));
	    
	   
	    A=FDEy;B=FDHy;
	    ExAs3[f]= sf_cadd(ExAs3[f],sf_cmul(FDEx[f],sf_conjf(A[f])));
	    HyBs3[f]= sf_cadd(HyBs3[f],sf_cmul(FDHy[f],sf_conjf(B[f])));
	    ExBs3[f]= sf_cadd(ExBs3[f],sf_cmul(FDEx[f],sf_conjf(B[f])));
	    HyAs3[f]= sf_cadd(HyAs3[f],sf_cmul(FDHy[f],sf_conjf(A[f])));
	    HxAs3[f]= sf_cadd(HxAs3[f],sf_cmul(FDHx[f],sf_conjf(A[f])));
	    HxBs3[f]= sf_cadd(HxBs3[f],sf_cmul(FDHx[f],sf_conjf(B[f])));
	    EyAs3[f]= sf_cadd(EyAs3[f],sf_cmul(FDEy[f],sf_conjf(A[f])));
	    EyBs3[f]= sf_cadd(EyBs3[f],sf_cmul(FDEy[f],sf_conjf(B[f])));
	   
	 
	    A=FDHx;B=FDHy;
	    ExAs4[f]= sf_cadd(ExAs4[f],sf_cmul(FDEx[f],sf_conjf(A[f])));
	    HyBs4[f]= sf_cadd(HyBs4[f],sf_cmul(FDHy[f],sf_conjf(B[f])));
	    ExBs4[f]= sf_cadd(ExBs4[f],sf_cmul(FDEx[f],sf_conjf(B[f])));
	    HyAs4[f]= sf_cadd(HyAs4[f],sf_cmul(FDHy[f],sf_conjf(A[f])));
	    HxAs4[f]= sf_cadd(HxAs4[f],sf_cmul(FDHx[f],sf_conjf(A[f])));
	    HxBs4[f]= sf_cadd(HxBs4[f],sf_cmul(FDHx[f],sf_conjf(B[f])));
	    EyAs4[f]= sf_cadd(EyAs4[f],sf_cmul(FDEy[f],sf_conjf(A[f])));
	    EyBs4[f]= sf_cadd(EyBs4[f],sf_cmul(FDEy[f],sf_conjf(B[f]))); 
	}
  
        ExA1[f]=sf_crmul(ExAs1[f],1./N);
        HyB1[f]=sf_crmul(HyBs1[f],1./N);
        ExB1[f]=sf_crmul(ExBs1[f],1./N);
	HyA1[f]=sf_crmul(HyAs1[f],1./N);
	HxA1[f]=sf_crmul(HxAs1[f],1./N);
	HxB1[f]=sf_crmul(HxBs1[f],1./N);
	EyA1[f]=sf_crmul(EyAs1[f],1./N);
	EyB1[f]=sf_crmul(EyBs1[f],1./N);
	ExA2[f]=sf_crmul(ExAs2[f],1./N);
	HyB2[f]=sf_crmul(HyBs2[f],1./N);
	ExB2[f]=sf_crmul(ExBs2[f],1./N);
	HyA2[f]=sf_crmul(HyAs2[f],1./N);
	HxA2[f]=sf_crmul(HxAs2[f],1./N);
	HxB2[f]=sf_crmul(HxBs2[f],1./N);
	EyA2[f]=sf_crmul(EyAs2[f],1./N);
	EyB2[f]=sf_crmul(EyBs2[f],1./N);
	ExA3[f]=sf_crmul(ExAs3[f],1./N);
	HyB3[f]=sf_crmul(HyBs3[f],1./N);
	ExB3[f]=sf_crmul(ExBs3[f],1./N);
	HyA3[f]=sf_crmul(HyAs3[f],1./N);
	HxA3[f]=sf_crmul(HxAs3[f],1./N);
	HxB3[f]=sf_crmul(HxBs3[f],1./N);
	EyA3[f]=sf_crmul(EyAs3[f],1./N);
	EyB3[f]=sf_crmul(EyBs3[f],1./N);
	ExA4[f]=sf_crmul(ExAs4[f],1./N);
	HyB4[f]=sf_crmul(HyBs4[f],1./N);
	ExB4[f]=sf_crmul(ExBs4[f],1./N);
	HyA4[f]=sf_crmul(HyAs4[f],1./N);
	HxA4[f]=sf_crmul(HxAs4[f],1./N);
	HxB4[f]=sf_crmul(HxBs4[f],1./N);
	EyA4[f]=sf_crmul(EyAs4[f],1./N);
	EyB4[f]=sf_crmul(EyBs4[f],1./N);	  

	Zxys1[f]=sf_cdiv((sf_csub(sf_cmul(ExA1[f],HxB1[f]),
				  sf_cmul(ExB1[f],HxA1[f]))),
			 (sf_csub(sf_cmul(HyA1[f],HxB1[f]),
				  sf_cmul(HyB1[f],HxA1[f]))));
	Zyxs1[f]=sf_cdiv((sf_csub(sf_cmul(EyA1[f],HyB1[f]),
				  sf_cmul(EyB1[f],HyA1[f]))),
			 (sf_csub(sf_cmul(HxA1[f],HyB1[f]),
				  sf_cmul(HxB1[f],HyA1[f]))));
	Zxys2[f]=sf_cdiv((sf_csub(sf_cmul(ExA2[f],HxB2[f]),
				  sf_cmul(ExB2[f],HxA2[f]))),
			 (sf_csub(sf_cmul(HyA2[f],HxB2[f]),
				  sf_cmul(HyB2[f],HxA2[f]))));
	Zyxs2[f]=sf_cdiv((sf_csub(sf_cmul(EyA2[f],HyB2[f]),
				  sf_cmul(EyB2[f],HyA2[f]))),
			 (sf_csub(sf_cmul(HxA2[f],HyB2[f]),
				  sf_cmul(HxB2[f],HyA2[f]))));
	Zxys3[f]=sf_cdiv((sf_csub(sf_cmul(ExA3[f],HxB3[f]),
				  sf_cmul(ExB3[f],HxA3[f]))),
			 (sf_csub(sf_cmul(HyA3[f],HxB3[f]),
				  sf_cmul(HyB3[f],HxA3[f]))));
	Zyxs3[f]=sf_cdiv((sf_csub(sf_cmul(EyA3[f],HyB3[f]),
				  sf_cmul(EyB3[f],HyA3[f]))),
			 (sf_csub(sf_cmul(HxA3[f],HyB3[f]),
				  sf_cmul(HxB3[f],HyA3[f]))));
	Zxys4[f]=sf_cdiv((sf_csub(sf_cmul(ExA4[f],HxB4[f]),
				  sf_cmul(ExB4[f],HxA4[f]))),
			 (sf_csub(sf_cmul(HyA4[f],HxB4[f]),
				  sf_cmul(HyB4[f],HxA4[f]))));
	Zyxs4[f]=sf_cdiv((sf_csub(sf_cmul(EyA4[f],HyB4[f]),
				  sf_cmul(EyB4[f],HyA4[f]))),
			 (sf_csub(sf_cmul(HxA4[f],HyB4[f]),
				  sf_cmul(HxB4[f],HyA4[f]))));

 
     
       Zxy[f]=sf_crmul(sf_cadd(sf_cadd(Zxys1[f],Zxys2[f]),
			       sf_cadd(Zxys3[f],Zxys4[f])),0.25);
       Zyx[f]=sf_crmul(sf_cadd(sf_cadd(Zyxs1[f],Zyxs2[f]),
			       sf_cadd(Zyxs3[f],Zyxs4[f])),0.25);
       ExHyres[f]=0.2*sf_cabsf(Zxy[f])*sf_cabsf(Zxy[f])/a;
       EyHxres[f]=0.2*sf_cabsf(Zyx[f])*sf_cabsf(Zyx[f])/a;
       ExHypha[f]=sf_cargf(Zxy[f]);
       EyHxpha[f]=sf_cargf(Zyx[f]);
       
       if(phase) {
	   if(comp) {
	       outp[k]= ExHypha[f];
	   } else {
	       outp[k]= EyHxpha[f];
	   }
	   
       } else {
	   if(comp) {
	       outp[k]= ExHyres[f];
	   } else {
	       outp[k]= EyHxres[f];
	   }
       }
       b=b+0.1;
       
    }
      
    sf_floatwrite(outp,M,out);
        
    exit(0);
}

/* 	$Id$	 */

