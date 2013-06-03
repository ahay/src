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
    bool opt, comp, phase;
    double n1;
    float nw1, b, a, dw, pi;
    int nt, nw, i, j, N, m1, M, k, f;                    
    float *TDEx,*TDEy,*TDHx,*TDHy,*outp,*ExHyres,*EyHxres,*ExHypha,*EyHxpha;
    float *bb,*pp,*qq,*dd;
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
    sf_file in, out,Ey,Hx,Hy;
    in= sf_input("in");
    Ey=sf_input("Ey");
    Hx=sf_input("Hx");
    Hy=sf_input("Hy");
    out = sf_output("out");
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getbool("comp",&comp)) comp=true;
    /* component selection */

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
    M=48;
    outp =sf_floatalloc(M);
    pi=3.14;
    f=1; 
    b=-3;
    for(k=0;k<M;k++) {
	a=pow(10,b);
	nw=(int)(128000./(a*1000));
	nt = opt? 2*kiss_fft_next_fast_size((nw+1)/2): nw;
	
	if (nt%2) nt++;
	nw1 = nt/2+1;
	dw = 128./nt; 
/*	sf_warning("a=%f,dw=%f\n",a,dw); */ 
	m1=n1;
	if(m1%nw) m1++;
	N=m1/nw;
/* 	f=1; */
/* 	f=(int)((float)a/(float)dw); */
/*  	sf_warning("nw=%d,nt=%d\n",nw,nt); */ 
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
	    /*  if (j%100000==0) sf_warning("slice %d of %d;",j+1,N); */
	    
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
		else {
		    TDEy[i]=0.;
		} 
	    }
	    kiss_fftr(cfg,TDEy,FDEy);  
	    
	    
	    for(i=0;i<nt;i++){
		if((i<nw)&&((j*nw+i)<n1)) {
		    TDHx[i]=qq[j*nw+i]*0.5*(1-cos(2*pi*i/nt));
		}
		else {
		    TDHx[i]=0.;
		} 
	    }
	    kiss_fftr(cfg,TDHx,FDHx); 
	    
	    
	    for(i=0;i<nt;i++){
		if((i<nw)&&((j*nw+i)<n1)) {
		    TDHy[i]=dd[j*nw+i]*0.5*(1-cos(2*pi*i/nt));
		}
		else {
		    TDHy[i]=0.;
		} 
	    }
	    kiss_fftr(cfg,TDHy,FDHy); 
	    /*  for(i=0;i<nt;i++){ */
/* 	    if((i*dw<=a)&&((i+1)*dw>a)) {f=i;} */
/* 	    } */
	    A=FDEx;B=FDEy;
	    ExAs1[f]=sf_cmul(FDEx[f],sf_conjf(A[f]));
	    HyBs1[f]=sf_cmul(FDHy[f],sf_conjf(B[f]));
	    ExBs1[f]=sf_cmul(FDEx[f],sf_conjf(B[f]));
	    HyAs1[f]=sf_cmul(FDHy[f],sf_conjf(A[f]));
	    HxAs1[f]=sf_cmul(FDHx[f],sf_conjf(A[f]));
	    HxBs1[f]=sf_cmul(FDHx[f],sf_conjf(B[f]));
	    EyAs1[f]=sf_cmul(FDEy[f],sf_conjf(A[f]));
	    EyBs1[f]=sf_cmul(FDEy[f],sf_conjf(B[f]));
	    
	    
	    A=FDEx;B=FDHx;
	    ExAs2[f]=sf_cmul(FDEx[f],sf_conjf(A[f]));
	    HyBs2[f]=sf_cmul(FDHy[f],sf_conjf(B[f]));
	    ExBs2[f]=sf_cmul(FDEx[f],sf_conjf(B[f]));
	    HyAs2[f]=sf_cmul(FDHy[f],sf_conjf(A[f]));
	    HxAs2[f]=sf_cmul(FDHx[f],sf_conjf(A[f]));
	    HxBs2[f]=sf_cmul(FDHx[f],sf_conjf(B[f]));
	    EyAs2[f]=sf_cmul(FDEy[f],sf_conjf(A[f]));
	    EyBs2[f]=sf_cmul(FDEy[f],sf_conjf(B[f]));
	    
	    
	    A=FDEy;B=FDHy;
	    ExAs3[f]=sf_cmul(FDEx[f],sf_conjf(A[f]));
	    HyBs3[f]=sf_cmul(FDHy[f],sf_conjf(B[f]));
	    ExBs3[f]=sf_cmul(FDEx[f],sf_conjf(B[f]));
	    HyAs3[f]=sf_cmul(FDHy[f],sf_conjf(A[f]));
	    HxAs3[f]=sf_cmul(FDHx[f],sf_conjf(A[f]));
	    HxBs3[f]=sf_cmul(FDHx[f],sf_conjf(B[f]));
	    EyAs3[f]=sf_cmul(FDEy[f],sf_conjf(A[f]));
	    EyBs3[f]=sf_cmul(FDEy[f],sf_conjf(B[f]));
	    
	    
	    A=FDHx;B=FDHy;
	    ExAs4[f]=sf_cmul(FDEx[f],sf_conjf(A[f]));
	    HyBs4[f]=sf_cmul(FDHy[f],sf_conjf(B[f]));
	    ExBs4[f]=sf_cmul(FDEx[f],sf_conjf(B[f]));
	    HyAs4[f]=sf_cmul(FDHy[f],sf_conjf(A[f]));
	    HxAs4[f]=sf_cmul(FDHx[f],sf_conjf(A[f]));
	    HxBs4[f]=sf_cmul(FDHx[f],sf_conjf(B[f]));
	    EyAs4[f]=sf_cmul(FDEy[f],sf_conjf(A[f]));
	    EyBs4[f]=sf_cmul(FDEy[f],sf_conjf(B[f]));
	    
	    
	    ExA1[k]=sf_cadd(ExA1[k],ExAs1[f]);
	    HyB1[k]=sf_cadd(HyB1[k],HyBs1[f]);
	    ExB1[k]=sf_cadd(ExB1[k],ExBs1[f]);
	    HyA1[k]=sf_cadd(HyA1[k],HyAs1[f]);
	    HxA1[k]=sf_cadd(HxA1[k],HxAs1[f]);
	    HxB1[k]=sf_cadd(HxB1[k],HxBs1[f]);
	    EyA1[k]=sf_cadd(EyA1[k],EyAs1[f]);
	    EyB1[k]=sf_cadd(EyB1[k],EyBs1[f]);
	    
	    ExA2[k]=sf_cadd(ExA2[k],ExAs2[f]);
	    HyB2[k]=sf_cadd(HyB2[k],HyBs2[f]);
	    ExB2[k]=sf_cadd(ExB2[k],ExBs2[f]);
	    HyA2[k]=sf_cadd(HyA2[k],HyAs2[f]);
	    HxA2[k]=sf_cadd(HxA2[k],HxAs2[f]);
	    HxB2[k]=sf_cadd(HxB2[k],HxBs2[f]);
	    EyA2[k]=sf_cadd(EyA2[k],EyAs2[f]);
	    EyB2[k]=sf_cadd(EyB2[k],EyBs2[f]);
	    
	    
	    ExA3[k]=sf_cadd(ExA3[k],ExAs3[f]);
	    HyB3[k]=sf_cadd(HyB3[k],HyBs3[f]);
	    ExB3[k]=sf_cadd(ExB3[k],ExBs3[f]);
	    HyA3[k]=sf_cadd(HyA3[k],HyAs3[f]);
	    HxA3[k]=sf_cadd(HxA3[k],HxAs3[f]);
	    HxB3[k]=sf_cadd(HxB3[k],HxBs3[f]);
	    EyA3[k]=sf_cadd(EyA3[k],EyAs3[f]);
	    EyB3[k]=sf_cadd(EyB3[k],EyBs3[f]);
	    
	    
	    ExA4[k]=sf_cadd(ExA4[k],ExAs4[f]);
	    HyB4[k]=sf_cadd(HyB4[k],HyBs4[f]);
	    ExB4[k]=sf_cadd(ExB4[k],ExBs4[f]);
	    HyA4[k]=sf_cadd(HyA4[k],HyAs4[f]);
	    HxA4[k]=sf_cadd(HxA4[k],HxAs4[f]);
	    HxB4[k]=sf_cadd(HxB4[k],HxBs4[f]);
	    EyA4[k]=sf_cadd(EyA4[k],EyAs4[f]);
	    EyB4[k]=sf_cadd(EyB4[k],EyBs4[f]);
	    
	}
	
        ExA1[k]=sf_crmul(ExA1[k],1./N);
        HyB1[k]=sf_crmul(HyB1[k],1./N);
        ExB1[k]=sf_crmul(ExB1[k],1./N);
	HyA1[k]=sf_crmul(HyA1[k],1./N);
	HxA1[k]=sf_crmul(HxA1[k],1./N);
	HxB1[k]=sf_crmul(HxB1[k],1./N);
	EyA1[k]=sf_crmul(EyA1[k],1./N);
	EyB1[k]=sf_crmul(EyB1[k],1./N);
	ExA2[k]=sf_crmul(ExA2[k],1./N);
	HyB2[k]=sf_crmul(HyB2[k],1./N);
	ExB2[k]=sf_crmul(ExB2[k],1./N);
	HyA2[k]=sf_crmul(HyA2[k],1./N);
	HxA2[k]=sf_crmul(HxA2[k],1./N);
	HxB2[k]=sf_crmul(HxB2[k],1./N);
	EyA2[k]=sf_crmul(EyA2[k],1./N);
	EyB2[k]=sf_crmul(EyB2[k],1./N);
	ExA3[k]=sf_crmul(ExA3[k],1./N);
	HyB3[k]=sf_crmul(HyB3[k],1./N);
	ExB3[k]=sf_crmul(ExB3[k],1./N);
	HyA3[k]=sf_crmul(HyA3[k],1./N);
	HxA3[k]=sf_crmul(HxA3[k],1./N);
	HxB3[k]=sf_crmul(HxB3[k],1./N);
	EyA3[k]=sf_crmul(EyA3[k],1./N);
	EyB3[k]=sf_crmul(EyB3[k],1./N);
	ExA4[k]=sf_crmul(ExA4[k],1./N);
	HyB4[k]=sf_crmul(HyB4[k],1./N);
	ExB4[k]=sf_crmul(ExB4[k],1./N);
	HyA4[k]=sf_crmul(HyA4[k],1./N);
	HxA4[k]=sf_crmul(HxA4[k],1./N);
	HxB4[k]=sf_crmul(HxB4[k],1./N);
	EyA4[k]=sf_crmul(EyA4[k],1./N);
	EyB4[k]=sf_crmul(EyB4[k],1./N);	  

	Zxys1[k]=sf_cdiv((sf_csub(sf_cmul(ExA1[k],HxB1[k]),
				  sf_cmul(ExB1[k],HxA1[k]))),
			 (sf_csub(sf_cmul(HyA1[k],HxB1[k]),
				  sf_cmul(HyB1[k],HxA1[k]))));
	Zyxs1[k]=sf_cdiv((sf_csub(sf_cmul(EyA1[k],HyB1[k]),
				  sf_cmul(EyB1[k],HyA1[k]))),
			 (sf_csub(sf_cmul(HxA1[k],HyB1[k]),
				  sf_cmul(HxB1[k],HyA1[k]))));
	Zxys2[k]=sf_cdiv((sf_csub(sf_cmul(ExA2[k],HxB2[k]),
				  sf_cmul(ExB2[k],HxA2[k]))),
			 (sf_csub(sf_cmul(HyA2[k],HxB2[k]),
				  sf_cmul(HyB2[k],HxA2[k]))));
	Zyxs2[k]=sf_cdiv((sf_csub(sf_cmul(EyA2[k],HyB2[k]),
				  sf_cmul(EyB2[k],HyA2[k]))),
			 (sf_csub(sf_cmul(HxA2[k],HyB2[k]),
				  sf_cmul(HxB2[k],HyA2[k]))));
	Zxys3[k]=sf_cdiv((sf_csub(sf_cmul(ExA3[k],HxB3[k]),
				  sf_cmul(ExB3[k],HxA3[k]))),
			 (sf_csub(sf_cmul(HyA3[k],HxB3[k]),
				  sf_cmul(HyB3[k],HxA3[k]))));
	Zyxs3[k]=sf_cdiv((sf_csub(sf_cmul(EyA3[k],HyB3[k]),
				  sf_cmul(EyB3[k],HyA3[k]))),
			 (sf_csub(sf_cmul(HxA3[k],HyB3[k]),
				  sf_cmul(HxB3[k],HyA3[k]))));
	Zxys4[k]=sf_cdiv((sf_csub(sf_cmul(ExA4[k],HxB4[k]),
				  sf_cmul(ExB4[k],HxA4[k]))),
			 (sf_csub(sf_cmul(HyA4[k],HxB4[k]),
				  sf_cmul(HyB4[k],HxA4[k]))));
	Zyxs4[k]=sf_cdiv((sf_csub(sf_cmul(EyA4[k],HyB4[k]),
				  sf_cmul(EyB4[k],HyA4[k]))),
			 (sf_csub(sf_cmul(HxA4[k],HyB4[k]),
				  sf_cmul(HxB4[k],HyA4[k]))));
     
	Zxy[k]=sf_crmul(sf_cadd(sf_cadd(Zxys1[k],Zxys2[k]),
				sf_cadd(Zxys3[k],Zxys4[k])),0.25);
	Zyx[k]=sf_crmul(sf_cadd(sf_cadd(Zyxs1[k],Zyxs2[k]),
				sf_cadd(Zyxs3[k],Zyxs4[k])),0.25);
	ExHyres[k]=0.2*sf_cabsf(Zxy[k])*sf_cabsf(Zxy[k])/a;
	EyHxres[k]=0.2*sf_cabsf(Zyx[k])*sf_cabsf(Zyx[k])/a;
	ExHypha[k]=sf_cargf(Zxy[k]);
	EyHxpha[k]=sf_cargf(Zyx[k]);

	if(phase) {
	    if(comp) {
		outp[k]= ExHypha[k];
	    } else {
		outp[k]= EyHxpha[k];
	    }
	    
	} else {
	    if(comp) {
		outp[k]= ExHyres[k];
	    } else {
		outp[k]= EyHxres[k];
	    }
	}
	b=b+0.1;
    } 
        
    sf_floatwrite(outp,M,out);
        
    exit(0);
}

/* 	$Id$	 */

