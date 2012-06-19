/* Stationary Wavelet Transform */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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
# include <rsf.h>

/* db1 wavelet is used. other db wavelet can also be integrated. */

#include "swt.h"

static void conv(float p[],int m,float q[],int n,float r[])
{
	int k=m+n-1;
	int i,j;
    for (i=0; i<k; i++) r[i]=0.0;
    for (i=0; i<m; i++)
    for (j=0; j<n; j++) r[i+j]=r[i+j]+p[i]*q[j];
}

static void dyaddown(float x1[],int n,float x2[])
{
    int i,k;
	for(i=1,k=0;i<n;i+=2,k++) x2[k]=x1[i];

}

static void dyaddownodd(float x1[],int n,float x2[])
{
    int i,k;
	for(i=0,k=0;i<n;i+=2,k++) x2[k]=x1[i];

}

static void dyadup(float x1[],int n,float x2[])
{
	int i,k;
	for(i=0,k=1;i<n;i++,k+=2)
	{
		x2[k]=x1[i];
		x2[k-1]=0;
	}
	x2[k-1]=0;
}

static void dyadupodd(float x1[],int n,float x2[])
{
	int i,k;
	for(i=0,k=1;i<n;i++,k+=2)
	{
		x2[k-1]=x1[i];
		x2[k]=0;
	}
	x2[k]=0;
}

static void wkeep(float x1[],int n,int l,float x2[])
{
	int m,i,k;
	m=(n-l)/2;
	for(i=m,k=0;i<n-m;i++,k++) x2[k]=x1[i];
}



void fwavelet(float s[]    /* input signal [n_signal] */,
	      int n_signal /* signal length */,
	      float wcl[]  /* output low-frequency part  [n_signal] */,
	      float wch[]  /* output high-frequency part [n_signal] */)   
/*< forward SWT >*/
{
  	int i,len_filter,n_len;

	float *ca1,*cd1,*tempo,*caodd1,*cdodd1;

	float Lo_D[]={0.70711,0.70711};
	float Hi_D[]={-0.70711,0.70711};

	len_filter=2;
	n_len=n_signal+len_filter-1;

	ca1=sf_floatalloc(n_signal/2);
	cd1=sf_floatalloc(n_signal/2);
	tempo=sf_floatalloc(n_signal+len_filter-1);
	caodd1=sf_floatalloc(n_signal/2);
	cdodd1=sf_floatalloc(n_signal/2);

	conv(s,n_signal,Lo_D,len_filter,tempo);   /* the even low frequency */
	dyaddown(tempo,n_len,ca1);
	dyaddownodd(tempo,n_len,caodd1);

	conv(s,n_signal,Hi_D,len_filter,tempo);   /* the even high frequency */
	dyaddown(tempo,n_len,cd1);
	dyaddownodd(tempo,n_len,cdodd1);

	for(i=0;i<n_signal/2;i++)
	{
	  wcl[2*i]=caodd1[i];
	  wcl[2*i+1]=ca1[i];      

	}

	for(i=0;i<n_signal/2;i++)
	{
	  wch[2*i]=cdodd1[i];
	  wch[2*i+1]=cd1[i];      

	}

}

void iwavelet(float wcl[]  /* input low-frequncy part [n_signal] */,
	      float wch[]  /* input high-frequency part [n_signal] */,
	      int n_signal /* signal length */,
	      float s[]    /* output signal [n_signal] */) 
/*< inverse transform >*/
{
    int i,len_filter,n_len;
	float *ca1,*cd1,*tempo,*caodd1,*cdodd1;
	float *tempo2,*d1,*a1,*a0,*d2,*a2,*a00;

    float Lo_R[]={0.70711,0.70711};
	float Hi_R[]={0.70711,-0.70711};

	len_filter=2;
	n_len=n_signal+len_filter-1;

	ca1=sf_floatalloc(n_signal/2);
	cd1=sf_floatalloc(n_signal/2);
	tempo=sf_floatalloc((n_signal+len_filter-1));
	caodd1=sf_floatalloc(n_signal/2);
	cdodd1=sf_floatalloc(n_signal/2);

	tempo2=sf_floatalloc((n_signal+len_filter));
	d1=sf_floatalloc(n_signal);
	a1=sf_floatalloc(n_signal);
	a0=sf_floatalloc(n_signal);
	d2=sf_floatalloc(n_signal);
	a2=sf_floatalloc(n_signal);
	a00=sf_floatalloc(n_signal);	
	
	for(i=0;i<n_signal/2;i++)
	{
	  caodd1[i]=wcl[2*i];
	  ca1[i]=wcl[2*i+1];      

	}

	for(i=0;i<n_signal/2;i++)
	{
	  cdodd1[i]=wch[2*i];
	  cd1[i]=wch[2*i+1];      

	}

	dyadup(cd1,n_signal/2,tempo);                      /* The even */
	conv(tempo,n_len,Hi_R,len_filter,tempo2);
	wkeep(tempo2,n_len+len_filter-1,n_signal,d1);

	dyadup(ca1,n_signal/2,tempo);
	conv(tempo,n_len,Lo_R,len_filter,tempo2);
	wkeep(tempo2,n_len+len_filter-1,n_signal,a1);


	dyadupodd(cdodd1,n_signal/2,tempo);                /* the odd */
	conv(tempo,n_len,Hi_R,len_filter,tempo2);
	wkeep(tempo2,n_len+len_filter-1,n_signal,d2);

	dyadupodd(caodd1,n_signal/2,tempo);
	conv(tempo,n_len,Lo_R,len_filter,tempo2);
	wkeep(tempo2,n_len+len_filter-1,n_signal,a2);  

	for(i=0;i<n_signal;i++) 
		a0[i]=a1[i]+d1[i];

	for(i=0;i<n_signal;i++) 
		a00[i]=a2[i]+d2[i];

    s[0]=a0[0];
    for(i=1;i<n_signal-1;i++) 
		s[i]=(a0[i]+a00[i])/2;
	s[n_signal-1]=a0[n_signal-1];

}

void  multi_fwavelet(float ss[]   /* input signal [n_signal] */,
		     int n_signal /* signal length */,
		     int n_layer  /* layer number */,
		     float cl[]   /* low-frequency [n_signal] */,
		     float ch[]   /* high-frequency  [n_signal*n_layer] */) 
/*<  forward multi-layer SWT >*/
{
    int i,j;
	float *wcl,*wch;

	wcl=sf_floatalloc(n_signal);
	wch=sf_floatalloc(n_signal);

	for(i=0;i<n_layer;i++)
	{
      fwavelet(ss,n_signal,wcl,wch);

	  for(j=0;j<n_signal;j++)
          ch[i*n_signal+j]=wch[j];

	  for(j=0;j<n_signal;j++)
          ss[j]=wcl[j];
	}
    
	for(i=0;i<n_signal;i++)
          cl[i]=wcl[i];

	free(wcl);
	free(wch);
}

void  multi_iwavelet(float cl[]       /* low-frequency [n_signal] */,
		     float ch[]       /* high-frequency  [n_signal*n_layer] */,
		     int n_signal     /* signal length */,
		     int n_layer      /* layer number */,
		     float ainverse[] /* output signal [n_signal] */)
/*<  inverse multi-layer SWT >*/
{
	
	int i,j;
	float *wcl,*wch;

	wcl=sf_floatalloc(n_signal);
	wch=sf_floatalloc(n_signal);
    
	for(i=0;i<n_signal;i++)
          wcl[i]=cl[i];

	for(i=n_layer-1;i>=0;i--)
	{     
	  for(j=0;j<n_signal;j++)
          wch[j]=ch[i*n_signal+j];

	  iwavelet(wcl,wch,n_signal,ainverse);

	  for(j=0;j<n_signal;j++)
          wcl[j]=ainverse[j];
	}

	free(wcl);
	free(wch);
 
}

