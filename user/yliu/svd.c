/* matrix SVD decomposition and reconstruction */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include <stdio.h>
#include <math.h>
#include <rsf.h>
#include "svd.h"

static int m, n;
static float *s, *e, *w;
static double eps;

void svdinit( int n2, int n1, int ka, double epsilon)
/*< initiate svd and allocate memory >*/
{
    m=n2;
    n=n1;
    s = sf_floatalloc(ka); 
    e = sf_floatalloc(ka);
    w = sf_floatalloc(ka);
    eps = epsilon;

}

void svdclose( void )
/*< release memory >*/
{
    free(s);
    free(e);
    free(w);
}

void brmul(const float *a, const float *b,int m,int n,int k,float *c) 
/*<SVD reconstruction>*/
{
    int i,j,l,u;
    for (i=0; i<=m-1; i++)
	for (j=0; j<=k-1; j++)	{ 
	    u=i*k+j; c[u]=0.0;
	    for (l=0; l<=n-1; l++)
		c[u]=c[u]+a[i*n+l]*b[l*k+j];
	}
    return;
}

static void ppp(float *a,float *e,float *s,float *v,int m, int n)
{
    int i,j,p,q;
    float d;
    if (m>=n) i=n;
    else i=m;
    for (j=1; j<=i-1; j++) {
	a[(j-1)*n+j-1]=s[j-1];
        a[(j-1)*n+j]=e[j-1];
    }
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n) a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
	for (j=i+1; j<=n; j++) {
	    p=(i-1)*n+j-1; q=(j-1)*n+i-1;
	    d=v[p]; v[p]=v[q]; v[q]=d;
	}
    return;
}

static void sss(float fg[2], float cs[2] )
{ 
     float r,d;
     if ((fabsf(fg[0])+fabsf(fg[1]))==0.0)  {
	 cs[0]=1.0; cs[1]=0.0; d=0.0;}
     else  {
	 d=sqrtf(fg[0]*fg[0]+fg[1]*fg[1]);
	 if (fabsf(fg[0])>fabsf(fg[1])) {
	     d=fabsf(d);
	     if (fg[0]<0.0) d=-d;
	 }
	 if (fabsf(fg[1])>=fabsf(fg[0])) {
	     d=fabsf(d);
	     if (fg[1]<0.0) d=-d;
	 }
	 cs[0]=fg[0]/d; cs[1]=fg[1]/d;
     }
    r=1.0;
    if (fabsf(fg[0])>fabsf(fg[1])) r=cs[1];
    else
	if (cs[0]!=0.0) r=1.0/cs[0];
    fg[0]=d; fg[1]=r;
    return;
}



int svduav(float *a,float *u,float *v) 
/*<SVD decomposition, A->UAV >*/
{
    int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
    float d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];

    it=60; k=n;
    if (m-1<n) k=m-1;
    l=m;
    if (n-2<m) l=n-2;
    if (l<0) l=0;
    ll=k;
    if (l>k) ll=l;
    if (ll>=1) {
	for (kk=1; kk<=ll; kk++) {
	    if (kk<=k) {
		d=0.0;
                for (i=kk; i<=m; i++) {
		    ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];
		}
                s[kk-1]=sqrtf(d);
                if (s[kk-1]!=0.0) {
		    ix=(kk-1)*n+kk-1;
                    if (a[ix]!=0.0) {
			s[kk-1]=fabsf(s[kk-1]);
                        if (a[ix]<0.0) s[kk-1]=-s[kk-1];
		    }
                    for (i=kk; i<=m; i++) {
			iy=(i-1)*n+kk-1;
                        a[iy]=a[iy]/s[kk-1];
		    }
                    a[ix]=1.0+a[ix];
		}
                s[kk-1]=-s[kk-1];
	    }
            if (n>=kk+1) {
		for (j=kk+1; j<=n; j++) {
		    if ((kk<=k)&&(s[kk-1]!=0.0)) {
			d=0.0;
                        for (i=kk; i<=m; i++) {
			    ix=(i-1)*n+kk-1;
                            iy=(i-1)*n+j-1;
                            d=d+a[ix]*a[iy];
			}
                        d=-d/a[(kk-1)*n+kk-1];
                        for (i=kk; i<=m; i++) {
			    ix=(i-1)*n+j-1;
                            iy=(i-1)*n+kk-1;
                            a[ix]=a[ix]+d*a[iy];
			}
		    }
                    e[j-1]=a[(kk-1)*n+j-1];
		}
	    }
            if (kk<=k) {
		for (i=kk; i<=m; i++) {
		    ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
		}
	    }
            if (kk<=l) {
		d=0.0;
                for (i=kk+1; i<=n; i++)
		    d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrtf(d);
                if (e[kk-1]!=0.0) {
		    if (e[kk]!=0.0) {
			e[kk-1]=fabsf(e[kk-1]);
                        if (e[kk]<0.0) e[kk-1]=-e[kk-1];
		    }
                    for (i=kk+1; i<=n; i++)
			e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
		}
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0)) {
		    for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                    for (j=kk+1; j<=n; j++)
			for (i=kk+1; i<=m; i++)
			    w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
			for (i=kk+1; i<=m; i++) {
			    ix=(i-1)*n+j-1;
			    a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
		}
                for (i=kk+1; i<=n; i++)
		    v[(i-1)*n+kk-1]=e[i-1];
	    }
	}
    }
    mm=n;
    if (m+1<n) mm=m+1;
    if (k<n) s[k]=a[k*n+k];
    if (m<mm) s[mm-1]=0.0;
    if (l+1<mm) e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) nn=n;
    if (nn>=k+1) {
	for (j=k+1; j<=nn; j++) {
	    for (i=1; i<=m; i++)
		u[(i-1)*m+j-1]=0.0;
            u[(j-1)*m+j-1]=1.0;
	}
    }
    if (k>=1) {
	for (ll=1; ll<=k; ll++) {
	    kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0) {
		if (nn>=kk+1)
		    for (j=kk+1; j<=nn; j++) {
			d=0.0;
			for (i=kk; i<=m; i++) {
			    ix=(i-1)*m+kk-1;
			    iy=(i-1)*m+j-1;
			    d=d+u[ix]*u[iy]/u[iz];
                        }
			d=-d;
			for (i=kk; i<=m; i++) {
			    ix=(i-1)*m+j-1;
			    iy=(i-1)*m+kk-1;
			    u[ix]=u[ix]+d*u[iy];
                        }
                    }
		for (i=kk; i<=m; i++) {
		    ix=(i-1)*m+kk-1; u[ix]=-u[ix];
		}
		u[iz]=1.0+u[iz];
		if (kk-1>=1)
                    for (i=1; i<=kk-1; i++)
			u[(i-1)*m+kk-1]=0.0;
	    }
            else {
		for (i=1; i<=m; i++)
		    u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
	    }
	}
    }
    for (ll=1; ll<=n; ll++) {
	kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0)) {
	    for (j=kk+1; j<=n; j++) {
		d=0.0;
                for (i=kk+1; i<=n; i++) {
		    ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
		}
                d=-d;
                for (i=kk+1; i<=n; i++) {
		    ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
		}
	    }
	}
        for (i=1; i<=n; i++)
	    v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
    }
    for (i=1; i<=m; i++)
	for (j=1; j<=n; j++)
	    a[(i-1)*n+j-1]=0.0;
    m1=mm; it=60;
    while (1==1) {
	if (mm==0) {
	    ppp(a,e,s,v,m,n);
	    return(1);
	}
        if (it==0) { 
	    ppp(a,e,s,v,m,n);
	    return(-1);
	}
        kk=mm-1;
	while ((kk!=0)&&(fabsf(e[kk-1])!=0.0)) {
	    d=fabsf(s[kk-1])+fabsf(s[kk]);
            dd=fabsf(e[kk-1]);
            if (dd>eps*d) kk=kk-1;
            else e[kk-1]=0.0;
	}
        if (kk==mm-1) {
	    kk=kk+1;
            if (s[kk-1]<0.0) {
		s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++) {
		    ix=(i-1)*n+kk-1; v[ix]=-v[ix];}
	    }
            while ((kk!=m1)&&(s[kk-1]<s[kk])) {
		d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
		    for (i=1; i<=n; i++) {
			ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
			d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                    }
                if (kk<m)
		    for (i=1; i<=m; i++) {
			ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
			d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                    }
                kk=kk+1;
	    }
            it=60;
            mm=mm-1;
	}
        else {
	    ks=mm;
            while ((ks>kk)&&(fabsf(s[ks-1])!=0.0)) {
		d=0.0;
                if (ks!=mm) d=d+fabsf(e[ks-1]);
                if (ks!=kk+1) d=d+fabsf(e[ks-2]);
                dd=fabsf(s[ks-1]);
                if (dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
	    }
            if (ks==kk) {
		kk=kk+1;
                d=fabsf(s[mm-1]);
                t=fabsf(s[mm-2]);
                if (t>d) d=t;
                t=fabsf(e[mm-2]);
                if (t>d) d=t;
                t=fabsf(s[kk-1]);
                if (t>d) d=t;
                t=fabsf(e[kk-1]);
                if (t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0)) {
		    shh=sqrtf(b*b+c);
                    if (b<0.0) shh=-shh;
                    shh=c/(b+shh);
		}
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++) {
		    sss(fg,cs);
                    if (i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
			for (j=1; j<=n; j++) {
			    ix=(j-1)*n+i-1;
			    iy=(j-1)*n+i;
			    d=cs[0]*v[ix]+cs[1]*v[iy];
			    v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
			    v[ix]=d;
                        }
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
			if ((cs[0]!=1.0)||(cs[1]!=0.0))
			    for (j=1; j<=m; j++) {
				ix=(j-1)*m+i-1;
				iy=(j-1)*m+i;
				d=cs[0]*u[ix]+cs[1]*u[iy];
				u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
				u[ix]=d;
			    }
		}
                e[mm-2]=fg[0];
                it=it-1;
	    }
            else {
		if (ks==mm) {
		    kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++) {
			i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk) {
			    fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
			}
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
			    for (j=1; j<=n; j++) {
				ix=(j-1)*n+i-1;
				iy=(j-1)*n+mm-1;
				d=cs[0]*v[ix]+cs[1]*v[iy];
				v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
				v[ix]=d;
                            }
		    }
		}
                else {
		    kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++) {
			fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
			    for (j=1; j<=m; j++) {
				ix=(j-1)*m+i-1;
				iy=(j-1)*m+kk-2;
				d=cs[0]*u[ix]+cs[1]*u[iy];
				u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
				u[ix]=d;
                            }
		    }
		}
	    }
	}
    }
    return(1);
}

/* 	$Id$	 */
