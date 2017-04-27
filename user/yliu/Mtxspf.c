/* Streaming prediction filter in t-x domain. */
/*
  Copyright (C) 2017 Jilin University
  
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

int main(int argc, char* argv[])
{
    int i1,i2,it,ix,n1,n2,n12,dim,na,i,nst;
    int a[SF_MAX_DIM],n[SF_MAX_DIM];
    float dd,da,dn,rn,eps,eps1,eps2;
    float *d,*aa,*r,*st,*t;
    sf_file in,out;
    
    sf_init(argc,argv);
    
    in = sf_input("in");
    out = sf_output("out");
    
    dim = sf_filedims(in,n);
    if (2 < dim) dim = 2;
    
    if (!sf_getints("a",a,dim)) sf_error("Need a=");
    
    if (dim < 2) sf_error("Need at least two dimension");
    
    a[1]=a[1]*2;
    n12 = 1;
    na = 1;
    
    for (i=0; i < dim; i++) {
	n12 *= n[i];
	na *= a[i];
    }
    
    n1=n[0];
    n2=n[1];
    nst=na*n1;
    
    if (!sf_getfloat("eps1",&eps1)) sf_error("Need eps1=");
    /* regularization in t direction */
    eps1*=eps1;
    if (!sf_getfloat("eps2",&eps2)) sf_error("Need eps2=");
    /* regularization in x direction */
    eps2*=eps2;
    eps=eps1+eps2;
    
    d = sf_floatalloc(n12);
    t = sf_floatalloc(n12);
    r = sf_floatalloc(n12);
    aa = sf_floatalloc(na);
    st = sf_floatalloc(nst);
    
    for (i=0;i<na;i++){
	aa[i]=0.0f;
    }
    for (i=0;i<nst;i++){
	st[i]=0.0f;
    }
    
    sf_floatread(d,n12,in);
    
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    dd = 0.0f;
	    da = 0.0f;
	    i=0;
	    for (ix=-a[1]/2; ix < a[1]/2+1; ix++) {
		for (it=-a[0]/2; it < (a[0]+1)/2; it++) {
		    if(ix!=0){
			if(i2+ix<0 ||
			   i2+ix>=n2 ||
			   i1+it<0 ||
			   i1+it>=n1){
			    i++;
			    continue;
			} else{
			    dd += d[(i2+ix)*n1+i1+it]*
				d[(i2+ix)*n1+i1+it];
			    
			    da += d[(i2+ix)*n1+i1+it]*
				(eps1*aa[i]+eps2*st[i1*na+i])/eps;
			    i++;
			}
		    }
		}
	    }
	    
	    dn=d[i2*n1+i1];
	    rn = (dn+da)/(eps+dd);
	    r[i2*n1+i1] = eps*rn;

	    i=0;
	    for (ix=-a[1]/2; ix < a[1]/2+1; ix++) {
		for (it=-a[0]/2; it < (a[0]+1)/2; it++) {
		    if(ix!=0){
			if(i2+ix<0 ||
			   i2+ix>=n2 ||
			   i1+it<0 ||
			   i1+it>=n1){
			    i++;
			    continue;
			} else{

			    aa[i] = (eps1*aa[i]+eps2*st[i1*na+i])/
				eps-rn*d[(i2+ix)*n1+i1+it];
			    i++;
			}
		    }
		}
	    }
	    for (i=0; i < na; i++) {
		st[i1*na+i]=aa[i];
	    }
	}
    }


    sf_floatwrite(r,n12,out);

    exit(0);
}

/* 	$Id$	 */

