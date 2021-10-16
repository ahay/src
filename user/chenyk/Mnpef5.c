/* Estimate adaptive nonstationary PEF on aliased traces (5D). */
/*
  Copyright (C) 2020 University of Texas at Austin
  
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

#include <math.h>
#include "npef_mask.h"

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], m[SF_MAX_DIM], rect[SF_MAX_DIM], a[SF_MAX_DIM];
    int ndim, mdim, nd, ns, n12, i, j, niter;
    int i1, i2, i3, i4, i5, j1, j2, j3, j4, j5, jump, jumps[SF_MAX_DIM], i6, n6;
    float *d, *f, *g, mean, *ff, *outm;
    char key[6];
    sf_file flt, mat, pre, maskin, maskout;
    bool verb, *mask;

    sf_init(argc,argv);

    mat = sf_input("in");
    flt = sf_output("out");

    if (NULL != sf_getstring("pred")) {
	pre = sf_output("pred");
    } else {
	pre = NULL;
    }

    if (!sf_getint("jump",&jump)) jump=2;
    /* Jump parameter */

    mdim = sf_filedims(mat,m);
    
    if (!sf_getint("dim",&ndim)) ndim=mdim; /* number of dimensions */

    n6 = sf_leftsize(mat,ndim);

    if (!sf_getints("a",a,ndim)) sf_error("Need a=");

    if (!sf_getints("j",jump,ndim)) {jumps[0]=jump;jumps[1]=jump;jumps[2]=jump;jumps[3]=jump;jumps[4]=jump;}


    if (ndim < 5) {
	a[4] = 1;
	n[4] = 1;
	m[4] = 1;
    }
    
    if (ndim < 4) {
	a[3] = 1;
	n[3] = 1;
	m[3] = 1;
    }
    
    if (ndim < 3) {
	a[2] = 1;
	n[2] = 1;
	m[2] = 1;
    }
    if (ndim < 2) {
	a[1] = 1;
	n[1] = 1;
	m[1] = 1;
    }

    for (j=0; j < ndim; j++) {
	n[j] = m[j];
    }
    n[ndim] =1;
    for (j=0; j < ndim; j++) {
	n[ndim] *= a[j];
    }
/*    ndim = mdim+1; */

    nd = 1;
    for (j=0; j < ndim; j++) {
	nd *= m[j];
	snprintf(key,6,"rect%d",j+1);
	if (!sf_getint(key,rect+j)) rect[j]=1;
    }

    ns = n[ndim];
    n12 = nd*ns;

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (NULL != sf_getstring("maskin")) {
	/* optional input mask file */
	mask = sf_boolalloc(nd);
	maskin = sf_input("maskin");
    } else {
	mask = NULL;
	maskin = NULL;
    }
    if (NULL != sf_getstring("maskout")) {
	/* optional output mask file */
	maskout = sf_output("maskout");
	outm = sf_floatalloc(nd);
    } else {
	outm = NULL;
	maskout = NULL;
    }	


	sf_warning("jump=%d",jump);
	sf_warning("jumps=%d,%d,%d,%d,%d",jumps[0],jumps[1],jumps[2],jumps[3],jumps[4]);
	sf_warning("ndim=%d,mdim=%d",ndim,mdim);
	sf_warning("n[5]=%d",n[5]);
	sf_warning("nd=%d,ns=%d,n12=%d",nd,ns,n12);
	
    d = sf_floatalloc(n12);
    f = sf_floatalloc(n12);
    ff = sf_floatalloc(n12);
    g = sf_floatalloc(nd);

    sf_multidivn_init(ns, ndim, nd, m, rect, d, NULL, verb); 

    if (1==m[2] && 1==m[1]) {
	sf_shiftdim(mat,flt,1);
	sf_putint(flt,"n1",a[0]);
    } else if (1==m[2] && 1!=m[1]) {
	sf_shiftdim2(mat,flt,1);
	sf_putint(flt,"n1",a[0]);
	sf_putint(flt,"n2",a[1]);
    } else {
	sf_putint(flt,"n1",a[0]);
	sf_putint(flt,"n2",a[1]);
	sf_putint(flt,"n3",a[2]);
	sf_putint(flt,"n4",a[3]);
	sf_putint(flt,"n5",a[4]);
	
	sf_putint(flt,"n6",m[0]);
	sf_putint(flt,"n7",m[1]);
	sf_putint(flt,"n8",m[2]);
	sf_putint(flt,"n9",m[3]*m[4]);
// 	sf_putint(flt,"n10",m[4]);
// 	sf_putint(flt,"n11",m[5]); /*Madagascar does not support*/

    }

    for (i6=0; i6 < n6; i6++) {
    	if (NULL != maskin) {
	    sf_floatread(g,nd,maskin);
	    mask4apef (a, jumps, m, g, mask);
	    if (NULL != maskout) {
		for (i3=0; i3 < nd; i3++) {
		    if (!mask[i3]) {
			outm[i3] = 0.;
		    } else {
			outm[i3] = 1.;
		    }
		}
		sf_floatwrite (outm,nd,maskout);
	    }
	}
	int sum1=0, sum2=0;
	for(i1=0;i1<nd;i1++)
	{sum1=sum1+mask[i1];sum2=sum2+g[i1];}
	/*sf_warning("sum1=%g, sum2=%g",sum1,sum2);*///mask[i1]=g[i1];
	sf_warning("Available samples for NPEF=%d/%d, percentage=%g %",sum1,nd,(float)sum1/nd*100);
	
	sf_floatread(g,nd,mat);
	
	for (i5=0; i5 < a[4]; i5++) {
	for (i4=0; i4 < a[3]; i4++) {
	for (i3=0; i3 < a[2]; i3++) {
	    for (i2=0; i2 < a[1]; i2++) {
		for (i1=-a[0]/2; i1 < (a[0]+1)/2; i1++) {
		    for (j5=0; j5 < m[4]; j5++) {
		    for (j4=0; j4 < m[3]; j4++) {		
		    for (j3=0; j3 < m[2]; j3++) {
			for (j2=0; j2 < m[1]; j2++) {
			    for (j1=0; j1 < m[0]; j1++) {
			    
				if (0==i5&& 0==i4 && 0==i2 && 0==i3 && 0 >=i1) {
				    d[
				      i5*a[3]*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
				      i4*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
				      i3*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
				      i2*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
				      (i1+a[0]/2)*m[4]*m[3]*m[2]*m[1]*m[0]+
				      j5*m[3]*m[2]*m[1]*m[0]+j4*m[2]*m[1]*m[0]+j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
				      
				} else {
				    if (
					(j5+i5*jumps[4])<0 || (j5+i5*jumps[4])>=m[4] ||
					(j4+i4*jumps[3])<0 || (j4+i4*jumps[3])>=m[3] ||
					(j1+i1*jumps[0])<0 || (j1+i1*jumps[0])>=m[0] ||
					(j2+i2*jumps[1])<0 || (j2+i2*jumps[1])>=m[1] ||
					(j3+i3*jumps[2])<0 || (j3+i3*jumps[2])>=m[2]) {
					d[
					  i5*a[3]*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i4*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+ 
					  i3*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[4]*m[3]*m[2]*m[1]*m[0]+
					  j5*m[3]*m[2]*m[1]*m[0]+j4*m[2]*m[1]*m[0]+j3*m[1]*m[0]+j2*m[0]+j1] = 0.; 
				    
				    } else {
					d[i5*a[3]*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i4*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i3*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[4]*m[3]*m[2]*m[1]*m[0]+
					  j5*m[3]*m[2]*m[1]*m[0]+j4*m[2]*m[1]*m[0]+j3*m[1]*m[0]+j2*m[0]+j1] = 
					    g[(j5+i5*jumps[4])*m[3]*m[2]*m[1]*m[0]+
					      (j4+i4*jumps[3])*m[2]*m[1]*m[0]+
					      (j3+i3*jumps[2])*m[1]*m[0]+
					      (j2+i2*jumps[1])*m[0]+
					      (j1+i1*jumps[0])];      
				    }
				}
				
				
				
				
				
			    }
			}
		    }
		}
	    }
	}
	}
	}
	}
	}
	
	mean = 0.;
	for(i=0; i < n12; i++) {
	    mean += d[i]*d[i];
	}
	if (mean == 0.) {
	    sf_floatwrite(d,n12,flt);
	    continue;
	}
	
	mean = sqrtf (mean/n12);
	
	for(i=0; i < n12; i++) {
	    d[i] /= mean;
	}
	for(i=0; i < nd; i++) {
	    g[i] /= mean;
	}


    	if (NULL != maskin) {
	    for (j5=0; j5 < m[4]; j5++) {
	    for (j4=0; j4 < m[3]; j4++) {  	
	    for (j3=0; j3 < m[2]; j3++) {
		for (j2=0; j2 < m[1]; j2++) {
		    for (j1=0; j1 < m[0]; j1++) {
			if (!mask[j5*m[3]*m[2]*m[1]*m[0]+j4*m[2]*m[1]*m[0]+j3*m[1]*m[0]+j2*m[0]+j1]) {
			    for (i5=0; i5 < a[4]; i5++) {
			    for (i4=0; i4 < a[3]; i4++) {			
			    for (i3=0; i3 < a[2]; i3++) {
				for (i2=0; i2 < a[1]; i2++) {
				    for (i1=0; i1 < a[0]; i1++) {
					d[i5*a[3]*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i4*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i3*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i1*m[4]*m[3]*m[2]*m[1]*m[0]+
					  j5*m[3]*m[2]*m[1]*m[0]+j4*m[2]*m[1]*m[0]+j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
				}
				}
			    	}
				}
				}
			    g[j5*m[3]*m[2]*m[1]*m[0]+j4*m[2]*m[1]*m[0]+j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
			}
		    }
		}
	    } 
	}
	}
	}
	
	sf_multidivn (g,f,niter);

	for (i5=0; i5 < a[4]; i5++) {
	for (i4=0; i4 < a[3]; i4++) {	
	for (i3=0; i3 < a[2]; i3++) {
	    for (i2=0; i2 < a[1]; i2++) {
		for (i1=0; i1 < a[0]; i1++) {
		    for (j5=0; j5 < m[4]; j5++) {
			for (j4=0; j4 < m[3]; j4++) {		
		    for (j3=0; j3 < m[2]; j3++) {
			for (j2=0; j2 < m[1]; j2++) {
			    for (j1=0; j1 < m[0]; j1++) {
				if (0==i5 && 0==i4 && 0==i2 && 0==i3 && a[0]/2==i1) {
				    ff[j5*m[3]*m[2]*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j4*m[2]*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j3*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j2*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j1*a[4]*a[3]*a[2]*a[1]*a[0]+
				       i5*a[3]*a[2]*a[1]*a[0]+
				       i4*a[2]*a[1]*a[0]+
				       i3*a[1]*a[0]+
				       i2*a[0]+
				       i1] = 1.0;
				} else if (0==i5 && 0==i4&&0==i2 && 0==i3 && a[0]/2>i1) {
				    ff[j5*m[3]*m[2]*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j4*m[2]*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j3*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j2*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j1*a[4]*a[3]*a[2]*a[1]*a[0]+
				       i5*a[3]*a[2]*a[1]*a[0]+
				       i4*a[2]*a[1]*a[0]+
				       i3*a[1]*a[0]+
				       i2*a[0]+
				       i1] = 0.0;
				} else {
				    ff[j5*m[3]*m[2]*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j4*m[2]*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j3*m[1]*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j2*m[0]*a[4]*a[3]*a[2]*a[1]*a[0]+
				       j1*a[4]*a[3]*a[2]*a[1]*a[0]+
				       i5*a[3]*a[2]*a[1]*a[0]+
				       i4*a[2]*a[1]*a[0]+
				       i3*a[1]*a[0]+
				       i2*a[0]+
				       i1] = -1.*
					f[i5*a[3]*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i4*a[2]*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+ 
					  i3*a[1]*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[4]*m[3]*m[2]*m[1]*m[0]+
					  i1*m[4]*m[3]*m[2]*m[1]*m[0]+
					  j5*m[3]*m[2]*m[1]*m[0]+
					  j4*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+
					  j2*m[0]+
					  j1];
				}
			    }
			}
		    }
		}
	    }
	}
	}
	}
	}
	}

	sf_floatwrite(ff,n12,flt);	
	
	if (pre) {
	    for(i=0; i < n12; i++) {
		d[i] *= mean;
	    }
	    sf_weight2_lop(false,false,n12,nd,f,g);
	    sf_floatwrite(g,nd,pre);
	}
    }
  
    exit(0);
}


