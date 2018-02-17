/* Generalized thresholding operator
 */
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <complex.h>

int main(int argc, char* argv[])
{
    char *mode;
    int i, n, n1;
    float *dat=NULL, *adat=NULL, pclip, p, thr, a, b;
    sf_complex *cdat=NULL;
    sf_file in=NULL, out=NULL;

    /* input and output variables */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getfloat("pclip",&pclip)) 	pclip=99.;/* percentage to clip */
    if ( !(mode=sf_getstring("mode")) ) mode = "exp";
    /* thresholding mode='hard', 'soft','pthresh' or 'exp';
       'hard', hard thresholding;	'soft', soft thresholding; 
       'pthresh', generalized quasi-p; 'exp', exponential shrinkage */
    if (!sf_getfloat("p",&p)) 		p=0.35; /* norm=p, where 0<p<=1 */;
    if (strcmp(mode,"soft") == 0) 	p=1;
    else if (strcmp(mode,"hard") == 0) 	p=0;

    n = sf_filesize(in);
    adat= sf_floatalloc(n);

    if (SF_FLOAT == sf_gettype(in)) {
	dat = sf_floatalloc(n);
	sf_floatread(dat,n,in);
	for (i=0; i < n; i++) {
	    adat[i] = fabsf(dat[i]);
	}
    } else if (SF_COMPLEX == sf_gettype(in)) {
	cdat = sf_complexalloc(n);
	sf_complexread(cdat,n,in);
	for (i=0; i < n; i++) {
	    adat[i] = cabsf(cdat[i]);
	}
    } else {
	sf_error("Need float or complex input");
    }

    n1 = 0.5+n*(1.-0.01*pclip);
    if (n1 < 0) n1=0;
    if (n1 >= n) n1=n-1;
    thr=sf_quantile(n1,n,adat);

    if (thr>0.){// do thresholding if thr>0
	if (NULL != dat) {// floating numbers
	    if (strcmp(mode,"hard")==0){
		for (i=0; i < n; i++)  {
		    a=fabsf(dat[i]);
		    dat[i]*=(a>thr)?1.:0.;
		}
	    }else if(strcmp(mode,"soft")==0){
		for (i=0; i < n; i++)  {
		    a=fabsf(dat[i]);
		    b=(a==0.)?1.:0.;
		    a=1.0-thr/(a+b);
		    dat[i]*=(a>0.)?a:0.;
		}
	    }else if (strcmp(mode,"pthresh") == 0){
		for (i=0; i < n; i++)  {
		    a=fabsf(dat[i]);
		    b=(a==0.)?1.:0.;
		    a=1.0-powf((a+b)/thr, p-2.0);
		    dat[i]*=(a>0.)?a:0.;
		}
	    }else if (strcmp(mode,"exp") == 0){
		for (i=0; i < n; i++)  {
		    a=fabsf(dat[i]);
		    b=(a==0.)?1.:0.;
		    a=expf(-powf((a+b)/thr, p-2.0));
		    dat[i]*=(a>0.)?a:0.;
		}
	    }
	    sf_floatwrite(dat,n,out);
	} else {// complex numbers
	    if (strcmp(mode,"hard")==0){
		for (i=0; i < n; i++)  {
		    a=cabsf(cdat[i]);
		    cdat[i]*=(a>thr)?1.:0.;
		}
	    }else if(strcmp(mode,"soft")==0){
		for (i=0; i < n; i++)  {
		    a=cabsf(cdat[i]);
		    b=(a==0.)?1.:0.;
		    a=1.0-thr/(a+b);
		    cdat[i]*=(a>0.)?a:0.;
		}
	    }else if (strcmp(mode,"pthresh") == 0){
		for (i=0; i < n; i++)  {
		    a=cabsf(cdat[i]);
		    b=(a==0.)?1.:0.;
		    a=1.0-powf((a+b)/thr, p-2.0);
		    cdat[i]*=(a>0.)?a:0.;
		}
	    }else if (strcmp(mode,"exp") == 0){
		for (i=0; i < n; i++)  {
		    a=cabsf(cdat[i]);
		    b=(a==0.)?1.:0.;
		    a=expf(-powf((a+b)/thr, p-2.0));
		    cdat[i]*=(a>0.)?a:0.;
		}
	    }
	    sf_complexwrite(cdat,n,out);
	}
    }

    exit(0);
}
