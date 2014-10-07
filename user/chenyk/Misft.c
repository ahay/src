/* Iterative seislet frame thresholding  */
/*
  Copyright (C) 2014 University of Texas at Austin
   
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
#include <rsfpwd.h>
void thresholding(float *xx, int n, float thr, float thrtp, char* mode);

int main(int argc, char *argv[])
{
    int i1,i2,ik,n1,n2,nk,iter,niter, nthr, order;
    float pthr, thrtp, ps, thr, eps, m, sum;
    float *din, *dout, *dseis, *dips, **dip, *mask, *diff, *tmp, *misfit;
    char *type, *mode;
    bool verb, ifdecay;
    sf_file in, out, Fs, Fmis;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    Fs=sf_input("dips");
    Fmis=sf_output("misfit");

    if (!sf_histint(in,"n1",&n1)) 	sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) 	sf_error("No n2= in input");
    if (!sf_histint(Fs,"n3",&nk)) 	nk=1;
    if(!sf_getfloat("eps",&eps)) 	eps=0.01;

    if(!sf_getint("order",&order)) 	order=1;

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    if (!sf_getfloat("ps",&ps)) 	ps=3;
    /* starting percentage */

    if (!sf_getfloat("pthr",&pthr)) pthr=18;
    /* percentile thresholding */

    if(!sf_getbool("verb",&verb))    	verb=false;
    /* If output verbosity information? */

    if(!sf_getbool("ifdecay",&ifdecay))    	ifdecay=true;
    /* If decay the threshold */

    if (!sf_getint("niter",&niter)) 	niter=20;
    /* total number iterations */

    if ( !(mode=sf_getstring("mode")) ) mode = "e";
    /* thresholding type */

    if (!sf_getfloat("thrtype",&thrtp)) thrtp=0.5;
    if (mode[0]=='s') 	thrtp=1;
    else if (mode[0]=='h') 	thrtp=0;

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putint(out,"n3",nk);
    sf_putint(Fmis,"n1",niter);
    sf_putint(Fmis,"n2",1);
    sf_putfloat(Fmis,"d1",1);
    sf_putfloat(Fmis,"o1",1);

    din=sf_floatalloc(n1*n2);
    diff=sf_floatalloc(n1*n2);
    dout=sf_floatalloc(n1*n2*nk);
    dips=sf_floatalloc(n1*n2*nk);
    misfit=sf_floatalloc(niter);
    dip=sf_floatalloc2(n1, n2);	
    dseis=sf_floatalloc(n1*n2);
    tmp=sf_floatalloc(n1*n2);
    mask=sf_floatalloc(n1*n2);		
	
    memset(dout, 0, n1*n2*nk*sizeof(float));
    memset(dseis, 0, n1*n2*sizeof(float));
    memset(tmp, 0, n1*n2*sizeof(float));
    memset(diff, 0, n1*n2*sizeof(float));
    memset(dip[0], 0, n1*n2*sizeof(float));

    sf_file Fm;
    if (NULL != sf_getstring("mask")){
    	Fm=sf_input("mask");  
	sf_floatread(mask, n1*n2, Fm);
    }else{
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
		mask[i1+i2*n1]=1.0;
    }	

    sf_floatread(din, n1*n2, in);
    sf_floatread(dips, n1*n2*nk, Fs);

    seislet_init(n1, n2, true, false, eps, order, type[0]);
    seislet_set(dip);

    for(iter=0; iter<niter; iter++)
    {
	memset(diff, 0, n1*n2*sizeof(float));
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
	{	
		for(ik=0; ik<nk; ik++) diff[i1+n1*i2]+=dout[i1+n1*i2+n1*n2*ik];
		m=(mask[i1+i2*n1])?1.:0; 
		diff[i1+n1*i2]=din[i1+n1*i2]-m*diff[i1+n1*i2];
		sum+=fabs(diff[i1+n1*i2]);
	}
	sum=sum/n1/n2;
	sf_warning("Data misfit is %g",sum);
	misfit[iter]=sum;	
		
	for(ik=0; ik<nk; ik++)
	{
		memcpy(dip[0], &dips[ik*n1*n2], n1*n2*sizeof(float));

		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++) 
		{	
			dout[i1+n1*i2+n1*n2*ik]+=diff[i1+n1*i2];
		}
		seislet_lop(true, false, n1*n2, n1*n2, dseis, &dout[ik*n1*n2]);

		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++) 
		{
			if (i2>0.01*pthr*n2) dseis[i1+i2*n1]=0;
			tmp[i1+n1*i2]=fabsf(dseis[i1+n1*i2]);
		}
		nthr = 0.5+n1*n2*(1.-0.01*ps);  
		if (nthr < 0) nthr=0;
		if (nthr >= n1*n2) nthr=n1*n2-1;
		thr=sf_quantile(nthr, n1*n2, tmp);
		if(ifdecay) thr*=(float)(niter-iter)/niter;
		thresholding(dseis, n1*n2, thr, thrtp, mode);

		seislet_lop(false, false, n1*n2, n1*n2, dseis, &dout[ik*n1*n2]);
	}

	if (verb) sf_warning("iteration %d;",iter+1);	
    }
    sf_floatwrite(misfit,niter,Fmis);
    sf_floatwrite(dout, n1*n2*nk, out);

    free(din);
    free(dout);
    free(diff);
    free(dseis);
    free(dips);
    free(*dip); 
    free(dip);
    free(tmp);
    free(mask);
    exit(0);
}

void thresholding(float *xx, int nx, float thr, float thrtp, char* mode)
{
    float p1,p2;
    int i;

	for(i=0;i<nx;i++){
	  p1=fabsf(xx[i]);
	  p2=(p1==0.)?1.0:0.0;
	  if (mode[0]=='h') xx[i]=(xx[i])*(p1>thr?1.:0.);
	  if(mode[0]=='s') p1=1.0-thr/(p1+p2);
	  if(mode[0]=='g') p1=1.0-thr/(p1+p2);
	  if(mode[0]=='e') p1=expf(-powf((p1+p2)/thr, thrtp-2.0));
	  xx[i]=(xx[i])*(p1>0.0?p1:0.0);	
	}
}
