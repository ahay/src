/* 2D Digital Wavelet Transoform Denoising */
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

int main(int argc, char *argv[])
{
    int i,j, n1,n2,n3,n, nthr; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    float *pp, *qq, *data1, *data2, *adata;
    sf_file in, out; /*, outwav; */
    float pclip, t;
    char *type;
    bool unit=false, inv=true;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/

    if (!sf_getfloat("pclip",&pclip)) pclip=99.;
    /* data clip percentile (default is 99)*/
    if (pclip <=0. || pclip > 100.)
	sf_error("pclip=%g should be > 0 and <= 100",pclip);
    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */
 
    if(n1>n2)n=n1;else n=n2; 
    pp = sf_floatalloc(n);      /*allocate memory*/
    qq = sf_floatalloc(n);	/*allocate memory*/

    data1 = sf_floatalloc(n1*n2);
    data2 = sf_floatalloc(n1*n2);
    adata = sf_floatalloc(n1*n2);

    sf_wavelet_init(n1,inv,unit,type[0]);  /* unit=false inv=true ; the length of the first axis is n1 */

    for(i=0;i<n3;i++)  {
	sf_floatread(data1,n1*n2,in);
	for (j=0;j<n2;j++){
	    memcpy(pp,data1+j*n1,n1*sizeof(float));
	    sf_wavelet_lop(0,false,n1,n1,pp,qq);	/*qq -> output*/
	    memcpy(data1+j*n1,qq,n1*sizeof(float));
	}

	for(i=0;i<n2;i++)
	    for(j=0;j<n1;j++)
		data2[j*n2+i]=data1[i*n1+j];

    	sf_wavelet_init(n2,inv,unit,type[0]);  /* unit=false inv=true ; the length of the first axis is n2 */

	for (j=0;j<n1;j++){
	    memcpy(pp,data2+j*n2,n2*sizeof(float));
	    sf_wavelet_lop(0,false,n2,n2,pp,qq);	/*qq -> output*/
	    memcpy(data2+j*n2,qq,n2*sizeof(float));
	}
	/***********************************************************/
	/* percentile thresholding */
	for(i=0;i<n1*n2;i++)
	    adata[i]=fabs(data2[i]);	
	
   	nthr = 0.5+n1*n2*(1.-0.01*pclip);
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2) nthr=n1*n2-1;
	t=sf_quantile(nthr,n1*n2,adata);
	
	for(i=0;i<n1*n2;i++)
	    if(fabs(data2[i])<t) data2[i]=0;
	/***********************************************************/
	for (j=0;j<n1;j++){
	    memcpy(qq,data2+j*n2,n2*sizeof(float));
	    sf_wavelet_lop(1,false,n2,n2,pp,qq);    /*pp -> output*/
	    memcpy(data2+j*n2,pp,n2*sizeof(float));
	}
	for(i=0;i<n1;i++)
	    for(j=0;j<n2;j++)
		data1[j*n1+i]=data2[i*n2+j];
	
    	sf_wavelet_init(n1,inv,unit,type[0]);  /* unit=false inv=true ; the length of the first axis is n1 */

	for (j=0;j<n2;j++){
	    memcpy(qq,data1+j*n1,n1*sizeof(float));
	    sf_wavelet_lop(1,false,n1,n1,pp,qq);	/*pp -> output*/
	    memcpy(data1+j*n1,pp,n1*sizeof(float));
	}
	sf_floatwrite(data1,n1*n2,out);
    }
    sf_wavelet_close();
    exit(0);
}

