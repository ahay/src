/* Accelerated local orthogonalization */
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
#include "memcpy.h"
#include "win.h"
//#include "divnn.h"

int main(int argc, char* argv[])
{
    bool verb, ifsm, ifnm, ifns;
    int *sft[SF_MAX_DIM];	/* storing non-stationary shifting size */
    int i, id, dim, n[SF_MAX_DIM], dim1, nd, i1, i2, i3, n1, n2,n3, b, nrad;
    float *noi, *sig, *rat, *rat1, *noi2, *sig2, remove, alpha, beta, gamma;
    char key[8];
    int nbox[SF_MAX_DIM],box[SF_MAX_DIM];
    float *rct[SF_MAX_DIM]; /* storing non-stationary smoothing radii */
    sf_file fnoi, fsig, fnoi2, fsig2, fw;
	sf_file rect[SF_MAX_DIM], shift[SF_MAX_DIM]; 

    
    sf_init(argc,argv);
    fnoi = sf_input("in");
    fsig = sf_input("sig"); /* input signal */

    if (SF_FLOAT != sf_gettype(fnoi) ||
	SF_FLOAT != sf_gettype(fsig)) sf_error("Need float input");

    fw = sf_output("w");
    fnoi2 = sf_output("out");
    fsig2 = sf_output("sig2"); /* output signal */

    dim = sf_filedims (fnoi,n);
    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
    }

    if (!sf_getbool("ifns",&ifns)) ifns=true;
    /* if apply the non-stationary smoothing radius*/ 
    
    if(ifns)/*if apply the non-stationary smoothing radius*/
    {
	/*Calculate dim1*/
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (NULL != sf_getstring(key)) {
	    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
	    rect[i] = sf_input(key);
	    if (SF_FLOAT != sf_gettype(rect[i])) sf_error("Need float %s",key);
	    dim1 = i;
	    snprintf(key,8,"shift%d",i+1);
	    if (NULL != sf_getstring(key)) {
		/*( shift# shifting of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
		shift[i] = sf_input(key);
		if (SF_INT != sf_gettype(shift[i])) sf_error("Need int %s",key);
	    } else {
		shift[i] = NULL;
	    }
	} else {
	    rect[i] = NULL;
	    shift[i] = NULL;
	}
    }
    }/*end if*/
    
	/*Calculate n1*/
    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }
    
    noi = sf_floatalloc(nd);
    sig = sf_floatalloc(nd);
    rat = sf_floatalloc(nd);    
    rat1 = sf_floatalloc(nd);    
    
    noi2 = sf_floatalloc(nd);
    sig2 = sf_floatalloc(nd);

	if(ifns)/*if apply the non-stationary smoothing radius*/
	{
	/*reading the non-stationary smoothing radii*/
    for (i=0; i <= dim1; i++) {
	box[i] = 1;
	if (NULL != rect[i]) {
	    rct[i] = sf_floatalloc (n1);
	    sft[i] = sf_intalloc (n1);

	    sf_floatread(rct[i],n1,rect[i]);
	    sf_fileclose(rect[i]);

	    if (NULL != shift[i]) {
		sf_intread(sft[i],n1,shift[i]);
		sf_fileclose(shift[i]);
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    sft[i][i1] = 0;
		}
	    }

		
	    for (i1=0; i1 < n1; i1++) {
		b = ceilf(rct[i][i1])+SF_ABS(sft[i][i1]);
		if (b > box[i]) box[i] = b;
	    }	    
	} else {
	    rct[i] = NULL;
	    sft[i] = NULL;
	}
    }
    }/*end if*/

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getbool("ifsm",&ifsm)) ifsm=true;
    /* if smooth the weight*/

    if (!sf_getbool("ifnm",&ifnm)) ifnm=true;
    /* if normalize the weight*/

    if (!sf_getint("nrad",&nrad)) nrad=1;
    /* NO of smoothing*/
    
    if (!sf_getfloat("alpha",&alpha)) alpha=0.0;
    /* Regularization in t direction */
    alpha *= alpha;
    
    if (!sf_getfloat("beta",&beta)) beta=0.0;
    /* Regularization in x direction */
    beta *= beta;
    
    if (!sf_getfloat("gamma",&gamma)) gamma=0.0;
    /* Regularization in y direction */
    gamma *= gamma;
      
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rad%d",i+1);
	if (!sf_getint(key,nbox+i)) nbox[i]=1;
    }
    
	sf_warning("dim=%d, nd=%d",dim,nd);
	sf_warning("n[0]=%d, n[1]=%d, n[2]=%d",n[0],n[1],n[2]);	
	sf_warning("rad[0]=%d, rad[1]=%d, rad[2]=%d",nbox[0],nbox[1],nbox[2]);	
	sf_warning("nd=%d",nd);
	    
    sf_floatread(noi,nd,fnoi);
    sf_floatread(sig,nd,fsig);

    
    n3=n[2];
    n2=n[1];
    n1=n[0];
    
    int iiw1, iiw2, iiw3, s1, s2, s3;
    int n1win, n2win, n3win, iw1, iw2, iw3, nw1, nw2, nw3, nov1, nov2, nov3, ov1, ov2, ov3, n1pad, n2pad, n3pad; 
    float r1, r2, r3;
	/*windowing parameters*/
    if (!sf_getint("n1win",&n1win)) n1win=n1; /*first window length*/
    if (!sf_getint("n2win",&n2win)) n2win=n2; /*second window length*/	
    if (!sf_getint("n3win",&n3win)) n3win=n3; /*second window length*/	
    
	if (!sf_getfloat("r1",&r1)) r1=0.5;		  /*first overlapping ratio*/
	if (!sf_getfloat("r2",&r2)) r2=0.5;		  /*second overlapping ratio*/
	if (!sf_getfloat("r3",&r3)) r3=0.5;		  /*third overlapping ratio*/
    

    rat[0]=0.0f;
    for(i3=0; i3 < n3; i3++) {
	if(n3==1)
	{
    for(i2=0; i2 < n2; i2++) {
	for(i1=0; i1 < n1; i1++) {
	    if(i2-1<0 || i1-1<0) {  /*meaning that n[0] and n[1] must be larger than 1*/
		rat[i3*n1*n2+i2*n1+i1]=0.;
	    } else {
		rat[i3*n1*n2+i2*n1+i1]=(sig[i3*n1*n2+i2*n1+i1]*noi[i3*n1*n2+i2*n1+i1]+
			     alpha*rat[i3*n1*n2+i2*n1+i1-1]+beta*rat[i3*n1*n2+(i2-1)*n1+i1])/
		    (sig[i3*n1*n2+i2*n1+i1]*sig[i3*n1*n2+i2*n1+i1]+alpha+beta);
	    }
	}/*end i1*/
    }/*end i2*/ 
    }else{
    for(i2=0; i2 < n2; i2++) {
	for(i1=0; i1 < n1; i1++) {
	    if(i2-1<0 || i1-1<0 || i3-1<0) {  /*meaning that n[0] and n[1] must be larger than 1*/
		rat[i3*n1*n2+i2*n1+i1]=0.;
	    } else {
		rat[i3*n1*n2+i2*n1+i1]=(sig[i3*n1*n2+i2*n1+i1]*noi[i3*n1*n2+i2*n1+i1]+
			     alpha*rat[i3*n1*n2+i2*n1+i1-1]+beta*rat[i3*n1*n2+(i2-1)*n1+i1]+gamma*rat[(i3-1)*n1*n2+i2*n1+i1])/
		    (sig[i3*n1*n2+i2*n1+i1]*sig[i3*n1*n2+i2*n1+i1]+alpha+beta+gamma);
	    }
	}/*end i1*/
    }/*end i2*/   
    }/*end if*/
    }/*end i3*/
    
	
    
    if(ifsm) /*if smooth*/
    {
    sf_trianglen_init(dim,nbox, n);
	sf_trianglen_lop(false, false, nd, nd, rat, rat1);
    }

    if(ifns)
    {
    sf_ntrianglen_init(dim,box,n,rct,sft,1);
    sf_ntrianglen_lop (false, false, nd, nd, rat, rat1);
    }
    
     float mmax_s=0, mmax_r=0, mmax_n=0;   
    if(ifnm) /*if normalize*/
    {
/*    for(i1=0;i1<nd;i1++)
    	{
    	if(mmax_s<sig[i1]) mmax_s=fabs(sig[i1]);
    	if(mmax_n<noi[i1]) mmax_n=fabs(noi[i1]);
    	if(mmax_r<rat1[i1])mmax_r=fabs(rat1[i1]);
    	}
    if(mmax_r!=0 && mmax_s!=0)	
    for(i1=0;i1<nd;i1++)
    	rat[i1]=rat1[i1]/mmax_r/mmax_s*mmax_n;*/ /*old one*/
    	

	nov1=(1-r1)*n1win;nov2=(1-r2)*n2win;nov3=(1-r3)*n3win;/*non-overlapping size 1,2,3*/
	ov1=r1*n1win;ov2=r2*n2win;ov3=r3*n3win;		/*overlapping size 1,2,3*/
	
	/*padding */	
	n1pad=n1win;nw1=1;while(n1pad<n1){n1pad=n1pad+nov1;nw1++;	 }   
	n2pad=n2win;nw2=1;while(n2pad<n2){n2pad=n2pad+nov2;nw2++;	 }   	
	n3pad=n3win;nw3=1;while(n3pad<n3){n3pad=n3pad+nov3;nw3++;	 }  
	 	
	float *dsig, *dnoi, *drat, ***dtmp, *dout;
	
    dsig =                  sf_floatalloc  (n1pad*n2pad*n3pad); 
    dnoi =                  sf_floatalloc  (n1pad*n2pad*n3pad); 
    drat =         			sf_floatalloc  (n1pad*n2pad*n3pad); 
    dtmp =                  sf_floatalloc3  (n1win,n2win,n3win); 
    dout =                  sf_floatalloc  (n1pad*n2pad*n3pad); 
    
    mcp(dsig, sig, 0, 0, n1*n2*n3);
    mcp(dnoi, noi, 0, 0, n1*n2*n3);
	mcp(drat, rat1, 0, 0, n1*n2*n3);
    
    for(i3=0;i3<n3;i3++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<n1;i1++)
    {
    dsig[i3*n1pad*n2pad+i2*n1pad+i1]=sig[i3*n1*n2+i2*n1+i1];
    dnoi[i3*n1pad*n2pad+i2*n1pad+i1]=noi[i3*n1*n2+i2*n1+i1];
    drat[i3*n1pad*n2pad+i2*n1pad+i1]=rat[i3*n1*n2+i2*n1+i1];    
    }
    
 	for(iw3=0;iw3<nw3;iw3++)
	for(iw2=0;iw2<nw2;iw2++)
		for(iw1=0;iw1<nw1;iw1++)
			{
			s1=iw1*nov1;s2=iw2*nov2;s3=iw3*nov3;

			mcp3d1d(dtmp,drat,0,0,0,s1,s2,s3,n1win,n2win,n3win,n1pad,n2pad,n3pad);

	mmax_s=0;
	mmax_n=0;
	mmax_r=0;
	for(iiw3=0;iiw3<n3win;iiw3++)
	for(iiw2=0;iiw2<n2win;iiw2++)
		for(iiw1=0;iiw1<n1win;iiw1++)
			{
			i1=(s3+iiw3)*n1pad*n2pad+(s2+iiw2)*n1pad+s1+iiw1;
    	if(mmax_s<dsig[i1]) mmax_s=fabs(dsig[i1]);
    	if(mmax_n<dnoi[i1]) mmax_n=fabs(dnoi[i1]);
    	if(mmax_r<drat[i1]) mmax_r=fabs(drat[i1]);
			}
			
			
	if(mmax_r!=0 && mmax_s!=0)	
 	for(iiw3=0;iiw3<n3win;iiw3++)
	for(iiw2=0;iiw2<n2win;iiw2++)
		for(iiw1=0;iiw1<n1win;iiw1++)
			{
			dtmp[iiw3][iiw2][iiw1]=dtmp[iiw3][iiw2][iiw1]/mmax_r/mmax_s*mmax_n;
			}
			
    		win_weight3d(dtmp,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3);
    		mcp_ad1d3d(dout,dtmp,s1,s2,s3,0,0,0,n1pad,n2pad,n3pad,n1win,n2win,n3win);
			}
    	
    		
    for(i3=0;i3<n3;i3++)
    for(i2=0;i2<n2;i2++)
    for(i1=0;i1<n1;i1++)
    {
    rat[i3*n1*n2+i2*n1+i1]=dout[i3*n1pad*n2pad+i2*n1pad+i1];
    }
    	
    	 sf_warning("n1pad=%d, n2pad=%d, n3pad=%d",n1pad,n2pad,n3pad);
    	 sf_warning("nw1=%d, nw2=%d, nw3=%d",nw1, nw2, nw3);

    }
    

    sf_warning("ifsm=%d,ifnm=%d,ifns=%d,max(w)=%g",ifsm,ifnm,ifns,mmax_n);
    
    for (id=0; id < nd; id++) {
	remove = rat[id]*sig[id];
	noi2[id] = noi[id]-remove;
	sig2[id] = sig[id]+remove;
    }





    sf_floatwrite(rat,nd,fw);
    sf_floatwrite(noi2,nd,fnoi2);
    sf_floatwrite(sig2,nd,fsig2);

    exit(0);
}
