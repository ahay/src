/* Stack prespecified values. */
/*
 Copyright (C) 2011 King Abdullah University of Science & Technology
 
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
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int *im1, *im2;
    int i1,i2,i3, n1,n2,n3,n12, nm1,nm2=0,nm3, cnt; 
	
    float **dat, *min,*max, *stk;
    float o2,d2, maxval,thres;
	
    bool mean, minflag=false, maxflag=false;
    sf_file in, minin=NULL,maxin=NULL, out;
    
    sf_init (argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1))   sf_error("Need n1= for in");
    if (!sf_histint(in,"n2",&n2))   sf_error("Need n2= for in");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("Need d2= for in");
    if (!sf_histfloat(in,"o2",&o2)) o2 = 0.;
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;
	
    sf_putint(out,"n2",1);
	
    if (NULL != sf_getstring("min")) {
		/* file determining from which value to stack */
		minin = sf_input("min");
		if (SF_FLOAT != sf_gettype(minin)) sf_error("Need float input");
		if (!sf_histint(minin,"n1",&nm1))  sf_error("Need n1= for min=");
		minflag = true;
		
		if (n1!=nm1) sf_error("Wrong dimensions in min=.");
    }
	
    if (NULL != sf_getstring("max")) {
		/* file determining up to which value to stack */
		maxin = sf_input("max");
		if (SF_FLOAT != sf_gettype(maxin)) sf_error("Need float input");
		if (!sf_histint(maxin,"n1",&nm1))  sf_error("Need n1= for max=");
		maxflag = true;
		
		if (n1!=nm1) sf_error("Wrong dimensions in max=.");
    }
    if (minflag && maxflag) {
		nm2 = sf_leftsize(minin,1);
		nm3 = sf_leftsize(maxin,1);
		if (nm2!=nm3) sf_error("n2= of min= and max= do not agree");
    }
	
    if (!sf_getbool("mean",&mean)) mean=true;
    /* if n, sum; if y, average */
    if (!sf_getfloat("thres",&thres)) thres=0.;
    /* threshold (percentage of maxabs) */
    if (thres<0 || thres>1.) thres = 0.;
	
	
	
    dat = sf_floatalloc2(n1,n2);
    min = sf_floatalloc(n1);
    max = sf_floatalloc(n1);
    im1 = sf_intalloc(n1);
    im2 = sf_intalloc(n1);
    stk = sf_floatalloc(n1);
    
	
    sf_warning("minflag=%d, maxflag=%d",minflag,maxflag);
    for (i3=0; i3<n3; i3++) {
		sf_floatread(dat[0],n12,in);
		
		if (i3==0 || nm2>1) {	
		    if (minflag) {	/* read min values and transform them to indices */
				sf_floatread(min,n1,minin);
				for (i1=0; i1<n1; i1++) {
					min[i1] = SF_MIN(SF_MAX(o2,min[i1]), o2+(n2-1)*d2);
					im1[i1] = roundf((min[i1]-o2)/d2);
				}
			} else {
				for (i1=0; i1<n1; i1++) im1[i1] = 0;
			}
			
		    if (maxflag) {	/* read max values and transform them to indices */
				sf_floatread(max,n1,maxin);
				for (i1=0; i1<n1; i1++) {
					max[i1] = SF_MIN(SF_MAX(o2,max[i1]), o2+(n2-1)*d2);
					im2[i1] = roundf((max[i1]-o2)/d2);
				}
			} else {
				for (i1=0; i1<n1; i1++) im2[i1] = n2-1;
			}
			
			for (i1=0; i1<n1; i1++)
				if (min[i1]>max[i1]) min[i1]=max[i1];
		}	
		
		for (i1=0; i1<n1; i1++) {
			if (thres!=0.) {
				for (i2=im1[i1], maxval=0.; i2<=im2[i1]; i2++){
					if (maxval<fabsf(dat[i2][i1])) maxval = fabsf(dat[i2][i1]);
				}
				
				for (i2=im1[i1], stk[i1]=0., cnt=0; i2<=im2[i1]; i2++){
					if (fabsf(dat[i2][i1])>thres*maxval) {
						cnt++;
						stk[i1] += dat[i2][i1];
					}
				}
				if (mean) stk[i1] /= cnt;
			} else {
				for (i2=im1[i1], stk[i1]=0.; i2<=im2[i1]; i2++){
					stk[i1] += dat[i2][i1];
				}
				if (mean) stk[i1] /= (im2[i1]-im1[i1]+1);
			}
		}	
		
		sf_floatwrite(stk,n1,out);
    }
	
    exit(0);
}
