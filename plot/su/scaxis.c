/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

void scaxis (float x1     /* first x value */, 
	     float x2     /* second x value */, 
	     int *nxnum   /* number of numbered values */, 
	     float *dxnum /* increment between numbered values (dxnum>0.0) */, 
	     float *fxnum /* first numbered value */)
/*< compute a readable scale for use in plotting axes

scaxis attempts to honor the user-specified nxnum.  However, nxnum
will be modified if necessary for readability.  Also, fxnum and nxnum
will be adjusted to compensate for roundoff error; in particular, 
fxnum will not be less than xmin-eps, and fxnum+(nxnum-1)*dxnum 
will not be greater than xmax+eps, where eps = 0.0001*(xmax-xmin).
xmin is the minimum of x1 and x2.  xmax is the maximum of x1 and x2. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int n,i,iloga;
	float d,f,rdint[4],eps,a,b,xmin,xmax;

	/* set readable intervals */
	rdint[0] = 1.0;  rdint[1] = 2.0;  rdint[2] = 5.0;  rdint[3] = 10.0;

	/* handle x1==x2 as a special case */
	if  (x1==x2) {
		*nxnum = 1;
		*dxnum = 1.0;
		*fxnum = x1;
		return;
	}

	/* determine minimum and maximum x */
	xmin = (x1<x2)?x1:x2;
	xmax = (x1>x2)?x1:x2;
	
	/* get desired number of numbered values */
	n = *nxnum;
	n = (2>n)?2:n;
	
	/* determine output parameters, adjusted for roundoff */
	a = (xmax-xmin)/(float)(n-1);
	iloga = (int)log10(a);
	if (a<1.0) iloga = iloga - 1;
	b = a/pow(10.0,(double)iloga);
	for (i=0; i<3 && b>=sqrt(rdint[i]*rdint[i+1]); i++);
	d = rdint[i]*pow(10.0,(float)iloga);
	f = ((int)(xmin/d))*d-d;
	eps = 0.0001*(xmax-xmin);
	while(f<(xmin-eps))
		 f = f+d;
	n = 1+(int)((xmax+eps-f)/d); 
        
	/* set output parameters before returning */
	*nxnum = n;
	*dxnum = d;
	*fxnum = f;
}
