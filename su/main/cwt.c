/* succwt - Complex continuous wavelet transform of seismic traces */
/*
 * Credits: CWP: John Stockwell, Nov 2004
 * The code is borrowed from Mdwt.c under https://github.com/JohnWStockwellJr/SeisUnix
 * It was rewritten into Madagascar's environment
 *
 *
*/
/* math is included by rsf.h */
// #include <math.h>    /* pow, max, abs, sqrt */
#include <rsf.h>        /* hilbert_init, hilbert   */
#include <rsf_su.h>     /* convolve_cwp */

int
main(int argc, char *argv[])
{

	int i,j,k;		/* counters */
	int n1, n2;     /* size of input   */
	int ns=0;		/* number of samples in input data */
	int nwavelet=1024;	/* number of samples in mother wavelet */

	float base=0.0;		/* base */
	float first=0.0;	/* first exponent */
	float expinc=0.0;	/* exponent increment */
	float last=0.0;		/* last exponent */
	float exponent=0.0;	/* each exponent */
	float maxscale=0.0;	/* maximum scale value */
	float minscale=0.0;	/* minimum scale value */

	float x=0.0;
	float dx=0.0;		/* xvalues incr */

	float xmin=0.0;		/* last xvalues - first vval */
	float xcenter=0.0;	/* x value of center of wavelet */
	float xmax=0.0;		/* last xvalues - first vval */
	float sigma=1.0;	/* sharpening parameter */

	float waveletinc=0.0;		/* wavelet interval */
	float fmin=0.0;		/* min, max filt value (debug) */
	float *xvalues=NULL;	/* wavelet xvalues */
	float **filt=NULL;	/* filter used for each conv */

	float *f=NULL;		/* scratch for filter fliplr */
	float *sucwt_buff=NULL;	/* scratch for convolution */
	float *scales=NULL;	/* scales */
	float *waveletsum=NULL;	/* cumulative sum of wavelet */

	float *rt=NULL;		/* temp data storage */
	float *qt=NULL;		/* temp hilbert transformed data storage */
	float **tmpdata=NULL;	/* temp data storage */

	int wtype=0;		/* type of wavelet selected */
	float *wavelet=NULL;	/* pointer to data constituting the wavelet */

	int verbose=0;		/* verbose flag */
	int *index=NULL;	/* wavelet subscripts to use for filter */
	int *nconv=NULL;	/* length of each filter */
	int nscales=0;		/* number of scales */

	int holder=0;		/* =1 compute the Holder-Lipschitz regularity */
	float divisor=1.0;	/* divisor used in Holder exponent calculation*/
	
	float lrcoeff[4];
	int nh;             /* Hilbert transformer order */
	float ch;           /* Hilbert transformer reference (0.5 < ref <= 1) */
	
	float *xr;         /* compute the Holder regularity traces */
	float *yr;
    int icount;
	float maxalpha, interval;
	
	
	
	sf_file in,out;
	/*
	in: ns * n2
	out: 
	*/
	
	/* Initialize */
	sf_init(argc,argv);
	in = sf_input("in");
    out = sf_output("out");
	
	if (!sf_histint(in,"n1",&ns) sf_error("No n1= in input");
	n2 = sf_leftsize(in,1);
	
	
	
	/* Get parameters */
	if(!sf_getfloat ("base",&base))         base = 10;	
	if(!sf_getfloat("first",&first))		first = -1.0;
	if(!sf_getfloat("expinc",&expinc))		expinc = 0.01;
	if(!sf_getfloat("last",&last))			last = 1.5;

	if(verbose)
		sf_warning("base=%f, first=%f, expinc=%f, last=%f",base,first,expinc,last);
	
	// para for Ricker wavelet
	if(!sf_getint("wtype",&wtype))			wtype = 0;
	if(!sf_getint("nwavelet",&nwavelet))	nwavelet = 1024;
	if(!sf_getfloat("xmin",&xmin))			xmin = -20.0;
	if(!sf_getfloat("xcenter",&xcenter))	xmin = 0.0;
	if(!sf_getfloat("xmax",&xmax))			xmax = 20.0;
	if(!sf_getfloat("sigma",&sigma))		sigma = 1.0;

	if(!sf_getint("holder",&holder))		holder = 0;
	if(!sf_getfloat("divisor",&divisor))	divisor = 1.0;

	if(!sf_getint("verbose",&verbose))		verbose = 0;
	
	// para for Hilbert
	if (!sf_getint("order",&nh))            nh=100;    
    if (!sf_getfloat("ref",&ch))            ch=1.;
	
	xvalues = sf_floatalloc(nwavelet);
	wavelet = sf_floatalloc(nwavelet);
	
    
	/* Allocate space */	
	for (i=0; i< nwavelet; i++)
	{
		xvalues[i] = 0.0;
	    wavelet[i] = 0.0;
	}
	
	/* Compute wavelet */
	if (wtype == 0 ) 
	{ /* so far only Mex. Hat (Ricker) function */
		MexicanHatFunction(nwavelet, xmin, xcenter,
					xmax, sigma, wavelet);
	} else 		
	{
		sf_error("%d  type of wavelet not yet implemented",wtype); 
	}

	/* wavelet increment */
	waveletinc = (xmax - xmin)/(nwavelet - 1);
	
	/* verbose  warning */
	if(verbose)
	 sf_warning("xmin=%f, xmax=%f, nwavelet=%d, waveletinc=%f",
			xmin,xmax,nwavelet,waveletinc);

	/* form xvalues[] array */
	for(i=0,x=xmin; i<nwavelet; ++i,x+=waveletinc) xvalues[i] = x;

	xvalues[nwavelet-1] = xmax;

	/* compute scales */
	
	scales = sf_floatalloc(SHRT_MAX);   // SHRT_MAX = 32767
	for (i=0; i<SHRT_MAX; i++)  scales[i] = 0.0;

	exponent = first;
	x = 0;
	nscales = 0;
	minscale = pow(base,first);
	maxscale = pow(base,last);
	while(x <= maxscale) 
	{
		x = pow(base,exponent);
		scales[nscales] = x;
		exponent+=expinc;
		++nscales;

		if(nscales == SHRT_MAX)
			sf_error("Too many scales, change params and re-run\n");
	}
	--nscales;
	
	
	/* Allocate space */
	nconv      = sf_intalloc(nscales);
	index      = sf_intalloc(nwavelet);
	waveletsum = sf_floatalloc(nwavelet);
	filt       = sf_floatalloc2(nwavelet,nscales);
	f          = sf_floatalloc(nwavelet);

	/* Zero out arrays */
	for (i=0; i < nscales; i++)
		nconv[i] = 0.0;
	
	for (i=0; i < nwavelet; i++)
	{
		index[i] = 0;
		waveletsum[i] = 0.0;		
		f[i] = 0.0;
		for (j=0; j<nscales; j++) filt[i][j] = 0.0;
	}		
	

	/* Form difference of xvalues */
	for(i=nwavelet-1; i>=0; --i)
		xvalues[i] = xvalues[i] - xvalues[0];	

	dx = xvalues[1];
	xmax = xvalues[nwavelet-1];

	/* verbose warning */
	if(verbose) {
		sf_warning("first xvalues=%f, last xvalues=%f",
				xvalues[0],xvalues[nwavelet-1]);
		sf_warning("dx=%f, xmax=%f",dx,xmax);
	}
	
	/* waveletsum is cumulative sum of wavelet multipled by dx */
	fmin = 0;

	for(i=0; i<nwavelet; ++i) {
		fmin += wavelet[i];
		waveletsum[i] = fmin * dx;
	}

	/* Build filters from summed wavelet */
	for(i=0; i<nscales; ++i) {
		nconv[i] = 1 + (int)(scales[i] * xmax);

		for(j=0; j<nconv[i]; ++j)
			index[j] = 1 + j / (scales[i] * dx);

		for(j=0; j<nconv[i]; ++j)
			f[j] = waveletsum[index[j]-1];

		/* flip left right */
		for(j=0,k=nconv[i]-1; j<nconv[i]; ++j,--k)
			filt[i][j] = f[k];
	}

	/* Verbose warning */
	if(verbose) {
		sf_warning("Convolution Lengths");
		for(i=0; i<nscales; ++i) sf_warning("%d ",nconv[i]);
	}
	if(verbose) sf_warning("%d scales will be used for transforms",nscales);



	/* Allocate temporary storage space */
	rt = sf_floatalloc(ns);
	qt = sf_floatalloc(ns);
	tmpdata = sf_floatalloc(nscales,ns);

	/* Zero out rt and qt */
	for (i=0; i < ns; i++)
	{
		rt[i] = 0.0;
		qt[i] = 0.0;
	}

	/* Alloc sucwt_buffer for longest convolution */
	sucwt_buff = sf_floatalloc(ns+nconv[nscales-1]+1);
		
	sf_putint(out,"n2",nscales);
	sf_putfloat(out, "d2",waveletinc);
    sf_putfloat(out, "o2",minscale);
	sf_hilbert_init(ns, nh, ch);
		
	for (i2=0; i2<n2; i2++)   /* main loop over traces */
	{
		
		sf_floatread(trace,ns,in);

		

		
		/* Apply filters to produce wavelet transform */
		for(i=0; i<nscales; ++i) 
		{ /* loop over scales */

			for(j=0; j<ns+nconv[nscales-1]+1; ++j)
			sucwt_buff[j] = 0;

			/* convolve wavelet with data */
			convolve_cwp(ns,0,trace,nconv[i],0,
					filt[i],ns,0,sucwt_buff);

			for(j=0; j<ns; ++j) 
				rt[j] = sucwt_buff[j+nconv[i]/2-1];

			for(j=ns-1; j>0; --j) 
				rt[j] = rt[j] - rt[j-1];

			for(j=0; j<ns; ++j)
				rt[j] = -sqrt(scales[i]) * rt[j];

				/* form the hilbert transform of rt */			
			    /* inherit from https://reproducibility.org/blog/2011/11/05/program-of-the-month-sfenvelope/ */
			sf_hilbert(rt,qt);

			/* If not holder, then output envelope */
			if (!holder) 
			{
				
				for (j=0 ; j<ns; ++j) 
				{			
			  		output[i][j] = sqrt(rt[j]*rt[j] + qt[j]*qt[j]);
				}

			} 
			else 
			{
				/* compute the modulus */
				for (j=0 ; j<ns; ++j) {			
			  		tmpdata[j][i] = sqrt(rt[j]*rt[j] + qt[j]*qt[j]);
				}

			}
		}


		if (holder) 
		{ /* compute the Holder regularity traces */
			

			xr = sf_floatalloc(nscales);
			yr = sf_floatalloc(nscales);
			
	        /* Compute an estimate of the Lipschitz (Holder) regularity. Following Mallat (1992)	
            *				
            * ln | Wf(x,s)| <  ln C + alpha * ln|s|
            *					
            * alpha here is the Holder or Lipschitz exponent
            * s is the length scale, and Wf(x,s) is f in the
            * wavelet basis.			         
            *					         
			* Here we merely fit a straight line		 
			* through the log-log graph of the of our wavelet
			* transformed data and return the slope as      
			* the regularity measure. 			
			*/

          	for ( j =0 ; j< ns ; ++j ) 
        	{
        		icount=0;
                xr[0]=0;
                for ( i = 1 ; i < nscales ; ++i ) 
        		{
        
        	
        	          /* We stay away from values that will make */
          			  /*  NANs in the output */
                     if ((i>1) && (tmpdata[j][i-1] - tmpdata[j][1] > 0.0)) 
                     {
                        yr[icount] = log(abs(tmpdata[j][i] - tmpdata[j][1]));
                        xr[icount] = log(scales[i]-scales[1]);                      	
                      	++icount;
                      }

                   }
                 --icount;

               	/* straight line fit, return slope */
               	if ((icount> 10) && (divisor==1.0) ) 
        		{
                   linear_regression(yr, xr, icount, lrcoeff);
                           	   /* lrcoeff[0] is the slope of the line */
               	   /* which is the Holder (Lipschitz) */
               	   /* exponent */
        
                   output[i][j] = lrcoeff[0];
        
               	} else if ((icount> 10) && (divisor>1.0) ) 
        		{
        
               	   maxalpha=0.0;
               	   interval=icount/divisor;
        
               	   for ( k = interval; k < icount; k+=interval)
				   {
                   	   linear_regression(yr, xr, k, lrcoeff);
               		   maxalpha = max(lrcoeff[0],maxalpha);
               		}
               	   output[i][j] = maxalpha;		
        
               	} else if ((icount < 10) && (divisor>=1.0)) {
               	   output[i][j] = 0.0;		
               	} else if ( divisor < 1.0 ) {
               	   sf_error("divisor = %f < 1.0!", divisor);	
               	}

          	}

			 /* output holder regularity traces */
			sf_floatwrite(output[0],nscales*ns,out);
		}
	} 

	exit(0);
	

}


void
MexicanHatFunction(int nwavelet, float xmin, float xcenter, 
			float xmax, float sigma, float *wavelet) 
/***********************************************************************
MexicanHat - the classic Mexican hat function of length nwavelet 
************************************************************************
Input:
nwavelet	number of points total
xmin		minimum  x value
xcenter		central x value
xmax		maximum  x value
Output:
wavelet[]	array of floats of length nwavelet
************************************************************************
Notes:  

Gaussian:  g(x) = 1
		------------------- exp[ - ( x - xcenter)^2/2*sigma^2) ]
	       sigma*sqrt(2 * PI)

1st derivative of Gaussian:
g'(x) =  -(x - xcenter) 
	------------------- exp[ - (x - xcenter)^2/(2*sigma^2) ]
        (sigma^3)*sqrt(2 * PI)

Mexican Hat (2nd derivative of a Gaussian):
g''(x) = 1  
	--- 
      (sigma^3)* sqrt(2*PI)

           / ( x - xcenter)^2         \
	* |----------------   -  1  |
	   \ (sigma^2)              /
   
        * exp[ - (x - xcenter)^2/(2*sigma^2) ]

3rd derivative of Gaussian:
g'''(x) = 1  
	--- 
      (sigma^3)* sqrt(2*PI)

	* (x - xcenter)

           /      3 * ( x - xcenter)^3   \
	* |  1 -  ----------------        |
	   \         (sigma^2)           /
   
        * exp[ - (x - xcenter)^2/(2*sigma^2) ]

4th derivative of Gaussian (Witches' hat)
 (iv)
g  (x) =  1
	----
      (sigma^5)* sqrt(2*PI)
  
	   /      10 *( x - xcenter)^2         3 * ( x - xcenter)^4    \
        * | 1 -  -----------------------   +  ----------------------   |
           \           (sigma^2)                   (sigma^4)           /

        * exp[ - (x - xcenter)^2/(2*sigma^2) ]

************************************************************************
Author: CWP: John Stockwell (Nov 2004) 
************************************************************************/
{
	int i;			/* counter */
	double dxwavelet;	/* sampling interval in x */
	double mult;		/* multiplier on function */
	double twopi=2*SF_PI;	/* 2 PI */

	/* check for error in nwavelet */
	if (nwavelet<=1) sf_error("nwavelet must be greater than 1!");

	/* increment in x */
	dxwavelet = (xmax - xmin)/( nwavelet - 1);


	/* generate mexican hat function */
        /* ... multiplier ....*/
        mult = (1.0/(sigma*sigma*sigma * sqrt(twopi)));
	for(i=0; i<nwavelet; ++i) {
		float x = i*dxwavelet - xcenter; 

		wavelet[i] = mult * (x*x/(sigma*sigma) - 1.0) 
				* exp(- x*x/(2.0*sigma*sigma) );
	}
}
		
void linear_regression(float *y, float *x, int n, float coeff[4])
/*****************************************************************************
Compute linear regression of (y1,y2,...,yn) against (x1,x2,...,xn)
******************************************************************************
Input:
y		array of y values
x		array of x values
n		length of y and x arrays
Output:
coeff[4] where:

coeff[0]	slope of best fit line
coeff[1]	intercept of best fit line
coeff[2]	correlation coefficient of fit (1 = perfect) [dimensionless]
coeff[3]	standard error of fit (0 = perfect) [dimensions of y]
******************************************************************************
Notes: 

y(x) 
    |      *  .    fit is  y(x) = a x + b
    |       .          
    |     .  *
    | * .    
    | . *         
     ------------------- x
     
         n Sum[x*y] - Sum[x]*Sum[y]
     a = --------------------------
         n Sum[x*x] - Sum[x]*Sum[x]
         
         Sum[y] - a*Sum[x]
     b = -----------------
                n
                
     cc = std definition
     
     se = std definition
    
******************************************************************************
Author:  Chris Liner, UTulsa, 11/16/03
******************************************************************************/
{
        /* local variables */
	float den;	/* generic denomenator */
	float num;	/* generic numerator */
    float sx;	/* sum of x values */
    float sx2;	/* sum of x-squared values */
    float sy;	/* sum of y values */
    float sy2;	/* sum of y-squared values */
    float sxy;	/* sum of x*y values */
    float tmp;	/* temporary variable */
    int i;		/* counter */

	float a;	/* slope of best fit line */
	float b;	/* intercept of best fit line */
	float cc;	/* correlation coefficient of fit */
				/* (1 = perfect) [dimensionless] */
	float se;	/* standard error of fit (0 = perfect)  */
				/* [dimensions of y] */
    
    /* initialize sums */
    sx = x[0];
    sx2 = x[0]*x[0];
    sy = y[0];
    sy2 = y[0]*y[0];
    sxy = x[0]*y[0];
    
    /* calculate sums */
    for (i=1;i<n;++i) {
        sx += x[i]; 
        sx2 += x[i]*x[i];
        sy += y[i];
        sy2 += y[i]*y[i];
        sxy += x[i]*y[i];
    }
    
    /* slope */
    num = n*sxy - sx*sy;
    den = n*sx2 - sx*sx;
    a = num/den;

	coeff[0] = a;
        
    /* intercept */
    b = (sy - a*sx)/n;

	coeff[1] = b;
        
   /* correlation coefficient */
   num = abs(n*sxy - sx*sy);
   den = (n*sx2 - sx*sx)*(n*sy2 - sy*sy);
   den = sqrt(den);

   if (den != 0.0) 
   {
       cc = num/den;
   } else 
   {
       cc = 999;
   }

	coeff[2] = cc;
        
    /* standard error */
    tmp = 0.0;
    for (i=0;i<n;++i) 
	{
      tmp += (y[i] - (a*x[i] + b) )*(y[i] - (a*x[i] + b) ); 
    }

    se = sqrt( tmp / (n - 2) );
	coeff[3] = se;
}		
