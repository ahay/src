/* Gaussian wavelet estimation in 2-D.

Takes: < data.rsf ma=ma.rsf > monof.rsf
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "monof2.h"

int main(int argc, char* argv[])
{
    int n2, i2, nx, ny, niter;
    float x0, dx, y0, dy, a[3], **data;

    bool verb;
    sf_file in, out, ma;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    ma = sf_output("ma");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ny)) sf_error("No n2= in input");
    n2 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d2",&dy)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&y0)) sf_error("No o2= in input");

    if (!sf_getfloat("a0",a)) a[0]=1.;
    /* starting sharpness in xx */
    if (!sf_getfloat("b0",a+1)) a[1]=0.;
    /* starting sharpness in xy */
    if (!sf_getfloat("c0",a+2)) a[2]=1.;
    /* starting sharpness in yy */

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(ma,"n1",3);
    sf_putint(ma,"n2",1);
    sf_putint(ma,"nx",nx);
    sf_putfloat(ma,"dx",dx);
    sf_putfloat(ma,"x0",x0);
    sf_putint(ma,"ny",ny);
    sf_putfloat(ma,"dy",dy);
    sf_putfloat(ma,"y0",y0);
    sf_fileflush(ma,in);

    data = sf_floatalloc2(nx,ny);

    for (i2=0; i2 < n2; i2++) {
	sf_read(data[0],sizeof(float),nx*ny,in);

	monof2(data,niter,a,nx,dx,x0,ny,dy,y0,verb);
         
	sf_write(a,sizeof(float),3,ma);
        
	sf_write (data[0],sizeof(float),nx*ny,out);
    }
    
    exit (0);
}

/* 	$Id: Mmonof2.c,v 1.1 2004/02/25 16:16:27 fomels Exp $	 */
