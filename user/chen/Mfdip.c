/* 3D fast dip estimation by plane wave destruction */

#include <rsf.h>
#include "fdip.h"
#include <time.h>

int main(int argc,char**argv)
{
	int dim,n[SF_MAX_DIM], rect[3], n123, n4, nr, ir, j, liter;
	bool verb;
	float *u,*p;
	char key[4];
	sf_file in, out;
	clock_t start,end;

	sf_init(argc,argv);
	in = sf_input ("in");
	out = sf_output ("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	dim = sf_filedims(in,n);
	if (dim < 2) n[1]=1;
	if (dim < 3) n[2]=1;
	n123 = n[0]*n[1]*n[2];
	nr = 1;
	for (j=3; j < dim; j++) 
		nr *= n[j];
	
	if (1 == n[2]) {
		n4=0;
	} else {
		if (!sf_getint("n4",&n4)) n4=2;
		/* what to compute in 3-D. 0: in-line, 1: cross-line, 2: both */ 
		if (n4 > 2) n4=2;
		if (2==n4) {
			sf_putint(out,"n4",2);
			for (j=3; j < dim; j++) {
				snprintf(key,4,"n%d",j+2);
				sf_putint(out,key,n[j]);
			}
		}
	}

	if (!sf_getint("liter",&liter)) liter=20;
	/* number of linear iterations */

	if (!sf_getint("rect1",&rect[0])) rect[0]=1;
	/* dip smoothness on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=1;
	/* dip smoothness on 2nd axis */
	if (!sf_getint("rect3",&rect[2])) rect[2]=1;
	/* dip smoothness on 3rd axuis */
	if (!sf_getbool("verb",&verb)) verb = false;
	/* verbosity flag */


	u=sf_floatalloc(n123);
	if(n4==0 || n4==1)	p=sf_floatalloc(n123);
	else	p=sf_floatalloc(n123*2);

	fdip3_init(n[0], n[1], n[2], rect, liter, verb);
	
	start=clock();
	for(ir=0;ir<nr;ir++){
		sf_warning(" %d/%d;", ir, nr);
		sf_floatread(u,n123,in);
		fdip3( u, p, n4);
		if( n4 == 2 )	sf_floatwrite(p, n123*2, out);
		else	sf_floatwrite(p, n123, out);
	}
	end=clock();

	sf_warning("time cost=%f", (double)(end-start)/CLOCKS_PER_SEC);
	fdip3_close();
	exit(0);
}



