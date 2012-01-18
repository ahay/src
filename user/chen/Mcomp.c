/* Compare 2 data set */

#include <rsf.h>
#include "comp.h"

int main(int argc,char*argv[])
{
	sf_file in,out,ref;
	int dim, dimr, n[SF_MAX_DIM], nr[SF_MAX_DIM];
	int n1, n2, id, i2,  mode;
	float *u,*ur,*p;

	sf_init(argc,argv);

	in  = sf_input("in");
	out = sf_output("out");

	ref = sf_input("ref");

	dimr = sf_filedims(ref,nr);
	dim  = sf_filedims(in,n);
	if(dim < dimr) sf_error("dim < dim_ref");

	for(id=0, n1=1; id<dimr; id++)
	{
		if(n[id] != nr[id]) 
			sf_error("Dimension %d not match!", id);
		n1 *= nr[id];
	}
	sf_unshiftdim(in,out,1);
	for(id=1; id<dimr; id++)	sf_unshiftdim(out,out,1);

	for(id=dimr, n2=1; id<dim; id++)	n2 *= n[id];

	if (!sf_getint("mode",&mode)) mode=0;
	/* compare method: 0 - unit xcorrelation; 1 - relative error */

	u  = sf_floatalloc(n1);
	ur = sf_floatalloc(n1);
	p  = sf_floatalloc(n2);

	sf_floatread(ur,n1,ref);

	comp_init(n1,ur);

	for(i2=0; i2<n2; i2++)
	{
		sf_floatread(u,n1,in);
		switch(mode)
		{
		case 1:
		default:
			p[i2]=comp_xcor(u);
		}
	}
	
	sf_floatwrite(p, n2, out);

	exit(0);
}

