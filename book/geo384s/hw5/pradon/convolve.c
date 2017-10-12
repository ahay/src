/* Match filtering */
#include <rsf.h>

int main(int argc, char* argv[])
{
    bool adj;
    int n1, n2, i1, i2, i, j, nf;
    float *data, *target, *filter;
    sf_file inp, out, other;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    other = sf_input("data");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	/* input data, output filter */
	if (!sf_histint(inp,"n1",&n1)) sf_error("No n1=");
	if (!sf_histint(inp,"n2",&n2)) sf_error("No n2=");
	if (!sf_getint("nf",&nf)) sf_error("Need nf=");
	/* filter size */

	sf_putint(out,"n1",nf);
	sf_putint(out,"n2",1);
    } else {
	/* input filter, output data */
	if (!sf_histint(inp,"n1",&nf)) sf_error("No n1=");
	if (!sf_histint(other,"n1",&n1)) sf_error("No n1=");
	if (!sf_histint(other,"n2",&n2)) sf_error("No n2=");

	sf_fileflush(out,other);  /* copy data dimensions */
    }

    filter = sf_floatalloc(nf);
    target = sf_floatalloc(n1);
    data = sf_floatalloc(n1);

    if (adj) {
	for (i=0; i < nf; i++) filter[i]=0.0f;
    } else {
	sf_floatread(filter,nf,inp);
    }

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,other);

	if (adj) {
	    sf_floatread(target,n1,inp);
	} else {
	    for (i1=0; i1 < n1; i1++) target[i1] = 0.0f;
	}
	
	for (i1=0; i1 < n1; i1++) {		
	    for (i=0; i < nf; i++) {
		
		j=i1-i-1; 

		/* zero value boundary conditions */
		if (j < 0 || j >= n1) continue; 
		    
		if (adj) {
		    filter[i] += data[j]*target[i1];
		} else {
		    target[i1] += data[j]*filter[i];
		}
	    }
	}

	if (!adj) sf_floatwrite(target,n1,out);
    }

    if (adj) sf_floatwrite(filter,nf,out);
		 
    exit(0);
}
