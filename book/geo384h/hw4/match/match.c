/* Match filtering */
#include <rsf.h>

int main(int argc, char* argv[])
{
    bool adj;
    int n1, n2, i1, i2, i, j, nf;
    float *data, *noiz, *filt;
    sf_file inp, out, oth;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    oth = sf_input("other");

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
	if (!sf_histint(oth,"n1",&n1)) sf_error("No n1=");
	if (!sf_histint(oth,"n2",&n2)) sf_error("No n2=");

	sf_fileflush(out,oth);  /* copy data dimensions */
    }

    filt = sf_floatalloc(nf);
    data = sf_floatalloc(n1);
    noiz = sf_floatalloc(n1);

    if (adj) {
	for (i=0; i < nf; i++) /* !!! COMPLETE LINE !!! */
    } else {
	sf_floatread(filt,nf,inp);
    }

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(noiz,n1,oth);

	if (adj) {
	    sf_floatread /* !!! COMPLETE LINE !!! */
	} else {
	    for (i1=0; i1 < n1; i1++) data[i1] = 0.;
	}

	for (i=0; i < nf; i++) {
	    for (i1=0; i1 < n1; i1++) {
		j=i1-i+nf/2; /* symmetric filter */

		/* zero value boundary conditions */
		if (j < 0 || j >= n1) continue; 
		    
		if (adj) {
		    filt[i] += /* !!! COMPLETE LINE !!! */
		} else {
		    data[i1] += noiz[j]*filt[i];
		}
	    }
	}

	if (!adj) sf_floatwrite(data,n1,out);
    }

    if (adj) sf_floatwrite  /* !!! COMPLETE LINE !!! */
		 
    exit(0);
}
