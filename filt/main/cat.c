#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <unistd.h>

#include <rsf.h>

static void check_compat (int esize, size_t nin, sf_file *in, int axis, int dim, 
			  const int *n, /*@out@*/ int *naxis);

int main (int argc, char* argv[])
{
    int i, axis, *naxis, n[SF_MAX_DIM], dim, n1, n2, i2, esize, nspace;
    size_t j, nin, ni, nbuf;
    sf_file *in, out;
    char* prog, key[3], buf[BUFSIZ];
    bool space;
    
    sf_init(argc,argv);
    
    in = (sf_file*) sf_alloc ((size_t) argc,sizeof(sf_file));

    if (0 != isatty(fileno(stdin))) { /* no input file in stdin */
	nin=0;
    } else {
	in[0] = sf_input("in");
	nin=1;
    }

    for (i=1; i< argc; i++) { /* collect inputs */
	if (NULL != strchr(argv[i],'=')) continue; /* not a file */
	in[nin] = sf_input(argv[i]);
	nin++;
    }
    if (0==nin) sf_error ("no input");

    out = sf_output ("out");

    if (!sf_getbool("space",&space)) {
	prog = sf_getprog();
	if (NULL != strstr (prog, "merge")) {
	    space = true;
	} else if (NULL != strstr (prog, "cat")) {
	    space = false;
	} else {
	    sf_warning("%s is neither merge nor cat, assume merge",prog);
	    space = true;
	}
    }

    dim = sf_filedims(in[0],n);
    if (!sf_getint("axis",&axis)) axis=3;
    if (1 > axis) sf_error("axis=%d < 1",axis);
    if (axis > dim) {
	while (dim < axis) {
	    n[dim++] = 1;
	}
    }

    n1=1;
    n2=1;
    for (i=1; i <= dim; i++) {
	if      (i < axis) n1 *= n[i-1];
	else if (i > axis) n2 *= n[i-1];
    }

    naxis = sf_intalloc(nin);
    
    if (!sf_histint(in[0],"esize",&esize)) {
	esize=4;
    } else if (0>=esize) {
	sf_error("wrong esize=%d",esize);
    }
    check_compat(esize,nin,in,axis,dim,n,naxis);

    /* figure out the length of extended axis */
    ni = 0;
    for (j=0; j < nin; j++) {
	ni += naxis[j];
	sf_setformat(in[j],"raw");
    }

    if (space) {
	if (!sf_getint("nspace",&nspace)) nspace = (int) (ni/(20*nin) + 1);
	ni += nspace*(nin-1);
    } 

    (void) snprintf(key,3,"n%d",axis);
    sf_putint(out,key,(int) ni);
    sf_fileflush(out,in[0]);
    sf_setformat(out,"raw");

    for (i2=0; i2 < n2; i2++) {
	for (j=0; j < nin; j++) {
	    for (ni = (size_t) n1*naxis[j]*esize; ni > 0; ni -= nbuf) {
		nbuf = (BUFSIZ < ni)? BUFSIZ: ni;
		sf_read (buf,1,nbuf,in[j]);
		sf_write (buf,1,nbuf,out);
	    }
	    if (!space || j == nin-1) continue;
	    /* Add spaces */
	    memset(buf,0,BUFSIZ);
	    for (ni = (size_t) n1*nspace*esize; ni > 0; ni -= nbuf) {
		nbuf = (BUFSIZ < ni)? BUFSIZ: ni;
		sf_write (buf,1,nbuf,out);
	    }
	}
    }
    
    sf_fileclose(out);
    exit(0);
}

static void check_compat (int esize, size_t nin, sf_file *in, int axis, int dim, 
			  const int *n, /*@out@*/ int *naxis) 
{
    size_t i;
    int ni, id;
    float f, fi;
    char key[3];
    const float tol=1.e-5;
    
    naxis[0] = n[axis-1];
    for (i=1; i < nin; i++) {
	if (!sf_histint(in[i],"esize",&ni) || ni != esize)
	    sf_error ("esize mismatch: need %d",esize);
	for (id=1; id <= dim; id++) {
	    (void) snprintf(key,3,"n%d",id);
	    if (!sf_histint(in[i],key,&ni) || (id != axis && ni != n[id-1]))
		sf_error("%s mismatch: need %d",key,n[id-1]);
	    if (id == axis) naxis[i] = ni;
	    (void) snprintf(key,3,"d%d",id);
	    if (sf_histfloat(in[0],key,&f) && 
		(!sf_histfloat(in[i],key,&fi) || 
		 (id != axis && fabsf(fi-f) > tol*fabsf(f))))
		sf_warning("%s mismatch: need %g",key,f);
	    (void) snprintf(key,3,"o%d",id);
	    if (sf_histfloat(in[0],key,&f) && 
		(!sf_histfloat(in[i],key,&fi) || 
		 (id != axis && fabsf(fi-f) > tol*fabsf(f))))
		sf_warning("%s mismatch: need %g",key,f);
	}
    }
}
