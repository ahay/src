/* Mathematical operations on data files.

Example:

sfmath x=file1.rsf y=file2.rsf power=file3.rsf output='sin((x+2*y)^power)' > out.rsf
sfmath < file1.rsf tau=file2.rsf output='exp(tau*input)' > out.rsf

Known functions: cos, sin, tan, acos, asin, atan, exp, log, sqrt, abs.

See also: sfheadermath.
*/

#include <string.h>
#include <math.h>

#include <unistd.h>

#include <rsf.h>

static void check_compat (size_t nin, sf_file *in, int dim, const int *n);

int main(int argc, char* argv[])
{
    int nin, i, n[SF_MAX_DIM], dim, nbuf, nsiz;
    size_t len;
    sf_file *in, out;
    char *eq, *output, *key, *arg;
    float **fbuf, **fst;

    sf_init (argc,argv);
    
    in = (sf_file*) sf_alloc ((size_t) argc-1,sizeof(sf_file));    
    out = sf_output ("out");
    
    if (0 != isatty(fileno(stdin))) { /* no input file in stdin */
	nin=0;
    } else {
	in[0] = sf_input("in");
	sf_putint(out,"input",0);
	nin=1;
    }

    for (i=1; i< argc; i++) { /* collect inputs */
	arg = argv[i];
	eq =  strchr(arg,'=');
	if (NULL == eq) continue; /* not a parameter */
	if (0 == strncmp(arg,"output",6)) continue; /* not a file */
	
	len = (size_t) (eq-arg);
	key = sf_charalloc(len+1);
	strncpy(key,arg,len);
	key[len]='\0';
	 
	in[nin] = sf_input(key);
	sf_putint(out,key,nin);
	nin++;
	free(key);
    }

    dim = sf_filedims(in[0],n);
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }
    
    check_compat(nin,in,dim,n);

    if (NULL == (output = sf_getstring("output"))) sf_error("Need output=");
    /* Mathematical description of the output */

    len = sf_math_parse (output,out);
    
    nbuf = BUFSIZ/sizeof(float);

    fbuf = sf_floatalloc2(nbuf,nin);
    fst  = sf_floatalloc2(nbuf,len+2);

    sf_setformat(out,sf_histstring(in[0],"data_format"));    
    sf_fileflush(out,in[0]);

    for (; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf = nsiz;
	for (i=0; i < nin; i++) {
	    sf_floatread(fbuf[i],nbuf,in[i]);
	}

	sf_math_evaluate (len, nbuf, fbuf, fst);

	sf_floatwrite(fst[1],nbuf,out);
    }
    
    exit(0);
}

static void check_compat (size_t nin, sf_file *in, int dim, const int *n) 
{
    int ni, id;
    size_t i;
    float d, di, o, oi;
    char key[3];
    const float tol=1.e-5;
    
    if (SF_FLOAT != sf_gettype(in[0])) 
	sf_error("Need float input");
    for (i=1; i < nin; i++) {
	if (SF_FLOAT != sf_gettype(in[i])) 
	    sf_error("Need float input");
	for (id=1; id <= dim; id++) {
	    (void) snprintf(key,3,"n%d",id);
	    if (!sf_histint(in[i],key,&ni) || ni != n[id-1])
		sf_error("%s mismatch: need %d",key,n[id-1]);
	    (void) snprintf(key,3,"d%d",id);
	    if (sf_histfloat(in[0],key,&d)) { 
		if (!sf_histfloat(in[i],key,&di) || 
		    (fabsf(di-d) > tol*fabsf(d)))
		    sf_warning("%s mismatch: need %g",key,d);
	    } else {
		d =1.;
	    }
	    (void) snprintf(key,3,"o%d",id);
	    if (sf_histfloat(in[0],key,&o) && 
		(!sf_histfloat(in[i],key,&oi) || (fabsf(oi-o) > tol*fabsf(d))))
		sf_warning("%s mismatch: need %g",key,o);
	}
    }
}

/* 	$Id: math.c,v 1.9 2004/06/23 18:30:00 fomels Exp $	 */
