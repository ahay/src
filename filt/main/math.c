#include <string.h>
#include <math.h>

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

    nin=0;
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

    len = sf_math_parse (output,out);
    
    nbuf = BUFSIZ/sizeof(float);

    fbuf = sf_floatalloc2(nbuf,nin);
    fst  = sf_floatalloc2(nbuf,len+2);
    
    sf_fileflush(out,in[0]);
    sf_setformat(out,sf_histstring(in[0],"data_format"));

    for (; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf = nsiz;
	for (i=0; i < nin; i++) {
	    sf_read(fbuf[i],sizeof(float),nbuf,in[i]);
	}

	sf_math_evaluate (len, nbuf, fbuf, fst);

	sf_write(fst[1],sizeof(float),nbuf,out);
    }
    
    exit(0);
}

static void check_compat (size_t nin, sf_file *in, int dim, const int *n) 
{
    int ni, id;
    size_t i;
    float f, fi;
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
	    if (sf_histfloat(in[0],key,&f) && 
		(!sf_histfloat(in[i],key,&fi) || (fabsf(fi-f) > tol*fabsf(f))))
		sf_warning("%s mismatch: need %g",key,f);
	    (void) snprintf(key,3,"o%d",id);
	    if (sf_histfloat(in[0],key,&f) && 
		(!sf_histfloat(in[i],key,&fi) || (fabsf(fi-f) > tol*fabsf(f))))
		sf_warning("%s mismatch: need %g",key,f);
	}
    }
}
