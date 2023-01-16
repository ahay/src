/* Mathematical operations on data files.
   
Known functions: 
cos,  sin,  tan,  acos,  asin,  atan, 
cosh, sinh, tanh, acosh, asinh, atanh,
exp,  log,  sqrt, abs,
erf,  erfc, sign (for float data),
arg,  conj, real, imag (for complex data).

sfmath will work on float or complex data, but all the input and output
files must be of the same data type.

An alternative to sfmath is sfadd, which may be more efficient, but is
less versatile.

Examples:

sfmath x=file1.rsf y=file2.rsf power=file3.rsf \
output='sin((x+2*y)^power)' > out.rsf
sfmath < file1.rsf tau=file2.rsf output='exp(tau*input)' > out.rsf
sfmath n1=100 type=complex output="exp(I*x1)" > out.rsf

Arguments which are not treated as variables in mathematical expressions:
datapath=, type=, out=

See also: sfheadermath.*/

/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* 	$Id$	 */

#include <string.h>
#include <ctype.h>
#include <math.h>

#include <unistd.h>

#include <rsf.h>

static void check_compat (size_t       nin, 
			  sf_file*     in, 
			  int          dim, 
			  const off_t* n, 
			  float*       d, 
			  float*       o, 
			  sf_datatype  type);

int main (int argc, char* argv[])
{
    int nin, i, i1, k, dim, nbuf;
    off_t j, nsiz, ii[SF_MAX_DIM], n[SF_MAX_DIM]; 
    size_t len, bufsiz;
    sf_file *in, out;
    char *eq, *output, *key, *arg, xkey[7], *ctype, *label, *unit;
    float **fbuf, **fst, d[SF_MAX_DIM], o[SF_MAX_DIM];
    bool nostdin;
    sf_complex **cbuf, **cst;
    sf_datatype type;


    /* init RSF */
    sf_init (argc,argv);
    
    in = (sf_file*) sf_alloc ((size_t) argc-1,sizeof(sf_file));    
    out = sf_output ("out");
    bufsiz = sf_bufsiz(out);

    if (!sf_getbool("nostdin", &nostdin)) nostdin=false;
    /* y - ignore stdin */

    /* find number of input files */
    if (nostdin || !sf_stdin()) { /* no input file in stdin */
	nin=0;
    } else {
	in[0] = sf_input("in");
	sf_putint(out,"input",0);
	nin=1;
    }
    
    for (i=1; i< argc; i++) { /* collect inputs */
	arg = argv[i];
	eq  = strchr(arg,'=');
	if (NULL == eq) continue; /* not a parameter */
    if (0 == strncmp(arg,"datapath=",9) ||
        0 == strncmp(arg, "nostdin=",8) ||
        0 == strncmp(arg,  "output=",7) ||
        0 == strncmp(arg,    "type=",5) ||
        0 == strncmp(arg,     "out=",4) ||
        0 == strncmp(arg,       "--",2) ||
	    (eq-arg == 2 &&
	     (arg[0] == 'n' || 
	      arg[0] == 'd' || 
	      arg[0] == 'o') &&
	     isdigit((int) arg[1])) ||
	    (0 == strncmp(arg,"label",5) &&
	     ((eq-arg == 5) || (eq-arg == 6 && isdigit((int) arg[5])))) ||
	    (0 == strncmp(arg,"unit",4) && 
	     ((eq-arg == 4) || (eq-arg == 5 && isdigit((int) arg[4]))))
	    ) continue; /* not a file */
	
	len = (size_t) (eq-arg);
	key = sf_charalloc(len+1);
	strncpy(key,arg,len);
	key[len]='\0';
	 
	in[nin] = sf_input(key);
	sf_putint(out,key,nin);
	nin++;
	free(key);
    }

    if (nin) {
	type = sf_gettype(in[0]);
	if (SF_FLOAT != type && SF_COMPLEX != type) 
	    sf_error("Need float or complex input");
	
	dim = sf_largefiledims(in[0],n);
	for (i=0; i < dim; i++) {
	    i1 = (i+1)%10u;
	    (void) snprintf(xkey,3,"d%d",i1);
	    if (!sf_histfloat(in[0],xkey,d+i)) d[i] = 1.;
	    (void) snprintf(xkey,3,"o%d",i1);
	    if (!sf_histfloat(in[0],xkey,o+i)) o[i] = 0.; 
	}
    } else { /* get type and size from parameters */
	
	ctype = sf_getstring("type");
	/* output data type [float,complex] */
	if (NULL != ctype) {	    
	    type = ('c'==ctype[0])? SF_COMPLEX:SF_FLOAT;
	} else {
	    type = SF_FLOAT;
	}
	
	dim = 1;
	for (i=0; i < SF_MAX_DIM; i++) {
	    i1 = (i+1)%10u;
	    
	    (void) snprintf(xkey,3,"n%d",i1);
	    if (!sf_getlargeint(xkey,n+i)) break;
            /*( n# size of #-th axis )*/
	    if (n[i] > 0) dim=i+1;
	    sf_putint(out,xkey,n[i]);

	    (void) snprintf(xkey,3,"d%d",i1);
	    if (!sf_getfloat(xkey,d+i)) d[i] = 1.;
            /*( d#=(1,1,...) sampling on #-th axis )*/
	    sf_putfloat(out,xkey,d[i]);

	    (void) snprintf(xkey,3,"o%d",i1);
	    if (!sf_getfloat(xkey,o+i)) o[i] = 0.;
            /*( o#=(0,0,...) origin on #-th axis )*/
	    sf_putfloat(out,xkey,o[i]);

	    (void) snprintf(xkey,7,"label%d",i1);
	    if (NULL != (label = sf_getstring(xkey))) {
                /*( label# label on #-th axis )*/
		sf_putstring(out,xkey,label);
		free(label);
	    }
	    if (NULL != (label = sf_getstring("label"))) {
		/* data label */
		sf_putstring(out,"label",label);
		free(label);
	    }
	    (void) snprintf(xkey,6,"unit%d",i1);
	    if (NULL != (unit = sf_getstring(xkey))) {
		/*( unit# unit on #-th axis )*/ 
		sf_putstring(out,xkey,unit);
		free(unit);
	    }
	    if (NULL != (unit = sf_getstring("unit"))) {
		/* data unit */
		sf_putstring(out,"unit",unit);
		free(unit);
	    }
	}
    }
    
    for (nsiz=1, i=0; i < dim; i++) {
	snprintf(xkey,3,"x%d",(i+1)%10u);
	
	if (NULL != sf_histstring(out,xkey)) 
	    sf_error("illegal use of %s parameter",xkey);
	sf_putint(out,xkey,nin+i);
	
	nsiz *= n[i];
    }
    
    if (SF_COMPLEX == type) sf_putint(out,"I",nin+dim);
    
    if (nin) check_compat(nin,in,dim,n,d,o,type);

    if (NULL == (output = sf_getstring("output"))) sf_error("Need output=");
    /* Mathematical description of the output */
    
    len = sf_math_parse (output,out,type);
    
    if (SF_FLOAT == type) { /* float type */
	nbuf = bufsiz/sizeof(float);
	
	fbuf = sf_floatalloc2(nbuf,nin+dim);
	fst  = sf_floatalloc2(nbuf,len+2);
	cbuf = NULL;
	cst  = NULL; 
    } else {                /* complex type */
	nbuf = bufsiz/sizeof(sf_complex);
	
	fbuf = NULL;
	fst  = NULL;
	cbuf = sf_complexalloc2(nbuf,nin+dim+1);
	cst  = sf_complexalloc2(nbuf,len+2);
    }

    if (nin) {
	sf_setformat(out,sf_histstring(in[0],"data_format"));    
	sf_fileflush(out,in[0]);
    } else {
	sf_settype(out,type);
	sf_setform(out,SF_NATIVE);
    }
    
    if (SF_FLOAT == type) {
	for (j=0; nsiz > 0; nsiz -= nbuf) {
	    if (nbuf > nsiz) nbuf = nsiz;
	    for (i=0; i < nin; i++) {
		sf_floatread(fbuf[i],nbuf,in[i]);
	    }
            /* fill x1, x2, etc */
	    for (k=0; k < nbuf; k++, j++) {
		sf_large_line2cart(dim,n,j,ii);
		for (i=0; i < dim; i++) {
		    fbuf[nin+i][k] = o[i]+ii[i]*d[i];
		}
	    }
	    
	    sf_math_evaluate (len, nbuf, fbuf, fst);
	    
	    sf_floatwrite(fst[1],nbuf,out);
	}
    } else {
	for (j=0; nsiz > 0; nsiz -= nbuf) {
	    if (nbuf > nsiz) nbuf = nsiz;
	    for (i=0; i < nin; i++) {
		sf_complexread(cbuf[i],nbuf,in[i]);
	    }
	    /* fill x1, x2, etc */
	    for (k=0; k < nbuf; k++, j++) {
		sf_large_line2cart(dim,n,j,ii);
		for (i=0; i < dim; i++) {
		    cbuf[nin+i][k] = sf_cmplx(o[i]+ii[i]*d[i],0.);
		}
		/* fill I */
		cbuf[nin+dim][k] = sf_cmplx(0.,1.);
	    }
	    
	    sf_complex_math_evaluate (len, nbuf, cbuf, cst);
	    
	    sf_complexwrite(cst[1],nbuf,out);
	}
    }
    

    exit(0);
}

/*------------------------------------------------------------*/
static void check_compat (size_t       nin, 
			  sf_file*     in, 
			  int          dim, 
			  const off_t* n, 
			  float*       d, 
			  float*       o, 
			  sf_datatype  type) 
{
    int ni, id;
    size_t i;
    float di, oi;
    char key[3];
    const float tol=1.e-5;
    
    for (i=1; i < nin; i++) {
	if (type != sf_gettype(in[i])) 
	    sf_error("Need %s input",(type==SF_FLOAT)?"float":"complex");
	for (id=0; id < dim; id++) {
	    (void) snprintf(key,3,"n%d",id+1);
	    if (!sf_histint(in[i],key,&ni) || ni != n[id])
		sf_error("%s mismatch: need %d",key,n[id]);
	    (void) snprintf(key,3,"d%d",id+1);
	    if (!sf_histfloat(in[i],key,&di) || 
		(fabsf(di-d[id]) > tol*fabsf(d[id])))
		sf_warning("%s mismatch: need %g",key,d[id]); 
	    (void) snprintf(key,3,"o%d",id+1);
	    if (!sf_histfloat(in[i],key,&oi) || 
		(fabsf(oi-o[id]) > tol*fabsf(d[id])))
		sf_warning("%s mismatch: need %g",key,o[id]);
	}
    }
}
