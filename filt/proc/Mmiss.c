#include <rsf.h>

#include "helix.h"
#include "mis2.h"
#include "bound.h"

int main(int argc, char* argv[])
{
    int dim, na,ia, niter, n1, n2, maxlag, padin, padout, p1, p2, j, i;
    float a0, *mm, *kk;
    bool prec, *known;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], a[SF_MAX_DIM];
    filter aa;
    char *lagfile;
    sf_file in, out, filt, lag, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    out = sf_output("out");

    dim = sf_filedims(in,n);
    sf_putint (out,"dim",dim);

    if (!sf_getbool("prec",&prec)) prec=true;
    if (!sf_getint("niter",&niter)) niter=100;

    if (!sf_getint("padin",&padin)) padin=0;   
    if (!sf_getint("padout",&padout)) padout=0;
    n[dim-1] += padin + padout;

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");
    aa = allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
    if (!sf_histints(filt,"n",m,dim)) sf_error("No n= in filt");

    if (!sf_histints(filt,"a",a,dim)) {
	for (j=0; j < dim; j++) {
	    a[j]=1;
	}
    }
    if (NULL == (lagfile = sf_histstring(filt,"lag"))) {
	if (NULL == (lagfile = sf_getstring("lag"))) {
	    for (ia=0; ia < na; ia++) {
		aa->lag[ia]=ia;
	    }
	    lag = NULL;
	} else {
	    lag = sf_input("lag");
	}
	lag = sf_input(lagfile);
    }

    if (NULL != lag) {
	if (SF_INT != sf_gettype(lag)) 
	    sf_error("Need int data in lag file '%s'",lagfile);
	sf_read(aa->lag,sizeof(int),na,lag);
	sf_fileclose(lag);
    }

    bound (dim, m, n, a, aa);

    sf_read(aa->flt,sizeof(float),na,filt);
    sf_fileclose(filt);

    for (ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    n1=1;
    for (j=0; j < dim-1; j++) {
	n1 *= n[j];
    }
    n2 = n1*n[dim-1]; 

    p1 = padin*n1; 
    p2 = n2 - (padin+padout)*n1;

    mm = sf_floatalloc(n2);
    kk = sf_floatalloc(n2);
    known = sf_boolalloc(n2);

    for (i=0; i < n2; i++) {
	mm[i]=0.;
	kk[i]=0.;
	known[i]=false;
    }

    sf_read(mm+p1,sizeof(float),p2,in);

    if (NULL != sf_getstring("mask")) {
	mask = sf_input("mask");
	sf_read(kk+p1,sizeof(float),p2,mask);
	sf_fileclose(mask);
	
	for (i=p1; i < p1+p2; i++) {
	    known[i] = (kk[i] != 0.);
	}
    } else {
	for (i=p1; i < p1+p2; i++) {
	    known[i] = (mm[i] != 0.);
	}
    }

    maxlag=0;
    for (ia=0; ia < na; ia++) {
	if (aa->lag[ia] > maxlag) maxlag=aa->lag[ia];
    }

    if (padout*n1 >= maxlag) {
	for (i=n2-maxlag; i < n2; i++) {
	    known[i] = true;
	}
    }

    if (prec) {
	for (i=p1; i < p1+p2; i++) {
	    if (known[i]) kk[i]=mm[i];
	}
    }

    mis2 (niter, n2, mm, aa, known, prec);

    if (prec) {
	for (i=p1; i < p1+p2; i++) {
	    if (known[i]) mm[i]=kk[i];
	}
    }

    sf_write(mm+p1,sizeof(float),p2,out);

    exit (0);
}

