/* Pad a dataset with zeros.

Takes: < in.rsf > padded.rsf 

Parameters: 

[beg1= beg2= ... end1= end2=... | n1=  n2 = ... | n1out= n2out= ...]

*/

#include <stdio.h>
#include <string.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int i, j, nj, dim, ntr, itr, esize;
    int n[SF_MAX_DIM], n2[SF_MAX_DIM], beg[SF_MAX_DIM], end[SF_MAX_DIM];
    size_t n0, n20, beg0, end0;
    sf_file in, out;
    float o, d;
    char key[10], key2[10], *tr, *zero;
    bool inside;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    dim = sf_filedims(in,n);
    
    ntr=1;
    for (j=0; j < SF_MAX_DIM; j++) {
	i=j+1;
	if (j>=dim) {
	    n[j]=1;
	}
	sprintf(key,"beg%d",i);
	if(!sf_getint(key,beg+j)) {
	    beg[j]=0;
	} else if (beg[j]<0) {
	    sf_error("negative beg=%d",beg[j]);
	}
	sprintf(key,"end%d",i);
	sprintf(key2,"n%d",i);
	if(!sf_getint(key,end+j)) {
	    sprintf(key,"n%dout",i);
	    if (sf_getint(key,&nj) || sf_getint(key2,&nj)) {
		if (0==nj) for (nj++; nj < n[j]; nj *= 2) ;
		end[j]=nj-n[j]-beg[j];
		if (end[j]<0)
		    sf_error("negative end=%d",end[j]);
	    } else {
		end[j]=0;
	    }
	}
	n2[j]=n[j]+beg[j]+end[j];
	if (j>0) ntr *= n2[j];

	if (n2[j] != n[j]) sf_putint(out,key2,n2[j]);

	if (beg[j] > 0) {
	    sprintf(key,"o%d",i);
	    if (sf_histfloat(in,key,&o)) {
		sprintf(key,"d%d",i);
		if (sf_histfloat(in,key,&d))
		    sf_putfloat(out,key2,o-d*beg[j]);
	    }
	}
	if (n2[j] > 1 && j >= dim) dim=i;
    }

    if (!sf_histint(in,"esize",&esize)) {
	esize=4;
    } else if (0>=esize) {
	sf_error("wrong esize=%d",esize);
    }

    sf_fileflush(out,in);
    sf_setformat(in,"raw");
    sf_setformat(out,"raw");

    n[0]   *= esize; n0   = (size_t) n[0];
    n2[0]  *= esize; n20  = (size_t) n2[0];
    beg[0] *= esize; beg0 = (size_t) beg[0];
    end[0] *= esize; end0 = (size_t) end[0];

    tr   = sf_charalloc(n20);
    zero = sf_charalloc(n20);
    memset(zero,0,n20);

    if (beg0>0) memcpy(tr,zero,beg0);
    if (end0>0) memcpy(tr+beg0+n0,zero,end0);

    for (itr=0; itr < ntr; itr++) {
	inside = true;
	for (nj=j=1; j < dim; nj *= n2[j], j++) {
	    i = (itr/nj)%n2[j];
	    if (i < beg[j] || i >= beg[j]+n[j]) {
		inside = false;
		break;
	    }
	}
	if (inside) {
	    sf_read (tr+beg0,1,n0,in);
	    sf_write(tr,1,n20,out);
	} else {
	    sf_write(zero,1,n20,out);
	}
    }

    exit (0);
}

/* 	$Id: pad.c,v 1.3 2003/09/29 14:34:56 fomels Exp $	 */
