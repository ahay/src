/* Graph plot.

Takes: < data.rsf > plot.vpl
*/

#include <string.h>
#include <float.h>

#include <rsf.h>
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    bool transp;
    int n1, n2, n3, i1, i2, i3, len;
    float min1, max1, min2, max2, o3, d3, o1, d1;
    float **x, **y, **tmp, f, *symbolsz=NULL, symsize, xc, yc;    
    complex float** data=NULL;
    char* symbol, sym[2]=" ";
    sf_datatype type;
    sf_file in;

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);
    if (n3 > 1) {
	if (!sf_histfloat(in,"o3",&o3)) o3=0.;
	if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    }

    x = sf_floatalloc2(n1,n2);
    y = sf_floatalloc2(n1,n2);

    type = sf_gettype(in);
    switch (type) {
	case SF_FLOAT:
	    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
	    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	    
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    x[i2][i1] = o1 + i1*d1;
		}
	    }
	    break;
	case SF_COMPLEX:
	    data = sf_complexalloc2(n1,n2);
	    break;
	default:
	    sf_error("Wrong data type (need float or complex)");
    }

    vp_plot_init(n2);

    symbol = sf_getstring("symbol");
    /* if set, plot with symbols instead of lines */
    if (NULL != symbol) {
	len = strlen(symbol);
	if (len < n2) {
	    symbol = (char*) sf_realloc(symbol,n2,sizeof(char));
	    for (i2=len; i2 < n2; i2++) {
		symbol[i2] = symbol[i2 % len];
	    }
	}

	symbolsz = sf_floatalloc(n2);
	if (!sf_getfloats("symbolsz",symbolsz,n2)) {
	    for (i2 = 0; i2 < n2; i2++)
		symbolsz[i2] = 2./33.;
	} else {
	    for (i2 = 0; i2 < n2; i2++)
		symbolsz[i2] /= 33.;
	}
    }

    if (!sf_getbool ("transp",&transp)) transp=false;
    /* if y, transpose the axes */
 
    for (i3 = 0; i3 < n3; i3++) {
	min2 = +FLT_MAX;
	max2 = -FLT_MAX;
	if (SF_COMPLEX == type) {
	    min1 = min2;
	    max1 = max2;
	    sf_read(data[0],sizeof(float complex),n1*n2,in);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    f = crealf(data[i2][i1]);
		    if (f > max1) max1=f;
		    if (f < min1) min1=f;
		    x[i2][i1] = f;
		    f = cimagf(data[i2][i1]);
		    if (f > max2) max2=f;
		    if (f < min2) min2=f;
		    y[i2][i1] = f;
		}
	    }
	} else {
	    sf_read(y[0],sizeof(float),n1*n2,in);
	    min1=o1;
	    max1=o1+(n1-1)*d1;
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    f = y[i2][i1];
		    if (f > max2) max2=f;
		    if (f < min2) min2=f;
		}
	    }
	}

	if (transp) {tmp=x; x=y; y=tmp;}

	vp_stdplot_init (min1, max1, min2, max2,
			 transp,false,false,true);
	vp_frame_init(in,"blt");

	if (i3 > 0) vp_erase();
	vp_frame();

	for (i2=0; i2 < n2; i2++) {
	    vp_plot_set (i2);
	    
	    if (NULL != symbol) {
		sym[0] = symbol[i2];
		symsize = symbolsz[i2];

		for (i1=0; i1 < n1; i1++) {
		    vp_umove(x[i2][i1],y[i2][i1]);
		    vp_where (&xc, &yc);
		    vp_tjust (TH_SYMBOL, TV_SYMBOL);
		    vp_gtext (xc,yc,symsize,0.,0.,symsize,sym);
		}
	    } else {
		vp_umove(x[i2][0],y[i2][0]);	    
		for (i1=1; i1 < n1; i1++) {
		    vp_udraw(x[i2][i1],y[i2][i1]);
		}
	    }
	}

	if (transp) {tmp=x; x=y; y=tmp;}
    } 
   
    sf_close();
    exit (0);
}

/* 	$Id: graph.c,v 1.10 2004/03/22 05:43:25 fomels Exp $	 */
