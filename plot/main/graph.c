#include <float.h>

#include <rsf.h>
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    int n1, n2, n3, i1, i2, i3;
    float min1, max1, min2, max2, o3, d3, o1, d1, **x, **y, f;    
    complex float** data;
    sf_datatype type;
    sf_file in;

    sf_init(argc,argv);
    in = sf_input("in");

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
		    if      (f > max1) max1=f;
		    else if (f < min1) min1=f;
		    x[i2][i1] = f;
		    f = cimagf(data[i2][i1]);
		    if      (f > max2) max2=f;
		    else if (f < min2) min2=f;
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
		    if      (f > max2) max2=f;
		    else if (f < min2) min2=f;
		}
	    }
	}

	vp_stdplot_init (min1, max1, min2, max2,
			 false,false,false,true);
	vp_frame_init(in,"blt");

	if (i3 > 0) vp_erase();
	vp_frame();

	for (i2=0; i2 < n2; i2++) {
	    vp_plot_set (i2);
	    vp_umove(x[i2][0],y[i2][0]);
	    
	    for (i1=1; i1 < n1; i1++) {
		vp_udraw(x[i2][i1],y[i2][i1]);
	    }
	}
    }
    
    exit (0);
}


#ifdef lkjghkjhg

	    gl_fat (plot.fat[js]);
	    gl_color (plot.col[js]);
	    ch = plot.symbolsz[js] / 33.;

	    gl_tjust ("s");
	    gl_penup ();
	    for (j = 0; j < data.n1[i]; j++)
	    {
		if (plot.symbol[i] != ' ')
		{
		    sprintf (string, "%c", plot.symbol[i]);
		    xpath = ch;
		    ypath = 0;
		    xup = 0;
		    yup = ch;
		    gl_umove (x[l], y[l]);
		    gl_where (&xc, &yc);
		    gl_gtext (xc, yc, xpath, ypath, xup, yup, string, "s");

		    if (plot.lineunder)
		    if(coordinate.transp)
		    {
			gl_umove (0., y[l]);
			gl_udraw (x[l], y[l]);
		    }
		    else
		    {
			gl_umove (x[l], 0.);
			gl_udraw (x[l], y[l]);
		    }
		}
		else
		{
		    gl_upendn (x[l], y[l]);
		}
		l = l + 1;
	    }
	    gl_purge ();
	}
	sti2 = i;
	dash.dash[0] = 0.;
	dash.dash[1] = 0.;
	dash.gap[0] = 0.;
	dash.gap[1] = 0.;
	gl_dash (&dash);
	gl_uclip (coordinate.min1, coordinate.min2, coordinate.max1, coordinate.max2);
	gl_stdplot (&data, &coordinate, &axis1, &axis2, &grid, &title, n3_loop, fastplt, wantframe,wantframenum);  

#endif
