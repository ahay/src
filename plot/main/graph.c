#include <float.h>

#include <rsf.h>
#include <rsfplot.h>


int main(int argc, char* argv[])
{
    bool wantframe, wantframenum, transp, yreverse, pad, npad; 
    int nnx, fastplt=0, i, j, k, l, i3, i2, i1;
    int movie, wantlegend, legendsz,legendfat, sti2 = 0, js;
    int n1, n2, n3;
    char string[80],titles[1024],title_temp[128],title_out[128];
    char curvelabel[1024],legendloc[3];
    float min1, max1, min2, max2, o3, d3, o2, d2, f;
    float ch, vs, **x, **y, **tmp, xc, yc, xup, yup, xpath, ypath;    
    complex float** data;
    sf_file in;
    sf_datatype type;

    sf_init(argc,argv);
    in = sf_input("in");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);
    if (n3 > 1) {
	if (!sf_histfloat(in,"o3",&o3)) o3=0.;
	if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    }

    type = sf_gettype(in);
    if (SF_FLOAT != type && SF_COMPLEX != type)
	sf_error("wrong data type %d",type);

    x = sf_floatalloc2(n1,n2);
    y = sf_floatalloc2(n1,n2);

    if (SF_FLOAT == type) {
	if (!sf_histfloat(in,"o2",&o2)) o2=0.;
	if (!sf_histfloat(in,"d2",&d2)) d2=1.;

	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		x[i2][i1] = o2 + i2*d2;
	    }
	}
    } else {
	data = sf_complexalloc2(n1,n2);
    }

    if (!sf_getbool ("wantframe",&wantframe)) wantframe=true;
    if (!sf_getbool ("wantframenum",&wantframenum)) wantframenum=true;
    if (!sf_getbool ("transp",&transp)) transp = false;
    if (!sf_getbool ("yreverse",&yreverse)) yreverse = false;

    vp_coord_init (transp, yreverse);
    vp_plot_init (n2);
    vp_title_init(in);
    vp_color_init ();
    vp_erase ();

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
	    sf_read(y[i2],sizeof(float),n1*n2,in);
	    min1=o2;
	    max1=o2+(n2-1)*d2;
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    f = y[i2][i1];
		    if      (f > max2) max2=f;
		    else if (f < min2) min2=f;
		}
	    }
	}
	vp_minmax (min1, min2, max1, max2);
	
	if (transp) {
	    tmp = x;
	    x = y;
	    y = tmp;
	}

	npad = sf_getbool ("pad",&pad);
	vp_pad_init(pad, npad);
	
    }
    
    exit (0);
}


#ifdef lkjghkjhg

	minmax (x,y);
	if (n3_loop != 0){
    gl_erase ();
 		gl_plotint (&plot, &dash); /*reinitialize plot dash type */
	}
	gl_vplotint (&position, &coordinate, &axis1, &axis2);
	if (color.fillcol[0] != color.backcol[0] || color.fillcol[1] != color.backcol[1] || color.fillcol[2] != color.backcol[2])
	    gl_fillin (&coordinate, &color);

/* next section of code is commented out. No need for the extra padding on the clip window */
#if 0
	/*
	 * need some extra padding  on the clip window if plot.fat and or
	 * axisfat are greater then 0. 
	 */
	tempfat = 0;
	for (i = 0; i < data.n2; i++)
	{
	    if (tempfat < plot.fat[i])
		tempfat = plot.fat[i];
	}
	tempfat = tempfat + axis1.fat[n3_loop];
	tempfat = tempfat / 200;
	tempfat1 = 1- tempfat ;
	tempfat2 = tempfat + 1;
	gl_uclip (tempfat1 * coordinatetemp.min1, tempfat1 * coordinatetemp.min2, tempfat2 * coordinatetemp.max1, tempfat2 * coordinatetemp.max2);
#endif         

	gl_uclip (coordinatetemp.min1, coordinatetemp.min2,  coordinatetemp.max1, coordinatetemp.max2); 

	/* setting the symbol size */
	l = 0;
	for (i = sti2; i < data.n2 + sti2; i++)
	{
		js=i-sti2; /*THIS RESTARTS COLOR WITH EACH FRAME -added by bob 10/99*/
	    gl_dashfig (&dash, js);
	    gl_dash (&dash);
	    gl_fat (plot.fat[js]);
	    gl_color (plot.col[js]);
	    ch = plot.symbolsz[js] / 33.;
/* next section of code is commented out, this is not really necessary.  Trying
to use actual font size instead of trying to make symbols look the same size */ 

#if 0

	    if (plot.symbol[i] == '.')
		ch = ch * 2.0 / 10.0;
	    else
		if ((int) (plot.symbol[i]) >= (int) ('a') &&
		    (int) (plot.symbol[i]) <= (int) ('z'))
		ch = ch * 4.0 / 10.0;
	    else
		ch = ch * 6.0 / 10.0;
#endif
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
	if(multi_t ==1) {
	  if(!search_title(n3_loop, title_out,titles)) {
	    strcpy (title.title, title_out);
	  } else strcpy (title.title, title_temp);
	}
	gl_stdplot (&data, &coordinate, &axis1, &axis2, &grid, &title, n3_loop, fastplt, wantframe,wantframenum);  

    }				/* end of Movie loop */



static void minmax (int n1, int n2, float *point1, float *point2)
{
    int  i, j, k;
    float min, max;

    if (coordinate.yreverse)
    {
	gl_rotate (y, coordinate.min2, coordinate.max2, &data);
    }
    if (coordinate.xreverse)
    {
	gl_rotate (x, coordinate.min1, coordinate.max1, &data);
    }
    gl_axisint (&axis1, &axis2, &coordinate, &position);
    gl_gridint (&grid, &coordinate, &axis1, &axis2);
}


#endif
