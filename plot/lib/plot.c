#include <math.h>

#include <rsf.h>

#include "vplot.h"

static float *dashtype;
static int *fat, *col;

static void set_dash (float type);

/*
 *	0 continuous   DEFAULT
 *	1 fine dash
 *	2 fine dot
 *	3 dash
 *	4 large dash
 *	5 dot dash
 *	6 large dash small dash
 *	7 double dot
 *	8 double dash
 *	9 loose dash  The part after the decimal point determines 
 *                     the pattern repetition interval
 */
static void set_dash (float type)
{
    float size;
    float dash[2], gap[2];
    
    switch ((int) type) {
	case 1:
	    size = .3;
	    dash[0] = dash[1] = 1.;
	    gap[0] = gap[1] = 2.;
	    break;
	case 2:
	    size = .2;
	    dash[0] = dash[1] = 1.;
	    gap[0] = gap[1] = 6.;
	    break;
	case 3:
	    size = .4;
	    dash[0] = dash[1] = 4.;
	    gap[0] = gap[1] = 2.;
	    break;
	case 4:
	    size = .6;
	    dash[0] = dash[1] = 3.;
	    gap[0] = gap[1] = 2.;
	    break;
	case 5:
	    size = .5;
	    dash[0] = .3;
	    dash[1] = 3.;
	    gap[0] = gap[1] = 1.;
	    break;
	case 6:
	    size = .6;
	    dash[0] = 4.;
	    dash[1] = 2.;
	    gap[0] = gap[1] = 1.;
	    break;
	case 7:
	    size = .4;
	    dash[0] = dash[1] = 1.;
	    gap[0] = 2.;
	    gap[1] = 4.;
	    break;
	case 8:
	    size = .8;
	    dash[0] = dash[1] = 5.;
	    gap[0] = 2.;
	    gap[1] = 4.;
	    break;
	case 9:
	    size = .6;
	    dash[0] = dash[1] = 1.;
	    gap[0] = gap[1] = 1.;
	    break;
	case 0:
	    dash[0] = 0.;
	    gap[0] = 0.;
	    dash[1] = 0.;
	    gap[1] = 0.;
	    break;
	default:
	    dash[0] = 0.;
	    gap[0] = 0.;
	    dash[1] = 0.;
	    gap[1] = 0.;
	    break;
    }
    if (dash[0] + dash[1] + gap[0] + gap[1] != 0.) {
	/*If not default case then find the decimal part of dash->dashtype*/
	type = fmodf(type,1.);
	if (type == 0.) type = .4;
	size *= type / (.4 * (dash[0] + dash[1] + gap[0] + gap[1]));
	dash[0] *= size;
	dash[1] *= size;
	gap[0] *= size;
	gap[1] *= size;
    }
    vp_setdash (dash, gap, 2);
}

void vp_plot_init(int n2)
{
    int i;

    dashtype = sf_floatalloc(n2);
    if (!sf_getfloats ("dash",dashtype,n2)) {
	for (i = 0; i < n2; i++) 
	    dashtype[i] = 0.;
    }

    fat = sf_intalloc(n2);
    if (!sf_getints("plotfat",fat,n2)) {
	for (i = 0; i < n2; i++)
	    fat[i] = 0;
    }

    col = sf_intalloc(n2);
    if (!sf_getints("plotcol",col,n2)) {
	for (i = 0; i < n2; i++)
	    col[i] = 6 - (i % 6);
    }
}

void vp_plot_set (int i2)
{
    vp_fat (fat[i2]);
    vp_color (col[i2]);
    set_dash (dashtype[i2]);
}

void vp_plot_close (void)
{
    free (dashtype);
    free (fat);
    free (col);
}
