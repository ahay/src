/*
 * Usage: Jontest color= | Tube
 *
 * color same as Ta2vplot, or if longer than 2 characters is a colfile
 * name. "IC" is your favorite color table.
 *
 * raster vplot example
 */

#include <rsf.h>
#include <rsfplot.h>

int main (int argc, char* argv[])
{
    unsigned char   **array;
    int             ii, jj;
    char            *color;
    int             offset, xpix, ypix, bit, blast;
    float           xll, yll, xur, yur, ppi;
    
    array = sf_ucharalloc2(100,100);

    sf_init(argc,argv);
    vp_init();

    if (NULL == (color = sf_getstring("color"))) color="I";

    for (ii = 0; ii < 100; ii++)
	for (jj = 0; jj < 100; jj++)
	    array[ii][jj] = (ii + 2 * jj) % 256;
    
/*
 * Reserve the 8 standard vplot colors, and fill up the rest with
 * the desired color table, suitably scrambled once on the bottom
 * and not scrambled once on the top.
 *
 * If you want to plot on an Envision, which has 8 colors, you'll
 * need those 8 for the raster so you should change this to 0... but
 * you change the background when you grab color 0, and you
 * lose the proper colors on the two diagonal lines.
 */
    vp_rascoltab (8, color);

/* This is to tell vp_raster to use the upper table */
    offset = 256;

    xpix = 100;
    ypix = 100;
    bit = 0;
    xll = 1.;
    yll = 1.;
    xur = 5.;
    yur = 5.;
    ppi = 0;
    blast = 0;

    vp_raster (array, bit, offset, xpix, ypix, xll, yll, xur, yur, 1); 

/*
 * To show we can put colors on top of raster.
 */
    vp_color (VP_CYAN);
    vp_move (xll, yll);
    vp_draw (xur, yur);
    vp_color (VP_GREEN);
    vp_move (xll, yur);
    vp_draw (xur, yll);

    exit(0);
}
