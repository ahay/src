/* Generate a surface plot. */
/*
  Copyright (C) 2010 The University of Texas at Austin
  
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

#include <math.h>
#include <string.h>

#include <rsf.h>
#include <rsfplot.h>

#include <plplot.h>

static void sf_plvpl_set_driver_dir (void);
static char* sf_plvpl_make_axis_title (char *label, char *unit);
static void* sf_plvpl_get_data_buffer (int n1, int n2, sf_datatype type);
static void sf_plvpl_read_data (int n1, int n2, void *buffer, sf_datatype type,
                                PLFLT **data, sf_file file);

int main (int argc, char *argv[]) {
    int n1, n2, n3;
    float o1, o2, o3, d1, d2, d3;
    int i, j;

    PLFLT *x, *y, **z;
    int nc;
    PLFLT *clevel = NULL;
    PLFLT zmin, zmax, step;

    char *color, *title = NULL;
    char *label1, *label2, *label3;
    char *unit1, *unit2, *unit3;
    char *axis1, *axis2, *axis3;
    bool wanttitle;

    float minval, maxval;
    float alt, az;
    int meshc;
    int opt = DRAW_LINEXY;
    bool sides, contour;

    void *buffer;
    sf_file in;
    sf_datatype type;

    /* Input data */    
    sf_init (argc, argv);
    in = sf_input ("in");
    type = sf_gettype (in);

    /* Dimensions */
    if (!sf_histint (in, "n1", &n1)) sf_error ("No n1= in input");
    if (!sf_histint (in, "n2", &n2)) sf_error ("No n2= in input");
    if (!sf_histint (in, "n3", &n3)) n3 = 1;
    if (!sf_histfloat (in, "o1", &o1)) o1 = 0.;
    if (!sf_histfloat (in, "o2", &o2)) o2 = 0.;
    if (!sf_histfloat (in, "o3", &o3)) o3 = 0.;
    if (!sf_histfloat (in, "d1", &d1)) d1 = 1.;
    if (!sf_histfloat (in, "d2", &d2)) d2 = 1.;
    if (!sf_histfloat (in, "d3", &d3)) d3 = 1.;

    /* Set up loadable driver directory */
    sf_plvpl_set_driver_dir ();

    /* Initialize plplot */
    plinit ();

    /* Axes */
    if (NULL == (label1 = sf_histstring (in, "label1"))) label1 = "";
    if (NULL == (label2 = sf_histstring (in, "label2"))) label2 = "";
    if (NULL == (label3 = sf_histstring (in, "label3"))) label3 = "";

    if (NULL == (unit1 = sf_histstring (in, "unit1"))) unit1 = "";
    if (NULL == (unit2 = sf_histstring (in, "unit2"))) unit2 = "";
    if (NULL == (unit3 = sf_histstring (in, "unit3"))) unit3 = "";

    axis1 = sf_plvpl_make_axis_title (label1, unit1);
    axis2 = sf_plvpl_make_axis_title (label2, unit2);
    axis3 = sf_plvpl_make_axis_title (label3, unit3);

    /* initialize color table */
    if (NULL == (color = sf_getstring ("color"))) color=NULL;
    /* color scheme (default is i) */
    if (color) {
        vp_rascoltab (VP_WHITE + 1, color);
        opt |= MAG_COLOR;
    }

    if (!sf_getbool ("wanttitle", &wanttitle)) wanttitle = true;
    /* if include title */
    if (wanttitle) {
        if (NULL == (title = sf_getstring ("title"))) title = NULL;
        /* title */
        if (NULL == title) wanttitle = false;
    }

    if (!sf_getfloat ("minval", &minval)) minval = SF_HUGE;
    /* minimum value for the vertical axis (default is data minimum) */
    if (!sf_getfloat ("maxval", &maxval)) maxval = -SF_HUGE;
    /* maximum value for the vertical axis (default is data maximum) */

    if (!sf_getfloat ("alt", &alt)) alt = 35.0;
    /* altitude [0;90] */
    if (!sf_getfloat ("az", &az)) az = 25.0;
    /* azimuth */
    if (alt > 90.0) alt = 0.0;
    if (alt < 0.0) alt = 0.0;

    if (!sf_getint ("meshc", &meshc)) meshc = VP_YELLOW;
    /* mesh color */
    if (meshc < VP_BLUE) meshc = VP_BLUE;
    if (meshc > VP_WHITE) meshc = VP_WHITE;

    if (!sf_getbool ("sides", &sides)) sides = false;
    /* draw sides */
    if (!sf_getbool ("contour", &contour)) contour = false;
    /* draw contour lines at the bottom */
    if (!sf_getint ("nc", &nc)) nc = 10;
    /* number of contour lines */
    if (nc < 1) nc = 1;
    if (contour) {
        opt |= BASE_CONT;
        clevel = (PLFLT*)calloc (nc, sizeof(PLFLT));
    }

    /* Set up x & y dimensions */
    x = (PLFLT*)calloc (n2, sizeof(PLFLT));
    y = (PLFLT*)calloc (n1, sizeof(PLFLT));
    plAlloc2dGrid (&z, n2, n1);
    for (i = 0; i < n1; i++)
        y[i] = o1 + i*d1;
    for (i = 0; i < n2; i++)
        x[i] = o2 + i*d2;

    buffer = sf_plvpl_get_data_buffer (n1, n2, type);

    /* Set font */
    plfontld (1);
    plfont (2);

    i = 0;
    while (i < n3) {
        sf_plvpl_read_data (n1, n2, buffer, type, z, in);

        plMinMax2dGrid (z, n2, n1, &zmax, &zmin);
        if (minval < zmin) zmin = minval;
        if (maxval > zmax) zmax = maxval;

        if (i > 0)
            vp_erase ();

        pladv (0);
        plcol0 (VP_WHITE);
        plvpor (0.0, 1.0, 0.0, 0.9);
        plwind (-0.85, 0.85, -0.6, 1.4);
        plw3d (1.0, 1.0, 1.2,
               o2, o2 + (n2 - 1)*d2,
               o1, o1 + (n1 - 1)*d1,
               zmin, zmax, alt, az);
        plbox3 ("bnstu", axis1, 0.0, 0,
	        "bnstu", axis2, 0.0, 0,
	        "bcdmnstuv", axis3, 0.0, 0);

        plcol0 (meshc);
        if (contour) {
            step = (zmax - zmin)/(nc + 1);
            for (j = 0; j < nc; j++)
                clevel[j] = zmin + step + step*j;
            plmeshc(x, y, z, n1, n1, opt, clevel, nc);
        } else {
            plot3d (x, y, z, n2, n1, opt, sides);
        }

        if (wanttitle) {
            plcol0 (VP_WHITE);
            plmtex ("t", 1.0, 0.5, 0.5, title);
        }

        i++;
    }

    /* Clean up */
    plend ();

    if (contour) free (clevel);
    free (x);
    free (y);
    plFree2dGrid (z, n2, n1);
    free (buffer);
    free (axis3);
    free (axis2);
    free (axis1);

    exit (0);
}

static void* sf_plvpl_get_data_buffer (int n1, int n2, sf_datatype type) {
    void *buffer;
    switch (type) {
        case SF_UCHAR:
            buffer = (void*)sf_ucharalloc (n1*n2);
            break;
        case SF_CHAR:
            buffer = (void*)sf_charalloc (n1*n2);
            break;
        case SF_INT:
            buffer = (void*)sf_intalloc (n1*n2);
            break;
        case SF_FLOAT:
            buffer = (void*)sf_floatalloc (n1*n2);
            break;
        case SF_COMPLEX:
            buffer = (void*)sf_complexalloc (n1*n2);
            break;
        case SF_SHORT:
            buffer = (void*)sf_shortalloc (n1*n2);
            break;
        case SF_DOUBLE:
            buffer = (void*)sf_alloc (n1*n2, sizeof (double));
            break;
    }
    return buffer;
}

static void sf_plvpl_read_data (int n1, int n2, void *buffer, sf_datatype type,
                                PLFLT **data, sf_file file) {
    int i, j;
    switch (type) {
        case SF_UCHAR: {
            unsigned char *buf = (unsigned char*)buffer;
            sf_ucharread (buf, n1*n2, file);
            for (i = 0; i < n2; i++) {
                for (j = 0; j < n1; j++) {
                    data[i][j] = buf[i*n1 + j];
                }
            }
            break;
        }
        case SF_CHAR: {
            char *buf = (char*)buffer;
            sf_charread (buf, n1*n2, file);
            for (i = 0; i < n2; i++) {
                for (j = 0; j < n1; j++) {
                    data[i][j] = buf[i*n1 + j];
                }
            }
            break;
        }
        case SF_INT: {
            int *buf = (int*)buffer;
            sf_intread (buf, n1*n2, file);
            for (i = 0; i < n2; i++) {
                for (j = 0; j < n1; j++) {
                    data[i][j] = buf[i*n1 + j];
                }
            }
            break;
        }
        case SF_FLOAT: {
            float *buf = (float*)buffer;
            sf_floatread (buf, n1*n2, file);
            for (i = 0; i < n2; i++) {
                for (j = 0; j < n1; j++) {
                    data[i][j] = buf[i*n1 + j];
                }
            }
            break;
        }
        case SF_COMPLEX: {
            sf_complex *buf = (sf_complex*)buffer;
            sf_complexread (buf, n1*n2, file);
            for (i = 0; i < n2; i++) {
                for (j = 0; j < n1; j++) {
                    data[i][j] = hypot (creal (buf[i*n1 + j]), cimag (buf[i*n1 + j]));
                }
            }
            break;
        }
        case SF_SHORT: {
            short *buf = (short*)buffer;
            sf_shortread (buf, n1*n2, file);
            for (i = 0; i < n2; i++) {
                for (j = 0; j < n1; j++) {
                    data[i][j] = buf[i*n1 + j];
                }
            }
            break;
        }
        case SF_DOUBLE: {
            double *buf = (double*)buffer;
            sf_ucharread ((unsigned char*)buf, n1*n2*sizeof(double), file);
            for (i = 0; i < n2; i++) {
                for (j = 0; j < n1; j++) {
                    data[i][j] = buf[i*n1 + j];
                }
            }
            break;
        }
    }
}

#define BUF_SIZ 2048
static void sf_plvpl_set_driver_dir (void) {
    char buf[BUF_SIZ];
    char *rsfroot = getenv ("RSFROOT");

    if (NULL == rsfroot)
        sf_error ("Need RSFROOT environment variable to run");

    snprintf (buf, BUF_SIZ, "%s/lib", rsfroot);
    setenv ("PLPLOT_DRV_DIR", buf, 1);
}

static char* sf_plvpl_make_axis_title (char *label, char *unit) {
    int i, j;
    char *title;

    i = label ? strlen (label) : 0;
    j = unit ? strlen (unit) : 0;
    title = (char*)calloc (i + j + 4, sizeof(char));
    if (i && j)
        sprintf (title, "%s (%s)", label, unit);
    else if (i)
        sprintf (title, "%s", label);
    else if (j)
        sprintf (title, "(%s)", unit);

    return title;
}

