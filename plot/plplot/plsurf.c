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
#include <stdio.h>
#include <stdlib.h>

#include <rsf.h>
#include <rsfplot.h>

#include "plvpl.h"

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
    char *axis1 = "", *axis2 = "", *axis3 = "";
    char *where;
    bool wanttitle, wantaxis;
    bool wantaxis1 = false, wantaxis2 = false, wantaxis3 = false;
    float labelsz, titlesz;
    int labelfat, titlefat, font;

    float minval, maxval;
    float alt, az;
    int meshc;
    int opt = 0;
    bool mesh, sides, bcontour, scontour, faceted;

    void *buffer;
    sf_file in;
    sf_datatype type;

    /* Input data */    
    sf_init (argc, argv);
    in = sf_input ("in");
    type = sf_gettype (in);

    vp_init();

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

    /* Axes; allow user redefine label and units over values from the history */
    if (!sf_getbool ("wantaxis", &wantaxis)) wantaxis = true;
    /* if generate axes with ticks and labels */
    if (wantaxis) {
        if (!sf_getbool ("wantaxis1", &wantaxis1)) wantaxis1 = true;
        if (false == wantaxis1 ||
            (NULL == (label1 = sf_getstring ("label1")) &&
             NULL == (label1 = sf_histstring (in, "label1")))) label1 = "";
        if (!sf_getbool ("wantaxis2", &wantaxis2)) wantaxis2 = true;
        if (false == wantaxis2 ||
            (NULL == (label2 = sf_getstring ("label2")) &&
             NULL == (label2 = sf_histstring (in, "label2")))) label2 = "";
        if (!sf_getbool ("wantaxis3", &wantaxis3)) wantaxis3 = true;
        if (false == wantaxis3 ||
            (NULL == (label3 = sf_getstring ("label3")) &&
             NULL == (label3 = sf_histstring (in, "label3")))) label3 = "";
        if (false == wantaxis1 ||
            (NULL == (unit1 = sf_getstring ("unit1")) &&
             NULL == (unit1 = sf_histstring (in, "unit1")))) unit1 = "";
        if (false == wantaxis2 ||
            (NULL == (unit2 = sf_getstring ("unit2")) &&
             NULL == (unit2 = sf_histstring (in, "unit2")))) unit2 = "";
        if (false == wantaxis3 ||
            (NULL == (unit3 = sf_getstring ("unit3")) &&
             NULL == (unit3 = sf_histstring (in, "unit3")))) unit3 = "";
        axis1 = sf_plvpl_make_axis_title (label1, unit1);
        axis2 = sf_plvpl_make_axis_title (label2, unit2);
        axis3 = sf_plvpl_make_axis_title (label3, unit3);
    }

    if (!sf_getbool ("wanttitle", &wanttitle)) wanttitle = true;
    /* if include title */
    if (wanttitle) {
        if (NULL == (title = sf_getstring ("title")) &&
            NULL == (title = sf_histstring (in, "title"))) title = NULL;
        if (NULL == title) wanttitle = false;
    }

    if (NULL == (where = sf_getstring ("wheretitle")) || 
	0 == strlen (where)) where = "top";
    /* where to put title (top,bottom) */

    if (!sf_getfloat ("labelsz", &labelsz)) labelsz = 8.0;
    /* label font size */
    if (!sf_getfloat ("titlesz", &titlesz)) titlesz = 10.0;
    /* title font size */
    if (!sf_getint ("labelfat", &labelfat)) labelfat = 1;
    /* label fatness */
    if (!sf_getint ("titlefat", &titlefat)) titlefat = 1;
    /* title fatness */

    if (!sf_getint ("font", &font)) font = 2;
    /* font */
    if (font < 1) font = 1;
    if (font > 4) font = 4;

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

    if (!sf_getbool ("mesh", &mesh)) mesh = true;
    /* what to draw: true - mesh, false - shaded surface */
    if (mesh) {
        opt |= DRAW_LINEXY;
        opt |= MESH;
    }
    if (!sf_getint ("meshc", &meshc)) meshc = VP_YELLOW;
    /* mesh color or surface contour color */
    if (meshc < VP_BLUE) meshc = VP_BLUE;
    if (meshc > VP_WHITE) meshc = VP_WHITE;

    if (!sf_getbool ("sides", &sides)) sides = false;
    /* draw sides */
    if (sides)
        opt |= DRAW_SIDES;
    if (!sf_getbool ("bcontour", &bcontour)) bcontour = false;
    /* draw contour lines at the bottom */
    if (bcontour)
        opt |= BASE_CONT;
    if (!sf_getbool ("scontour", &scontour)) scontour = false;
    /* draw contour lines on the surface (surface mode only) */
    if (scontour)
        opt |= SURF_CONT;
    if (!sf_getint ("nc", &nc)) nc = 10;
    /* number of contour lines */
    if (nc < 1) nc = 1;
    if (bcontour || scontour)
        clevel = (PLFLT*)calloc (nc, sizeof(PLFLT));
    else 
        nc = 0;
    if (!sf_getbool ("faceted", &faceted)) faceted = false;
    /* each cell is faceted on the surface (surface mode only) */
    if (faceted)
        opt |= FACETED;

    /* initialize color table */
    if (NULL == (color = sf_getstring ("color"))) color=NULL;
    /* color scheme (default is i) */
    if (color) {
        vp_rascoltab (VP_WHITE + 1, color);
        opt |= MAG_COLOR;
    } else if (false == mesh) { /* Set default b/w color palette for surface rendering */
        vp_rascoltab (VP_WHITE + 1, "b");
        opt |= MAG_COLOR;
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
    plfontld (1); /* Extended set */
    plfont (font);

    i = 0;
    while (i < n3) {
        sf_plvpl_read_data (n1, n2, buffer, type, z, in);

        plMinMax2dGrid (z, n2, n1, &zmax, &zmin);
        if (minval < zmin) zmin = minval;
        if (maxval > zmax) zmax = maxval;

        if (scontour || bcontour) {
            step = (zmax - zmin)/(nc + 1);
            for (j = 0; j < nc; j++)
                clevel[j] = zmin + step + step*j;
        }

        if (i > 0)
            vp_erase ();
        pladv (0);
        plcol0 (VP_WHITE);
        plvpor (0.0, 1.0,
                where[0] != 't' && wanttitle ? 0.1 : 0.0,
                where[0] != 't' || !wanttitle ? 1.0 : 0.9);
        plwind (-0.9 - fdimf (labelsz, 8.0)*0.025, 0.9 + fdimf (labelsz, 8.0)*0.025,
                -0.6 - fdimf (labelsz, 8.0)*0.02, 1.45);
        plw3d (1.0, 1.0, 1.2,
               o2, o2 + (n2 - 1)*d2,
               o1, o1 + (n1 - 1)*d1,
               zmin, zmax, alt, az);
        plschr (0, labelsz/10.0);
        vp_fat (labelfat);
        plbox3 (wantaxis && wantaxis2 ? "bnstu" : "b", axis2, 0.0, 0,
                wantaxis && wantaxis1 ? "bnstu" : "b", axis1, 0.0, 0,
                wantaxis && wantaxis3 ? "bcdmnstuv" : "bcduv", axis3, 0.0, 0);
        vp_fat (1);
        plschr (0, 1.0);

        plcol0 (meshc);
        if (opt & MESH)
            plot3dc (x, y, z, n2, n1, opt, clevel, nc);
        else
            plsurf3d (x, y, z, n2, n1, opt, clevel, nc);

        if (wanttitle) {
            plcol0 (VP_WHITE);
            plschr (0, titlesz/10.0);
            vp_fat (titlefat);
            plmtex (where[0] != 'b' ? "t" : "b", 1.0, 0.5, 0.5, title);
            vp_fat (1);
            plschr (0, 1.0);
        }

        i++;
    }

    /* Clean up */
    plend ();

    if (clevel) free (clevel);
    free (x);
    free (y);
    plFree2dGrid (z, n2, n1);
    free (buffer);
    if (wantaxis) {
        free (axis3);
        free (axis2);
        free (axis1);
    }

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
                    data[i][j] = hypot (creal (buf[i*n1 + j]), 
					cimag (buf[i*n1 + j]));
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

    sprintf (buf, "%s/lib", rsfroot);
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

