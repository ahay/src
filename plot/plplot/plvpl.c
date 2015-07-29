/* vplot loadable module for PLPLOT. */
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
#include <stdio.h>
#include <unistd.h>
#include <rsfplot.h>

#include <plplot.h>
#include <plplotP.h>
#include <plstrm.h>
/*^*/

#include "plvpl.h"

/* Device info */

PLDLLIMPEXP_DRIVER const char* plD_DEVICE_INFO_plvpl = 
  "plvpl:VPLOT:0:plvpl:0:plvpl\n";

void plD_dispatch_init_plvpl (PLDispatchTable *pdt) 
/*< initialize >*/
{
    pdt->pl_MenuStr = strdup ("VPLOT");
    pdt->pl_DevName = strdup ("plvpl");
    pdt->pl_type = plDevType_FileOriented;
    pdt->pl_seq = 0;
    pdt->pl_init     = (plD_init_fp)     plD_init_plvpl;
    pdt->pl_line     = (plD_line_fp)     plD_line_plvpl;
    pdt->pl_polyline = (plD_polyline_fp) plD_polyline_plvpl;
    pdt->pl_eop      = (plD_eop_fp)      plD_eop_plvpl;
    pdt->pl_bop      = (plD_bop_fp)      plD_bop_plvpl;
    pdt->pl_tidy     = (plD_tidy_fp)     plD_tidy_plvpl;
    pdt->pl_state    = (plD_state_fp)    plD_state_plvpl;
    pdt->pl_esc      = (plD_esc_fp)      plD_esc_plvpl;
}

void plD_init_plvpl (PLStream *pls) 
/*< Initialize the device. >*/
{
    pls->color = 1;
    pls->colorset = 0;
    pls->termin = 0;
    pls->width = 1; /* pen width */
    pls->dev_text = 0; /* don't want to draw text */
    pls->dev_fill0 = 1; /* can do solid fills */
    pls->dev_fill1 = 0; /* don't do pattern pattern fills */
    pls->dev_fastimg = 0; /* give up raster */
    pls->ncol0 = VP_WHITE + 1;
    pls->ncol1 = 256;
    /* Dummy colormaps - VPLOT pens have their own defaults */
    pls->cmap0 = (PLColor *)calloc(1, pls->ncol0*sizeof(PLColor));
    pls->cmap1 = (PLColor *)calloc(1, pls->ncol1*sizeof(PLColor));
    pls->verbose = 0;
    pls->debug = 0;

    pls->graphx = GRAPHICS_MODE;

    if (isatty (fileno (stdout))) {
        fprintf (stderr, "You don't want to dump binary to terminal.\n");
        exit (-1);
    }

    pls->OutFile = stdout;
    pls->output_type = 1;
    vp_filep (pls->OutFile);
    vp_style (VP_STANDARD);
    vp_orig (0, 0);

    pls->xlength = (PLINT) (VP_STANDARD_HEIGHT/(float)VP_SCREEN_RATIO);
    pls->ylength = (PLINT) VP_STANDARD_HEIGHT;
    pls->xdpi = RPERIN;
    pls->ydpi = RPERIN;
    plP_setpxl ((PLFLT)RPERIN/25.4, (PLFLT)RPERIN/25.4);
    plP_setphy (0, (PLINT)(pls->xlength*pls->xdpi),
                0, (PLINT)(pls->ylength*pls->ydpi));
}

void plD_line_plvpl (PLStream *pls, short x1a, short y1a, short x2a, short y2a)
/*< Draw a line in the current color from (x1,y1) to (x2,y2). >*/
{
    float xp[2], yp[2];

    xp[0] = x1a/(PLFLT)RPERIN; xp[1] = x2a/(PLFLT)RPERIN;
    yp[0] = y1a/(PLFLT)RPERIN; yp[1] = y2a/(PLFLT)RPERIN;
    vp_pline (xp, yp, 2);
}

void plD_polyline_plvpl (PLStream *pls, short *xa, short *ya, PLINT npts) 
/*< Draw a polyline in the current color. >*/
{
    PLINT i;

    for (i = 0; i < npts - 1; i++)
        plD_line_plvpl (pls, xa[i], ya[i], xa[i + 1], ya[i + 1]);
}

void plD_eop_plvpl (PLStream *pls) 
/*< End of page. >*/
{
/*
    vp_purge ();
*/
}

void plD_bop_plvpl (PLStream *pls) 
/*< Set up for the next page. >*/
{
/*
    vp_erase ();
*/
}

void plD_tidy_plvpl (PLStream *pls) 
/*< Close output or otherwise clean up. >*/
{
/*
    vp_purge ();
*/
}

void plD_state_plvpl (PLStream *pls, PLINT op) 
/*< Handle change in PLStream state (color, pen width, fill attribute, etc). >*/
{
    int i;

    switch (op) {
        case PLSTATE_WIDTH:
            vp_fat (pls->width);
            break;
        case PLSTATE_COLOR0:
            vp_color (pls->icol0);
            break;
        case PLSTATE_COLOR1:
            vp_color (256 + pls->icol1);
            break;
        case PLSTATE_CMAP0:
            for (i = 0; i < pls->ncol0; i++)
                vp_coltab (i, pls->cmap0[i].r/(float)MAX_GUN,
                              pls->cmap0[i].g/(float)MAX_GUN,
                              pls->cmap0[i].b/(float)MAX_GUN);
            break;
        case PLSTATE_CMAP1:
            for (i = 0; i < pls->ncol1; i++)
                vp_coltab (256 + i,
                           pls->cmap1[i].r/(float)MAX_GUN,
                           pls->cmap1[i].g/(float)MAX_GUN,
                           pls->cmap1[i].b/(float)MAX_GUN);
            break;
    }
}

static void fill_polygon (PLStream *pls);

void plD_esc_plvpl (PLStream *pls, PLINT op, void *ptr) 
/*< Escape function. >*/
{
    switch (op) {
        case PLESC_FILL:
            fill_polygon (pls);
            break;
/*
        case PLESC_IMAGE:
            draw_image (pls);
            break;
*/
    }
}

static void fill_polygon (PLStream *pls) {
    int i;
    float *xp, *yp; 

    xp = sf_floatalloc (pls->dev_npts);
    yp = sf_floatalloc (pls->dev_npts);

    for (i = 0; i < pls->dev_npts; i++) {
        xp[i] = pls->dev_x[i]/(PLFLT)RPERIN;
        yp[i] = pls->dev_y[i]/(PLFLT)RPERIN;
    }

    vp_fill (xp, yp, pls->dev_npts);  

    free (yp);
    free (xp);
}

