/* Basic plot operations with cairo. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include <setjmp.h>

#include <cairo.h>

#include <gtk/gtk.h>
#include <gtkcairo.h>

#include <rsf.h>

#include "cplot.h"

#ifndef _cr_cplot_h

enum {
    CR_RPERIN=600, /* units per inch */
    CR_INITIAL_WIDTH=737,
    CR_INITIAL_HEIGHT=540
};
/*^*/

#define CR_SCALE 0.088
/* scaling fron units to screen */
#define CR_MAX 54.6           
/* absolute maximum x or y in inches */
/*^*/

#endif

static cairo_t *cr0;
static GtkWidget *win, *vbox, *frame, *gtkcairo;

static float fx=0.0, fy=0.0;     /* origin in inches */
static float ufx=0.0, ufy=0.0;   /* origin in user units */
static float xscl=1.0, yscl=1.0; /* scaling from user units to inches */
static float xold=0.0, yold=0.0; /* old pen position (in inches) */

static jmp_buf env; /* hacking gtk event loop */

static struct {
    int xmin, xmax, ymin, ymax;
} clip;

static void show (GtkWidget *widget, cairo_t *cr, gpointer data)
{
    gint width  = widget->allocation.width;
    gint height = widget->allocation.height;
    gint box_size = (width+height)/6;

    cairo_save (cr);
    cairo_default_matrix (cr);

    cairo_translate (cr, width/2, height/2);
    sf_warning("translating to %d %d",width/2, height/2);

    cairo_rectangle (cr, -box_size, -box_size, box_size, box_size);
    sf_warning("rectangle %d",box_size);

    cairo_set_rgb_color (cr, 1, 0, 0);
    cairo_fill (cr);
    
    cairo_restore (cr);
}

void cr_init(int* argc, char** argv[])
/*< Initialize >*/
{
    gtk_init (argc, argv);

    win = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title (GTK_WINDOW (win), sf_getprog());
    g_signal_connect (G_OBJECT (win), "delete-event",
                      G_CALLBACK (gtk_main_quit), NULL);

    vbox = gtk_vbox_new (FALSE, 6);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 12);

    frame = gtk_frame_new (NULL);
    gtk_frame_set_shadow_type (GTK_FRAME (frame), GTK_SHADOW_IN);

    gtkcairo = gtk_cairo_new ();
    gtk_widget_set_usize (GTK_WIDGET (gtkcairo), 
			  CR_INITIAL_WIDTH, CR_INITIAL_HEIGHT);

    g_signal_connect (G_OBJECT (gtkcairo), "paint",
		      G_CALLBACK (show), NULL);

    gtk_container_add (GTK_CONTAINER (frame), gtkcairo);
    gtk_box_pack_start (GTK_BOX (vbox), frame, TRUE, TRUE, 0);
    
    gtk_container_add (GTK_CONTAINER (win), vbox);
 
    cr0 = gtk_cairo_get_cairo(GTK_CAIRO(gtkcairo));
    cairo_set_rgb_color (cr0, 0, 1, 0);
}

void cr_main(void)
/*< go into main loop >*/
{
    if (0 == setjmp (env)) {
	gtk_widget_show_all (vbox);
	gtk_widget_show (win);
	gtk_main ();
    }
}

void cr_uorig (float x,float  y)
/*< set the origin in user coordinates >*/
{
    ufx = x;
    ufy = y;
}

static double snap(float x)
/* round to the specified resolution */
{
    x *= CR_RPERIN;
    return (CR_SCALE*((int) ((x < 0.0)? x-0.5 : x+0.5)));
}

void cr_clip (float xmin, float ymin, float xmax, float ymax)
/*< set rectangular clip >*/
{
    clip.xmin = snap(xmin);
    clip.ymin = snap(ymin);
    clip.xmax = snap(xmax);
    clip.ymax = snap(ymax);
}

void cr_uclip (float xmin, float ymin, float xmax, float ymax)
/*< set rectangular clip in user coordinates >*/
{
    xmin = fx + (xmin - ufx) * xscl;
    ymin = fy + (ymin - ufy) * yscl;
    xmax = fx + (xmax - ufx) * xscl;
    ymax = fy + (ymax - ufy) * yscl;
    cr_clip (xmin, ymin, xmax, ymax);
}

void cr_umove (float x,float  y)
/*< move to a point in user coordinates >*/
{
    x = fx + (x - ufx) * xscl;
    y = fy + (y - ufy) * yscl;
    cr_plot (x, y, false);
}

void cr_udraw (float x,float  y)
/*< line drawing step in user coordinates >*/
{
    x = fx + (x - ufx) * xscl;
    y = fy + (y - ufy) * yscl;
    cr_plot (x, y, true);
}

static void pout (float xp, float  yp, bool down)
{
    
    if      (xp >  CR_MAX) xp =  CR_MAX;
    else if (xp < -CR_MAX) xp = -CR_MAX;
    if      (yp >  CR_MAX) yp =  CR_MAX;
    else if (yp < -CR_MAX) yp = -CR_MAX;
    
    if (down) {
	cairo_line_to(cr0,snap(xp),snap(yp));
	sf_warning("drawing to %g,%g",snap(xp),snap(yp));
    } else {
	cairo_move_to(cr0,snap(xp),snap(yp));
	sf_warning("moving to %g,%g",snap(xp),snap(yp));
    }
}

void cr_plot (float x, float y, bool  down)
/*< line drawing >*/
{
    pout (x, y, down);	/* output a move or draw */
    xold = x;
    yold = y;		/* save old x and y */
}
