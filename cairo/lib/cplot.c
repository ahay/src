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

#include <cairo.h>

#include <gtk/gtk.h>
#include <gtkcairo.h>

#include <rsf.h>

#define INITIAL_WIDTH 737
#define INITIAL_HEIGHT 540

void cr_init(int* argc, char** argv[])
/*< Initialize >*/
{
    GtkWidget *win, *vbox, *frame, *gtkcairo;

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
			  INITIAL_WIDTH, INITIAL_HEIGHT);
/*    g_signal_connect (G_OBJECT (gtkcairo), "paint",
      G_CALLBACK (paint), slider); */

    gtk_container_add (GTK_CONTAINER (frame), gtkcairo);
    gtk_box_pack_start (GTK_BOX (vbox), frame, TRUE, TRUE, 0);
    
    gtk_container_add (GTK_CONTAINER (win), vbox);
    gtk_widget_show_all (vbox);

    gtk_widget_show (win);

    gtk_main ();
}


