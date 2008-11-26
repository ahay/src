/* -*- C -*-  (not really, but good for syntax highlighting) */

%module c_vplot

%{

#include <rsfplot.h>

%}

void vp_init(void);
void vp_uorig (float x,float  y);
void vp_uclip (float xmin, float ymin, float xmax, float ymax);
void vp_umove (float x,float  y);
void vp_udraw (float x,float  y);
void vp_move (float x,float  y);
void vp_draw (float x,float  y);
void vp_fat (int f);
void vp_color (int col);
void vp_penup (void);
void vp_pendn (float x, float y);
void vp_text (float x, float y    /* coordinate of the reference point */, 
	      int size            /* height of character */, 
	      int orient          /* text drawing direction ( in degrees counter-clockwise
				     from horizontal, right-facing) */, 
	      const char *string /* test */);
