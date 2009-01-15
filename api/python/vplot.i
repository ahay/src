/* -*- C -*-  (not really, but good for syntax highlighting) */

%module c_vplot

%include carrays.i
%array_functions(float,floatp);

%{

#include <rsfplot.h>

%}

void vp_init(void);
void vp_orig (float x,float  y);
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
void vp_utext (float x, float y    /* coordinate of the reference point */, 
	       int size            /* height of character */, 
	       int orient          /* text drawing direction ( in degrees counter-clockwise
				      from horizontal, right-facing) */, 
	       const char *string /* test */);
void vp_scale (float xscale, float  yscale);
void vp_uarrow (float x1, float y1, float x, float y, float r);
void vp_tjust (int xjust1, int yjust1);
void vp_clip (float xmin, float ymin, float xmax, float ymax);
void vp_dash (float dash1, float gap1, float dash2, float gap2);
void vp_upline (const float *xp /* [np] */, 
		const float *yp /* [np] */, 
		int np          /* number of points */);
void vp_upendn (float x, float y);
