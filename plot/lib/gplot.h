#ifndef _vp_gplot_h
#define _vp_gplot_h

#include <rsf.h>

void vp_coord_init (bool transp1, bool yreverse1);
void vp_plot_init(int n2);
void vp_title_init(sf_file file);
void vp_color_init (void);
void vp_minmax (float min1, float min2, float max1, float max2);
void vp_pad_init(bool pad, bool npad);
void vp_rotate (int n, float* x, float* y);
void vp_axis_init (const sf_file in);
void vp_vplot_init (void);
void vp_dash_fig (int type);

#endif
