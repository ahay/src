#ifndef _vp_gplot_h
#define _vp_gplot_h

#include <rsf.h>

void vp_simpleaxis (float x1, float y1, 
		    float x2, float y2, 
		    float num1, float num2,
		    float dnum, float ltic, char* label, float labelsz);

void vp_coord_init (bool transp1, bool yreverse1);
void vp_plot_init(int n2);
void vp_title_init(sf_file file);
void vp_color_init (void);
void vp_minmax (float min1, float min2, float max1, float max2);
void vp_pad_init(bool pad, bool npad);

#endif
