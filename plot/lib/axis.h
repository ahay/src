#ifndef _vp_axis_h
#define _vp_axis_h

float vp_opttic(float min, float max, float inch, float num0, float size);
void vp_simple_axis (float x1, float y1, 
		     float x2, float y2, 
		     float num1, float num2,
		     float onum, float dnum, 
		     float ltic, char* label, float size);

int vp_optimal_scale(float chars, float min, float max, 
		     float *onum, float *dnum);

#endif
