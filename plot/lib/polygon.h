#ifndef _vp_polygon_h
#define _vp_polygon_h

#include "device.h"

void vp_polyfix (int x, int y, bool *first, bool *allgone);

void vp_ymaxclip (int xin, int yin, int *first, vp_device dev);
void vp_xmaxclip (int xin, int yin, int *first, vp_device dev);
void vp_yminclip (int xin, int yin, int *first, vp_device dev);
void vp_xminclip (int xin, int yin, int *first, vp_device dev);

#endif
