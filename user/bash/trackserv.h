/* 
  Additional service for the trackball: window coordinates->rotation, etc
*/
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#ifndef __TRACK_BALL_SERVICE_H__
#define __TRACK_BALL_SERVICE_H__

#if defined(__APPLE__)&& defined(__MACH__)
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include "trackball.h"

typedef struct _RSFTrackballData RSFTrackballData;

struct _RSFTrackballData {
    float     begin_x, begin_y;  /* position of mouse */
    float     zoom;
    float     dx, dy;
    float     quat[4]; /* orientation of object */
    float     dquat[4];
    GLdouble  shift_wx, shift_wy;
    GLdouble  last_wx, last_wy;
    GLfloat   near_plane, far_plane;
    GLfloat   translate;
    GLfloat   world_width, world_height;
    GLint     screen_width, screen_height;
    GLfloat   m[4][4]; /* Rotation matrix */
    void(*motion_callback)(RSFTrackballData*, int, int); 
    int       redraw; /* if 1 - needs redraw after motion_callback); */
    int       reshape; /* if 1 - needs reshape after mouse_callback); */
};

void rsf_trackball_mouse_event (RSFTrackballData *data, int button, int state, int x, int y);
void rsf_trackball_init (RSFTrackballData *data);

#endif /* __TRACK_BALL_SERVICE_H__ */
