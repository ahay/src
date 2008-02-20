/* 
 * Additional service for the trackball: window coordinates->rotation, etc
 *
 * Copyright (C) 2007 The University of Texas at Austin
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more av.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * Author:  Vladimir Bashkardin  <vovizmus@users.sourceforge.net>
 */

#include "trackserv.h"

#if defined(__APPLE__)&& defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define ZOOM_DEF 1.0f
#define ZOOM_MIN 0.25f
#define ZOOM_MAX 4.0f
#define ZOOM_STEP 0.125f

#if !defined(GLUT_WHEEL_UP)
#define GLUT_WHEEL_UP   3
#define GLUT_WHEEL_DOWN 4
#endif


static void rsf_trackball_motion_rotate_event (RSFTrackballData *data, int x,int y) {

    if (NULL == data)
        return;

    trackball (data->dquat,
               (2.0 * data->begin_x - data->screen_width) / data->screen_width * data->zoom,
               (data->screen_height - 2.0 * data->begin_y) / data->screen_height * data->zoom,
               (2.0 * x - data->screen_width) / data->screen_width * data->zoom,
               (data->screen_height - 2.0 * y) / data->screen_height * data->zoom);

    data->dx = x - data->begin_x;
    data->dy = y - data->begin_y;

    add_quats (data->dquat,
               data->quat,
               data->quat);
    build_rotmatrix (data->m,
                     data->quat);

    data->begin_x = x;
    data->begin_y = y;
    data->redraw = 1;
}

static void rsf_trackball_get_zero_z_position (RSFTrackballData *data, int x,int y,
                                               GLdouble *wx, GLdouble *wy) {
    GLint    viewport_matrix[4];
    GLdouble projection_matrix[16];
    GLdouble modelview_matrix[16];
    GLdouble vx, vy, vz, wz;

    if (NULL == wx || NULL == wy)
        return;

    glPushMatrix ();
    glLoadIdentity ();
    glTranslatef (0, 0, data->translate);

    glGetIntegerv (GL_VIEWPORT, viewport_matrix);
    glGetDoublev (GL_PROJECTION_MATRIX, projection_matrix);
    glGetDoublev (GL_MODELVIEW_MATRIX, modelview_matrix);

    vx = x;
    vy = data->screen_height - y;
    vz = 0;
    gluUnProject (vx, vy, vz, modelview_matrix, projection_matrix, viewport_matrix, wx, wy, &wz);

    *wx *= 1.5;
    *wy *= 1.5;

    glPopMatrix ();
}

static void rsf_trackball_motion_move_event (RSFTrackballData *data, int x,int y) {
    GLdouble wx, wy;

    if (NULL == data)
        return;

    rsf_trackball_get_zero_z_position (data, x, y, &wx, &wy);

    data->shift_wx += wx - data->last_wx;
    data->shift_wy += wy - data->last_wy;

    data->last_wx = wx;
    data->last_wy = wy;

    data->redraw = 1;
}

void rsf_trackball_mouse_event (RSFTrackballData *data, int button, int state, int x, int y) {
    if (GLUT_LEFT_BUTTON == button) {
        if (GLUT_DOWN == state) {
            data->dquat[0] = 0.0;
            data->dquat[1] = 0.0;
            data->dquat[2] = 0.0;
            data->dquat[3] = 1.0;
            data->begin_x = x;
            data->begin_y = y;
            data->motion_callback = rsf_trackball_motion_rotate_event;
        } else {
            data->motion_callback = NULL;
            data->dx = 0.0;
            data->dy = 0.0;
        }
    } else if (GLUT_WHEEL_UP == button) {
        if (GLUT_DOWN == state) {
            data->zoom += ZOOM_STEP;
            if (data->zoom > ZOOM_MAX)
                data->zoom = ZOOM_MAX;
            data->reshape = 1;
            data->redraw = 1;
        }
    } else if (GLUT_WHEEL_DOWN == button) {
        if (GLUT_DOWN == state) {
            data->zoom -= ZOOM_STEP;
            if (data->zoom < ZOOM_MIN)
                data->zoom = ZOOM_MIN;
            data->reshape = 1;
            data->redraw = 1;
        }
    } else if (button == GLUT_MIDDLE_BUTTON) {
        if (GLUT_DOWN == state) {
            rsf_trackball_get_zero_z_position (data, x, y, &data->last_wx, &data->last_wy);
            data->motion_callback = rsf_trackball_motion_move_event;
        } else {
            data->motion_callback = NULL;
            data->last_wx = 0;
            data->last_wy = 0;
        }
    }
}

void rsf_trackball_init (RSFTrackballData *data) {

    if (NULL == data)
        return;

    data->begin_x = 0;
    data->begin_y = 0;
    data->zoom = ZOOM_DEF;
    data->near_plane = 5.0;
    data->far_plane = 60.0;
    data->translate = -8.5;
    data->dx = 0;
    data->dy = 0;
    data->quat[0] = 0;
    data->quat[1] = 0;
    data->quat[2] = 0;
    data->quat[3] = 1;
    data->dquat[0] = 0;
    data->dquat[1] = 0;
    data->dquat[2] = 0;
    data->dquat[3] = 1;
    data->motion_callback = NULL;
    data->shift_wx = 0;
    data->last_wx = 0;
    data->shift_wy = 0;
    data->last_wy = 0;
    data->redraw = 0;
    data->reshape = 0;
    trackball (data->quat, 0.0, 0.0, 0.0, 0.0);
    build_rotmatrix (data->m, data->quat);
}
