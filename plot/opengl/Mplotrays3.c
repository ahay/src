/* Plot rays in 3D with OpenGL.

Takes: < rays.bin

rays.bin should be in the special binary format for variable length rays

Run "sfdoc stdplot" for more parameters.
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


#include <rsf.h>
#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#ifdef FREEGLUT
#ifdef __APPLE__
#include <GLUT/freeglut_ext.h>
#else
#include <GL/freeglut_ext.h>
#endif
#endif

#include "trackserv.h"
#include "gl2ps.h"

RSFTrackballData trackball_data;

GLuint rays_list_base = 0;
GLuint rays_list_num = 0;

GLint num_aux_buffers = 0;
GLuint highlighted_ray = 0;

float scale_x = 1.0, scale_y = 1.0, scale_z = 1.0;

int is_rendering_eps = 0;
int has_drawn_cross = 0;
GLdouble cross_x = -10.0;
GLdouble cross_y = -10.0;
GLdouble cross_z = -10.0;

GLint    viewport_matrix[4];
GLdouble projection_matrix[16];
GLdouble modelview_matrix[16];

void plot_rays_draw_bounding_box () {
    glPushMatrix ();
    glScalef (scale_x, scale_y, scale_z);

    /* Bounding box */
    glBegin (GL_QUAD_STRIP);

    glVertex3f (-1.0, -1.0, 1.0);
    glVertex3f (-1.0, 1.0, 1.0);

    glVertex3f (1.0, -1.0, 1.0);
    glVertex3f (1.0, 1.0, 1.0);

    glVertex3f (1.0, -1.0, -1.0);
    glVertex3f (1.0, 1.0, -1.0);

    glVertex3f (-1.0, -1.0, -1.0);
    glVertex3f (-1.0, 1.0, -1.0);

    glVertex3f (-1.0, -1.0, 1.0);
    glVertex3f (-1.0, 1.0, 1.0);

    glEnd ();

    GLdouble vz;
    GLdouble wx, wy, wz;
    GLfloat xy_start[3] = { 0, -1, 0 };
    GLfloat x_end[3] = { 1, -1, 0 };
    GLfloat y_end[3] = { 0, -1, 1 };

    glGetDoublev (GL_MODELVIEW_MATRIX, modelview_matrix);
    gluProject (0, -1, 0, modelview_matrix, projection_matrix, viewport_matrix, &wx, &wy, &wz);
    vz = wz;
printf ("%f ", wz);
    gluProject (1, -1, 0, modelview_matrix, projection_matrix, viewport_matrix, &wx, &wy, &wz);
    if (wz < vz) {
        vz = wz;
        xy_start[0] = 1;
        xy_start[2] = 0;
        x_end[0] = 0;
        x_end[2] = 0;
        y_end[0] = 1;
        y_end[2] = 1;
    }
printf ("%f ", wz);
    gluProject (1, -1, 1, modelview_matrix, projection_matrix, viewport_matrix, &wx, &wy, &wz);
    if (wz < vz) {
        vz = wz;
        xy_start[0] = 1;
        xy_start[2] = 1;
        x_end[0] = 1;
        x_end[2] = 0;
        y_end[0] = 0;
        y_end[2] = 1;
    }
printf ("%f ", wz);
    gluProject (0, -1, 1, modelview_matrix, projection_matrix, viewport_matrix, &wx, &wy, &wz);
    if (wz < vz) {
        vz = wz;
        xy_start[0] = 0;
        xy_start[2] = 1;
        x_end[0] = 1;
        x_end[2] = 1;
        y_end[0] = 1;
        y_end[2] = 1;
    }
printf ("%f\n", wz);

    glPopMatrix ();
}

void plot_rays_redraw (void) {
    int i;

    glLoadIdentity ();
    glTranslatef (trackball_data.shift_wx, trackball_data.shift_wy,
                  trackball_data.translate);

    glMultMatrixf (&trackball_data.m[0][0]);

    /* Render rays in the auxiliary buffer for further mouse selection */
    if (num_aux_buffers > 0 && 0 == is_rendering_eps) {
        glDrawBuffer (GL_AUX0);
        glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLineWidth (3);
        for (i = 0; i < rays_list_num; i++)
            glCallList ((rays_list_base + i) * 2 + 1);
        glDrawBuffer (GL_BACK);
        glLineWidth (1);
    }

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (is_rendering_eps)
        glColor3f (0.0, 0.0, 0.0);
    else
        glColor3f (1.0, 1.0, 1.0);


    plot_rays_draw_bounding_box ();

    glPushMatrix ();
    glScalef (scale_x, scale_y, scale_z);
/*
glPushMatrix();
glLoadIdentity ();
	glTranslatef(-1,0,-8.5);
	glScalef(0.005,0.005,0.005);
//	glLineWidth(4.0);
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'A'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'B'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'C'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'D'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'E'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'F'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'G'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'H'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'I'); 
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'J'); 
	glPopMatrix ();*/
    if (is_rendering_eps)
        glColor3f (0.0, 0.0, 0.0);
    else
        glColor3f (1.0, 1.0, 1.0);

    /* Rays */
    for (i = 0; i < rays_list_num; i++) {
        if ((i + 1) == highlighted_ray) {
            glLineWidth (2);
            glColor3f (1, 0, 0);
        }
        glCallList ((rays_list_base + i) * 2);
        if ((i + 1) == highlighted_ray) {
            if (is_rendering_eps)
                glColor3f (0.0, 0.0, 0.0);
            else
                glColor3f (1.0, 1.0, 1.0);
            glLineWidth (1);
        }
    }

    glPopMatrix ();

    glFlush ();
    glutSwapBuffers ();

    has_drawn_cross = 0;
}

void plot_rays_reshape (int width, int height) {
    GLdouble wh, world_width, world_height;

    if (width <= height) {
        wh = (GLfloat)height / (GLfloat)width;
        world_width = 2.0;
        world_height = wh * 2.0;
    } else {
        wh = (GLfloat)width / (GLfloat)height;
        world_width = wh * 2.0;
        world_height = 2.0;
    }

    glViewport (0, 0, width, height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glFrustum (-world_width / 2.0 * trackball_data.zoom,
                world_width / 2.0 * trackball_data.zoom,
               -world_height / 2.0 * trackball_data.zoom,
                world_height / 2.0 * trackball_data.zoom,
               trackball_data.near_plane, trackball_data.far_plane);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    glGetIntegerv (GL_VIEWPORT, viewport_matrix);
    glGetDoublev (GL_PROJECTION_MATRIX, projection_matrix);

    trackball_data.world_width = world_width;
    trackball_data.world_height = world_height;
    trackball_data.screen_width = width;
    trackball_data.screen_height = height;
}

void plot_rays_mouse (int button, int state, int x, int y) {
    GLdouble vx, vy;
    GLubyte col[3] = { 0, 0, 0 };

    if (GLUT_LEFT_BUTTON == button && GLUT_DOWN == state && GLUT_ACTIVE_CTRL == glutGetModifiers () &&
        num_aux_buffers > 0) {

        glGetDoublev (GL_MODELVIEW_MATRIX, modelview_matrix);

        vx = x;
        vy = trackball_data.screen_height - y;

        glReadBuffer (GL_AUX0);
        glReadPixels (vx, vy, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, col);
        highlighted_ray = col[0] * 256 * 256 + col[1] * 256 + col[2];
        glReadBuffer (GL_BACK);

        plot_rays_redraw ();
    } else {
        rsf_trackball_mouse_event (&trackball_data, button, state, x, y);
        if (trackball_data.reshape) {
            plot_rays_reshape (trackball_data.screen_width, trackball_data.screen_height);
            trackball_data.reshape = 0;
        }
        if (trackball_data.redraw) {
            plot_rays_redraw ();
            trackball_data.redraw = 0;
        }
    }
}

void plot_rays_motion (int x,int y) {

    if (trackball_data.motion_callback) {
        trackball_data.motion_callback (&trackball_data, x, y);
        if (trackball_data.redraw) {
            plot_rays_redraw ();
            trackball_data.redraw = 0;
        }
    }
}

void plot_rays_draw_cross (GLdouble x, GLdouble y, GLdouble z) {
    if (x > 1.0 || x < -1.0 ||
        y > 1.0 || y < -1.0 ||
        z > 1.0 || z < -1.0)
        return;

    glDisable (GL_DEPTH_TEST);
    glDrawBuffer (GL_FRONT);
    glPushAttrib (GL_LINE_BIT);
    glLineWidth (2);
    glEnable (GL_COLOR_LOGIC_OP);
    glLogicOp (GL_XOR);

    glColor3f (1, 1, 1);

    glBegin (GL_LINES);
    glVertex3d (-1.0, y, z);
    glVertex3d (1.0, y, z);
    glVertex3d (x, -1.0, z);
    glVertex3d (x, 1.0, z);
    glVertex3d (x, y, -1.0);
    glVertex3d (x, y, 1.0);
    glEnd ();

    glPopAttrib ();
    glLogicOp (GL_SET);
    glDisable (GL_COLOR_LOGIC_OP);
    glDrawBuffer (GL_BACK);
    glEnable (GL_DEPTH_TEST);
    glFlush ();

    has_drawn_cross = 1;
}


void plot_rays_passive_motion (int x, int y) {
    if (0 == num_aux_buffers)
        return;

    GLdouble vx, vy, vz;
    GLdouble wx, wy, wz;
    GLfloat vvz;
    GLubyte col[3] = { 0, 0, 0 };

    if (has_drawn_cross)
        plot_rays_draw_cross (cross_x, cross_y, cross_z);
    has_drawn_cross = 0;

    vx = x;
    vy = trackball_data.screen_height - y;

    glReadBuffer (GL_AUX0);
    glReadPixels (vx, vy, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, col);
    glReadPixels (vx, vy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &vvz);
    vz = vvz;
    glReadBuffer (GL_BACK);

    if (0 == col[0] && 0 == col[1] && 0 == col[2])
        return;

    glGetDoublev (GL_MODELVIEW_MATRIX, modelview_matrix);
    gluUnProject (vx, vy, vz, modelview_matrix, projection_matrix, viewport_matrix, &wx, &wy, &wz);

    cross_x = wx;
    cross_y = wy;
    cross_z = wz;
    plot_rays_draw_cross (cross_x, cross_y, cross_z);
}

void write_eps(void) 
{
    char * eps;

    eps = sf_getstring("eps");
    /* name of EPS file for hard copy */
    if (NULL == eps) {
	eps = sf_charalloc(12);
	strncpy(eps,"sfrays.eps",12);
    }
    
    FILE *fp = fopen (eps, "wb");
    free(eps);
    
    GLint buffsize = 0, state = GL2PS_OVERFLOW;
    GLint viewport[4];
    
    glGetIntegerv(GL_VIEWPORT, viewport);
    
    is_rendering_eps = 1;
    while (state == GL2PS_OVERFLOW) { 
	buffsize += 1024*1024;
	gl2psBeginPage ("3D Rays", "Madagascar 3D ray visualizer", 
			viewport,
			GL2PS_EPS, GL2PS_BSP_SORT, GL2PS_SILENT |
			GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
			GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
			GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
			fp, "sfrays");
	plot_rays_redraw ();
	state = gl2psEndPage ();
    }
    is_rendering_eps = 0;
}    

void plot_rays_menu (int value) {
    switch (value) {
        case 1: 
	    write_eps ();
	    break;
        case -1:
	    write_eps ();
            exit (0);
            break;
        case 0:
            glutFullScreen ();
            break;
        default:
            break;    
    }
    glutPostRedisplay ();
}

int main (int argc, char **argv) {

    int n1, n2, n3, ir, i, traj_len = 0;
    float o1, o2, o3, d1, d2, d3, *traj = NULL;
    int nz, nx, ny;
    float oz, dz, ox, dx, oy, dy;
    float curr1, curr2, curr3, prev1, prev2, prev3;
    float scale_max;
    int nsr;
    GLubyte ri, gi, bi;
    sf_file in=NULL;
    sf_file frame=NULL;

    sf_init(argc,argv);
    in = sf_input ("in");

    if (NULL != sf_getstring ("frame"))
        frame = sf_input ("frame");
    else
        frame = NULL;

    if (SF_FLOAT != sf_gettype (in))
        sf_error("Need float input");

    if (!sf_histint (in, "n1", &n1) && !sf_getint ("n1", &n1))
        sf_error ("Need n1=");
    if (!sf_histint (in, "n2", &n2) && !sf_getint ("n2", &n2))
        sf_error ("Need n2=");
    if (!sf_histint (in, "n3", &n3) && !sf_getint ("n3", &n3))
        sf_error ("Need n3=");

    if (!sf_histfloat (in, "d1", &d1) && !sf_getfloat ("d1", &d1))
        sf_error ("Need d1="); /* sampling on the first axis */
    if (!sf_histfloat (in, "d2", &d2) && !sf_getfloat ("d2", &d2))
        sf_error ("Need d2=");
    if (!sf_histfloat (in, "d3", &d3) && !sf_getfloat ("d3", &d3))
        sf_error ("Need d3=");

    if (!sf_histfloat (in, "o1", &o1) && !sf_getfloat ("o1", &o1))
        sf_error ("Need o1=");
    if (!sf_histfloat (in, "o2", &o2) && !sf_getfloat ("o2", &o2))
        sf_error ("Need o2=");
    if (!sf_histfloat (in, "o3", &o3) && !sf_getfloat ("o3", &o3))
        sf_error ("Need o3=");

    if (n1 != 3)
        sf_error ("n1 != 3 - does not seem like a 3D rays data");

    if (!sf_getint ("n1", &nz) && 
        (NULL == frame || !sf_histint (frame,"n1", &nz)))
        sf_error ("Need nz=");
    if (!sf_getint ("n2", &nx) &&
        (NULL == frame || !sf_histint (frame,"n2", &nx)))
        sf_error ("Need nx=");
    if (!sf_getint ("n3", &n3) &&
        (NULL == frame || !sf_histint (frame,"n3", &ny)))
        sf_error ("Need ny=");

    if (!sf_getfloat ("d1", &dz) &&
        (NULL == frame || !sf_histfloat (frame, "d1", &dz)))
        sf_error ("Need dz=");
    if (!sf_getfloat ("d2", &dx) &&
        (NULL == frame || !sf_histfloat (frame, "d2", &dx)))
        sf_error ("Need dx=");
    if (!sf_getfloat ("d3", &dy) &&
        (NULL == frame || !sf_histfloat (frame, "d3", &dy)))
        sf_error ("Need dy=");

    if (!sf_getfloat ("o1", &oy) &&
        (NULL == frame || !sf_histfloat (frame, "o1", &oz)))
        sf_error ("Need oz=");
    if (!sf_getfloat ("o2", &ox) &&
        (NULL == frame || !sf_histfloat (frame, "o2", &ox)))
        sf_error ("Need ox=");
    if (!sf_getfloat ("o3", &oy) &&
        (NULL == frame || !sf_histfloat (frame, "o3", &oy)))
        sf_error ("Need oy=");

    traj_len = n2 * 3;
    nsr = n3;

    traj = sf_floatalloc (traj_len);

    if (0 == nsr)
        sf_error ("Zero rays read");

    rsf_trackball_init (&trackball_data);

    glutInit (&argc, argv);
    glutInitWindowSize (768, 768);

#ifdef GLUT_AUX1
    glutInitDisplayMode (GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_AUX1);
#else
    glutInitDisplayMode (GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    sf_warning ("GLUT library on this system does not support aux buffers - mouse selection is disabled");
#endif
    glutCreateWindow ("3D Rays Plotter");
    glutDisplayFunc (plot_rays_redraw);
    glutIdleFunc (NULL);

    trackball (trackball_data.quat, -0.35, 0.35, -0.1, 0);
    build_rotmatrix (trackball_data.m,
                     trackball_data.quat);

    glutReshapeFunc (plot_rays_reshape);
    glutMouseFunc (plot_rays_mouse);
    glutMotionFunc (plot_rays_motion);
    glutPassiveMotionFunc (plot_rays_passive_motion);
    glutCreateMenu (plot_rays_menu);
    glutAddMenuEntry ("Full Screen",0);
    glutAddMenuEntry ("Save view as EPS",1);
    glutAddMenuEntry ("Quit",-1);
/*  glutAddMenuEntry ("White/Color Rays",1);
    glutAddMenuEntry ("Plot Rays",2);
    glutAddMenuEntry ("Surface Traveltimes",3);
    glutAddMenuEntry ("Wired or Solid WFs",4);
    glutAddMenuEntry ("Plot Wavefronts",5);
    glutAddMenuEntry ("TRI or TETRA or LAYER or HORZ",6);*/
    glutAttachMenu (GLUT_RIGHT_BUTTON);

    glShadeModel (GL_FLAT);
    glEnable (GL_DEPTH_TEST);
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glDisable (GL_LIGHTING);
    glPixelStorei (GL_PACK_ALIGNMENT, 1);
    glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
    glPixelTransferi (GL_MAP_COLOR, GL_FALSE);

#ifdef GLUT_AUX1
    glGetIntegerv (GL_AUX_BUFFERS, &num_aux_buffers);
#endif

    rays_list_base = glGenLists (2 * nsr);
    /* Read rays */
    for (ir = 0; ir < nsr; ir++) {
        prev1 = -10.0; prev2 = -10.0; prev3 = -10.0;
        sf_floatread (traj, n2, in);
        glNewList ((rays_list_base + ir) * 2, GL_COMPILE);
        glBegin (GL_LINE_STRIP);
        for (i = 0; i < n2; i++) {
            curr1 = -1.0 + 2.0*(traj[i * 3] - oz)/(float)(nz*dz);
            curr2 = -1.0 + 2.0*(traj[i * 3 + 1] - ox)/(float)(nx*dx);
            curr3 = -1.0 + 2.0*(traj[i * 3 + 2] - oy)/(float)(ny*dy);
            /* If exactly the same point as the previous one is encountered, then - end of the ray */
            if ((curr1 == prev1 && curr2 == prev2 && curr3 == prev3) ||
                curr1 < -1.0 || curr1 > 1.0 || curr2 < -1.0 || curr2 > 1.0 ||
                curr3 < -1.0 || curr3 > 1.0)
                break;
            glVertex3f (curr2, -curr1, curr3);
            prev1 = curr1;
            prev2 = curr2;
            prev3 = curr3;
        }
        glEnd ();
        glEndList ();
        if (num_aux_buffers > 0) {
            /* Encode the index of the ray into color */
            ri = gi = bi = 0;
            bi = (ir + 1) % 256;
            if ((ir + 1) >= 256) {
                if ((ir + 1) >= (256 * 256)) {
                    ri = ((ir + 1) / 256) % 256;
                    gi = ((ir + 1) / 256) / 256;
                } else
                    gi = (ir + 1) / 256;
            }
            glNewList ((rays_list_base + ir) * 2 + 1, GL_COMPILE);
            glColor3ub (ri, gi, bi);
            glCallList ((rays_list_base + ir) * 2);
            glEndList ();
        }
    }
    rays_list_num = nsr;

    scale_x = nx * dx;
    scale_y = nz * dz;
    scale_z = ny * dy;

    scale_max = scale_x;
    if (scale_y > scale_max)
        scale_max = scale_y;
    if (scale_z > scale_max)
        scale_max = scale_z;

    scale_x /= scale_max;
    scale_y /= scale_max;
    scale_z /= scale_max;

    scale_x *= 0.5;
    scale_y *= 1.0;
    scale_z *= 0.7;

    glutMainLoop ();


    return 0;
}
