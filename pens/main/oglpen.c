/* vplot filter for OpenGL. */
/*
  Copyright (C) 2009 The University of Texas at Austin
  
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

#if defined(__APPLE__)&& defined(__MACH__)
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "../include/attrcom.h"
#include "../include/extern.h"
#include "../include/params.h"
#include "../include/err.h"
#include "../include/erasecom.h"
#include "../include/closestat.h"
#include "../include/enum.h"

#include "../genlib/genpen.h"
#include "../utilities/util.h"

#include "dovplot.h"
#include "init_vplot.h"

#include "oglpen.h"

char            name[] = "oglpen";
#include "ogldoc.h"

#include "_device.h"

#define PPI 100.0  /* pixels per inch */
#define NCOLOR 256 /* number of colors */
#define MAXVERT 1000

static float *color_table = NULL;
static bool light = false;
static int oglcolor, delay, nx, ny;
static unsigned int ogllist;

void opendev (int argc, char* argv[])
/*< open >*/
{
    char *color;
    int dwidth, dheight, mwidth, mheight;

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    /* control routines */
    dev.reset = oglreset;
    dev.erase = oglerase;
    dev.close = oglclose;

    dev.area = oglarea;
    dev.attributes = oglattr;
    dev.plot = oglplot;
    dev.point = oglpoint;

    dev.need_end_erase = true;
    dev.smart_clip = false;
    dev.num_col = NCOLOR;

    color_table = (float*)malloc (3 * NCOLOR * sizeof (float));

    glutInit (&argc, argv);
    glutInitWindowSize (760 / VP_SCREEN_RATIO, 760);

    glutInitDisplayMode (GLUT_RGBA | GLUT_DOUBLE);
    glutCreateWindow ("OpenGL pen");
    glutDisplayFunc (oglredraw);
    glutIdleFunc (NULL);

    glutReshapeFunc (oglreshape);
    glutCreateMenu (oglmenu);
    glutAddMenuEntry ("Full Screen", 0);
    glutAddMenuEntry ("Quit", -1);
    glutAttachMenu (GLUT_RIGHT_BUTTON);

    dwidth = glutGet (GLUT_SCREEN_WIDTH);
    dheight = glutGet (GLUT_SCREEN_HEIGHT);
    mwidth = glutGet (GLUT_SCREEN_WIDTH_MM);
    mheight = glutGet (GLUT_SCREEN_HEIGHT_MM);
    dev.pixels_per_inch = ((float)dwidth / (float)mwidth) * 25.4 ;
    dev.aspect_ratio = (((float)dheight / (float)mheight) * 25.4) / dev.pixels_per_inch;
    /*
     * Since X rounds values to integer millimeters, the aspect ratio
     * has some error in it. Push aspect ratio to nearest round percent.
     */
    dev.aspect_ratio = ((int)((dev.aspect_ratio * 100.) + .5)) / 100.;

    nx = dheight * 10; ny = VP_SCREEN_RATIO * dheight * 10;

    sf_getfloat ("aspect", &dev.aspect_ratio);
    /* aspect ratio */
    sf_getfloat ("ppi", &dev.pixels_per_inch);
    /* pixels per inch */

    dev.xmax = nx - 1;
    dev.ymax = ny - 1;

    if (NULL == (color = sf_getstring ("bgcolor"))) color= "black";
    /* background color (black,white,light,dark) */

    switch (color[0]) {
        case 'b': /* black */
        case 'd': /* dark */
        default:
            light = false;
            break;
        case 'w': /* white */
        case 'l': /* light */
            light = true;
            break;
    }

    if (!sf_getint ("delay", &delay)) delay = 10;
    /* animation delay */
}

void oglreset (void)
/*< reset >*/
{
    int value, r, g, b;

    /* reset color table */
    color_table[0] = 1.0 * light;
    color_table[1] = 1.0 * light;
    color_table[2] = 1.0 * light;

    for (value = 1; value < 8; value++) {
	r = MAX_GUN * ((value & 2) / 2);
	g = MAX_GUN * ((value & 4) / 4);
	b = MAX_GUN * ((value & 1) / 1);

	if (light) {
	    color_table[value * 3] = (255 - r) / 255.0;
	    color_table[value * 3 + 1] = (255 - g) / 255.0;
	    color_table[value * 3 + 2] = (255 - b) / 255.0;
	} else {
	    color_table[value * 3] = r / 255.0;
	    color_table[value * 3 + 1] = g / 255.0;
	    color_table[value * 3 + 2] = b / 255.0;
	}
    }
 
    for (value = 8; value < NCOLOR; value++) {
        color_table[value * 3] = 1.0 * light;
        color_table[value * 3 + 1] = 1.0 * light;
        color_table[value * 3 + 2] = 1.0 * light;
    }

    glShadeModel (GL_SMOOTH);
    glDisable (GL_DEPTH_TEST);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glDisable (GL_LIGHTING);
    glDisable (GL_SCISSOR_TEST);
    glClearColor (color_table[0],
                  color_table[1],
                  color_table[2], 0);
}

void oglerase (int command)
/*< erase >*/
{
    switch (command) {
        case ERASE_START:
            ogllist = glGenLists (1);
            glNewList (ogllist, GL_COMPILE);
            glClear (GL_COLOR_BUFFER_BIT);
            break;
        case ERASE_MIDDLE:
            glClear (GL_COLOR_BUFFER_BIT);
            break;
        case ERASE_END:
            glFlush ();
            glEndList ();
            break;
        default:
            break;
    }
}

void oglclose (int status)
/*< close >*/
{
    switch (status) {
        case CLOSE_PAUSE:
            break;
        case CLOSE_ERROR:
        case CLOSE_INTERRUPT:
        case CLOSE_NORMAL:
            glutMainLoop ();
            break;
        case CLOSE_FLUSH:
            glFlush ();
        default:
            break;
    }
}

void oglattr (int command, int value, int v1, int v2, int v3)
/*< attr >*/
{
/*
    int xmin, ymin, xmax, ymax;
*/
    switch (command) {
        case SET_COLOR:
            oglcolor = value;
            glColor3fv (&color_table[oglcolor * 3]);
            break;
        case SET_COLOR_TABLE:
            color_table[value * 3] = v1 / 255.0;
            color_table[value * 3 + 1] = v2 / 255.0;
            color_table[value * 3 + 2] = v3 / 255.0;
            break;
/*
        case SET_WINDOW:
            xmin = value;
            ymin = dev.ymax - v3;
            xmax = v2;
            ymax = dev.ymax - v1;
            glEnable (GL_SCISSOR_TEST);
            glScissor (xmin, ymin, xmax - xmin, ymin - ymax);
            break;*/
        default:
            break;
    }
}

void oglplot(int x, int y, int draw)
/*< plot >*/
{
    static int oldx = 0, oldy = 0;

    if (draw) {
        glBegin (GL_LINES);
        glVertex2i (oldx, oldy);
        glVertex2i (x, y);
        glEnd ();
    } else {
        dev.lost = 0;
    }

    oldx = x;
    oldy = y;
}

void oglpoint(int x, int y)
/*< point >*/
{
    glBegin (GL_POINTS);
    glVertex2i (x, y);
    glEnd ();
}

void oglarea (int npts, struct vertex *head)
/*< area >*/
{
    int i;

    glBegin (GL_POLYGON);

    /* translate data structures */
    for (i = 0; i < npts && i < MAXVERT; i++) {
        head = head->next;
        glVertex2i (head->x, head->y);
    }

    glEnd ();
}

void oglredraw (void)
/*< OpenGL redraw >*/
{
    glCallList (ogllist);
    glFlush ();
    glutSwapBuffers ();
}

void oglreshape (int width, int height)
/*< OpenGL reshape >*/
{
    GLdouble wh = 1.0;

    glViewport (0, 0, width, height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    if (width <= height) {
        wh = (GLfloat)height / (GLfloat)width * VP_SCREEN_RATIO;
        glOrtho (0, nx,
                 0, ny * wh, -1, 1);
    } else {
        wh = (GLfloat)width / (GLfloat)height * VP_SCREEN_RATIO;
        glOrtho (0, nx * wh,
                 0, ny, -1, 1);
    }
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    glutPostRedisplay ();
}

void oglmenu (int value)
/*< OpenGL menu >*/
{
    switch (value) {
        case -1:
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
