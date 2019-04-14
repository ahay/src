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

#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
#include <GL/glut.h>
#ifdef FREEGLUT
#include <GL/freeglut_ext.h>
#endif
#endif

#include <limits.h>

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
#define TEX_SIZE 256 /* Texture size  - should be a power of two */
#define LIST_CHUNK 1000 /* Define how many display lists generate at once */
#define MAX_CHUNKS 100000 /* How many chunks of display lists to have */
#define WIN_HEIGHT 760 /* Default window height */
#define DEV_SCALE 10 /* Scale factor between OpenGL and window coordinate systems */
#define MIN_DELAY 50.0 /* Maximum animation delay */
#define MAX_DELAY 10000.0 /* Minimum animation delay */

static float *color_table = NULL;
static bool light = false;
static bool stretchy = false;
static bool aalias = false;
static int aawidth = 1;
static int oglcolor, nx, ny;
static unsigned int *ogllists = NULL; /* Bases of every LIST_CHUNK display lists */
static unsigned char *tex_buf = NULL; /* Temporary buffer for loading textures */
static unsigned int frames_num = 0;
static unsigned int curr_frame = 0;
static float delay = 50; /* milliseconds */
static int direction = 0; /* Animation direction */
static int dirway = 0; /* Animation direction for both ways */
static bool animate = false; /* Animatio flag */
static int menu_tag = 0; /* GLUT menu tag */
static bool has_menu = false;
static int window_id; /* GLUT window identificator */

#define MAX_TITLE_STR 64
static const char *pen_title = "OpenGL pen";
static char title_strbuf[MAX_TITLE_STR];
static int curr_title_frame = -1;

enum {
    MENU_NEXT,
    MENU_PREV,
    MENU_RESTART,
    MENU_RUN,
    MENU_STOP,
    MENU_DIRECTION,
    MENU_STRETCHY,
    MENU_FULLSCREEN,
    MENU_QUIT
};


static void oglquit (int status) {
    glutDestroyWindow (window_id);
    exit (status);
}

void opendev (int argc, char* argv[])
/*< open >*/
{
    const char *color;
    int dwidth, dheight, mwidth; /*, mheight; */
    int win_width, win_height;

    sf_parenv("SFPENOPTS");

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    /* control routines */
    dev.reset = oglreset;
    dev.erase = oglerase;
    dev.close = oglclose;

    dev.area = oglarea;
    dev.raster = oglraster;
    dev.attributes = oglattr;
    dev.plot = oglplot;
    dev.startpoly = oglstartpoly;
    dev.midpoly = oglmidpoly;
    dev.endpoly = oglendpoly;
    dev.point = oglpoint;

    dev.need_end_erase = true;
    dev.smart_clip = true;
    dev.smart_raster = true;
    dev.num_col = NCOLOR;

    color_table = sf_floatalloc (3*NCOLOR);
    tex_buf = sf_ucharalloc (TEX_SIZE*TEX_SIZE);
    ogllists = (unsigned int*) sf_alloc (MAX_CHUNKS,sizeof (unsigned int));

    glutInit (&argc, argv);

    dwidth = glutGet (GLUT_SCREEN_WIDTH);
    dheight = glutGet (GLUT_SCREEN_HEIGHT);
    mwidth = glutGet (GLUT_SCREEN_WIDTH_MM);
/*    mheight = glutGet (GLUT_SCREEN_HEIGHT_MM); */

    win_height = WIN_HEIGHT;
    if (win_height > 0.85*dheight)
        win_height = 0.85*dheight;
    win_width = win_height / VP_SCREEN_RATIO;

    ny = win_height*DEV_SCALE;
    nx = ny / VP_SCREEN_RATIO;

    dev.pixels_per_inch = (float)nx / 
                          ((float)win_width / (float)dwidth*mwidth / 25.4);
    dev.aspect_ratio = 1.0;
    sf_getfloat ("aspect", &dev.aspect_ratio);
    /* aspect ratio */
    sf_getfloat ("ppi", &dev.pixels_per_inch);
    /* pixels per inch */

    dev.xmin = 0;
    dev.ymin = 0;
    dev.xmax = nx;
    dev.ymax = ny;

    invras = false;

    glutInitWindowSize (win_width, win_height);

    glutInitDisplayMode (GLUT_RGBA | GLUT_DOUBLE | GLUT_STENCIL);
    if (!glutGet (GLUT_DISPLAY_MODE_POSSIBLE)) {
        fprintf (stderr, "OpenGL device does not support stencil buffer: using genpatarea ()\n");
        dev.area = genpatarea;
        glutInitDisplayMode (GLUT_RGBA | GLUT_DOUBLE);
        if (!glutGet (GLUT_DISPLAY_MODE_POSSIBLE)) {
            fprintf (stderr, "OpenGL device does not support double buffer: trying single buffer\n");
            glutInitDisplayMode (GLUT_RGBA | GLUT_SINGLE);
            if (!glutGet (GLUT_DISPLAY_MODE_POSSIBLE)) {
                fprintf (stderr, "OpenGL device device initialization failed\n");
                exit (-1);
            }
        }
    }

    window_id = glutCreateWindow (pen_title);
    glutDisplayFunc (oglredraw);
    glutKeyboardFunc (oglkeyboard);
    glutIdleFunc (NULL);
    glutReshapeFunc (oglreshape);

    if (NULL == (color = sf_getstring ("bgcolor"))) color = "black";
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

    if (!sf_getbool ("stretchy", &stretchy)) stretchy = false;
    /* y - stretch plot to all borders of the drawing area */

    if (!sf_getbool ("aalias", &aalias)) aalias = false;
    /* y - draw lines, points, polygons with antialialing */
    if (!sf_getint ("aawidth", &aawidth)) aawidth = 1;
    /* width of transition areas when aalias=y */
}

void oglbuildmenu (void)
/*< menu builder >*/
{
    if (has_menu)
        glutDestroyMenu (menu_tag);
    else
        has_menu = true;

    menu_tag = glutCreateMenu (oglmenu);
    if (frames_num > 1 || animate) {
        if (!animate) {
            glutAddMenuEntry ("Next", MENU_NEXT);
            glutAddMenuEntry ("Prev", MENU_PREV);
            glutAddMenuEntry ("Restart", MENU_RESTART);
            glutAddMenuEntry ("Run", MENU_RUN);
        } else
            glutAddMenuEntry ("Stop", MENU_STOP);
        glutAddMenuEntry (direction == 2 ? "Forwards" :
                          direction == 1 ? "Both ways" : "Backwards", MENU_DIRECTION);
    }
    glutAddMenuEntry (stretchy ? "Rigid" : "Stretchy", MENU_STRETCHY);
    glutAddMenuEntry ("Full Screen", MENU_FULLSCREEN);
    glutAddMenuEntry ("Quit", MENU_QUIT);
    glutAttachMenu (GLUT_RIGHT_BUTTON);
}

void oglreset (void)
/*< reset >*/
{
    int value, r, g, b;

    /* reset color table */
    color_table[0] = 1.0*light;
    color_table[NCOLOR] = 1.0*light;
    color_table[NCOLOR*2] = 1.0*light;

    for (value = 1; value < 8; value++) {
        r = MAX_GUN*((value & 2) / 2);
        g = MAX_GUN*((value & 4) / 4);
        b = MAX_GUN*((value & 1) / 1);

        if (light) {
            color_table[value] = (MAX_GUN - r) / (float)MAX_GUN;
            color_table[NCOLOR + value] = (MAX_GUN - g) / (float)MAX_GUN;
            color_table[NCOLOR*2 + value] = (MAX_GUN - b) / (float)MAX_GUN;
        } else {
            color_table[value] = r / (float)MAX_GUN;
            color_table[NCOLOR + value] = g / (float)MAX_GUN;
            color_table[NCOLOR*2 + value] = b / (float)MAX_GUN;
        }
    }
 
    for (value = 8; value < NCOLOR; value++) {
        color_table[value] = 1.0*light;
        color_table[NCOLOR + value] = 1.0*light;
        color_table[NCOLOR*2 + value] = 1.0*light;
    }

    glDisable (GL_DEPTH_TEST);
    glDepthFunc (GL_ALWAYS);
    glShadeModel (GL_FLAT);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glDisable (GL_CULL_FACE);
    glDisable (GL_LIGHTING);
    glDisable (GL_SCISSOR_TEST);
    glDisable (GL_STENCIL_TEST);
    glPixelStorei (GL_PACK_ALIGNMENT, 1);
    glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
    glClearStencil (0);
    if (aalias) {
        glEnable (GL_LINE_SMOOTH);
        glEnable (GL_POINT_SMOOTH);
        glEnable (GL_POLYGON_SMOOTH);
        glLineWidth (aawidth);
        glPointSize (aawidth);
    }
    glDisable (GL_CLIP_PLANE0);
    glDisable (GL_CLIP_PLANE1);
    glDisable (GL_CLIP_PLANE2);
    glDisable (GL_CLIP_PLANE3);
}

/* Process queued events from the event loop */
static void oglprocessevents (void) {
#if defined(__APPLE__) && (GLUT_MACOSX_IMPLEMENTATION >= 2)
    if (!animate && frames_num == 1) {
        oglstartanimation ();
        oglbuildmenu ();
    }
    glutCheckLoop ();
#elif defined(FREEGLUT)
    if (!animate && frames_num == 1) {
        oglstartanimation ();
        oglbuildmenu ();
    }
    glutMainLoopEvent ();
#endif
}

void oglerase (int command)
/*< erase >*/
{
    switch (command) {
        case ERASE_START:
            ogllists[0] = glGenLists (LIST_CHUNK*2);
            glNewList (ogllists[0], GL_COMPILE_AND_EXECUTE);
            glClear (GL_COLOR_BUFFER_BIT);
            frames_num++;
            break;
        case ERASE_MIDDLE:
            glFlush ();
            glEndList ();
            /* Generate another display list - for color background */
            glNewList (ogllists[frames_num / LIST_CHUNK] +
                       (frames_num % LIST_CHUNK)*2 - 1,
                       GL_COMPILE_AND_EXECUTE);
            glClearColor (color_table[0],
                          color_table[NCOLOR],
                          color_table[NCOLOR*2], 0.0f);
            glEndList ();
            if (1 == frames_num)
                fprintf (stderr, "Preparing animation\n");
            else if (frames_num % 50 == 0)
                fprintf (stderr, "Finished %d frames\n", frames_num);
            oglprocessevents ();
            if (!(frames_num % LIST_CHUNK))
                ogllists[frames_num / LIST_CHUNK] = glGenLists (LIST_CHUNK*2);
            /* Main display list - all plotting */
            glNewList (ogllists[frames_num / LIST_CHUNK] +
                       (frames_num % LIST_CHUNK)*2,
                       GL_COMPILE_AND_EXECUTE);
            glClear (GL_COLOR_BUFFER_BIT);
            frames_num++;
            break;
        case ERASE_END:
            glFlush ();
            glEndList ();
            /* Generate another display list - for color background */
            glNewList (ogllists[frames_num / LIST_CHUNK] +
                       (frames_num % LIST_CHUNK)*2 - 1,
                       GL_COMPILE_AND_EXECUTE);
            glClearColor (color_table[0],
                          color_table[NCOLOR],
                          color_table[NCOLOR*2], 0.0f);
            glEndList ();
            if (frames_num > 1)
                fprintf (stderr, "Finished all frames\n");
            break;
        default:
            break;
    }
    /* Reset clipping */
    glDisable (GL_CLIP_PLANE0);
    glDisable (GL_CLIP_PLANE1);
    glDisable (GL_CLIP_PLANE2);
    glDisable (GL_CLIP_PLANE3);
}

void oglstartanimation ()
{
/*< start animation >*/
    animate = true;
    glutTimerFunc (delay, oglanimate, 0);
}

void oglanimate (int value)
/*< animation >*/
{
    if (!animate)
        return;

    switch (direction) {
        case 0: /* Forwards */
            if (curr_frame == (frames_num - 1))
                curr_frame = 0;
            else
                curr_frame++;
            break;
        case 1: /* Backwards */
            if (curr_frame == 0)
                curr_frame = frames_num - 1;
            else
                curr_frame--;
            break;
        case 2: /* Both ways */
            if (dirway == 0 && curr_frame == (frames_num - 1)) {
                dirway = 1;
                curr_frame = frames_num - 1;
            } if (dirway == 1 && curr_frame == 0) {
                dirway = 0;
                curr_frame = 1;
            } else {
                curr_frame += dirway ? -1 : 1;
            }
            break;
    }

    glutPostRedisplay ();
    glutTimerFunc (delay, oglanimate, 0);
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
            if (frames_num > 1 && !animate)
                oglstartanimation ();
            oglbuildmenu ();
#if defined(__APPLE__) && (GLUT_MACOSX_IMPLEMENTATION >= 2)
            for (;;)
                glutCheckLoop ();
#else
            glutMainLoop ();
#endif
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
    GLdouble eq[4];

    switch (command) {
        case SET_COLOR:
            oglcolor = value;
            glColor3f (color_table[oglcolor],
                       color_table[NCOLOR + oglcolor],
                       color_table[NCOLOR*2 + oglcolor]);
            break;
        case SET_COLOR_TABLE:
            color_table[value] = v1 / (float)MAX_GUN;
            color_table[NCOLOR + value] = v2 / (float)MAX_GUN;
            color_table[NCOLOR*2 + value] = v3 / (float)MAX_GUN;
            if (value == oglcolor)
                glColor3f (color_table[oglcolor],
                           color_table[NCOLOR + oglcolor],
                           color_table[NCOLOR*2 + oglcolor]);
            break;
        case SET_WINDOW:
            /* Right */
            glEnable (GL_CLIP_PLANE0);
            eq[0] = -1; eq[1] = 0; eq[2] = 0; eq[3] = v2 + 0.5;
            glClipPlane (GL_CLIP_PLANE0, eq);
            /* Left */
            glEnable (GL_CLIP_PLANE1);
            eq[0] = 1; eq[1] = 0; eq[2] = 0; eq[3] = -value + 0.5;
            glClipPlane (GL_CLIP_PLANE1, eq);
            /* Up */
            glEnable (GL_CLIP_PLANE2);
            eq[0] = 0; eq[1] = -1; eq[2] = 0; eq[3] = v3 + 0.5;
            glClipPlane (GL_CLIP_PLANE2, eq);
            /* Down */
            glEnable (GL_CLIP_PLANE3);
            eq[0] = 0; eq[1] = 1; eq[2] = 0; eq[3] = -v1 + 0.5;
            glClipPlane (GL_CLIP_PLANE3, eq);
            break;
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

void oglstartpoly (int npts)
/*< startploy >*/
{
    glBegin (GL_LINE_STRIP);
}

void oglmidpoly (int x, int y)
/*< midploy >*/
{
    glVertex2i (x, y);
}

void oglendpoly (int last)
/*< midploy >*/
{
    glEnd ();
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
    int i, x, y;
    int xvmin = INT_MAX, xvmax = INT_MIN,
        yvmin = INT_MAX, yvmax = INT_MIN;

    /* Use stencil buffer to implement the even-odd rule */
    glColorMask (GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glEnable (GL_STENCIL_TEST);
    glClear (GL_STENCIL_BUFFER_BIT);
    glStencilMask (0x01);
    glStencilFunc (GL_ALWAYS, 0, 0);
    /* First - calculate self-intersections in the stencil buffer */
    glStencilOp (GL_KEEP, GL_KEEP, GL_INVERT);
    glBegin (GL_TRIANGLE_FAN);
    for (i = 0; i < npts && i < MAXVERT; i++) {
        x = head->x;
        y = head->y;
        glVertex2i (x, y);
        head = head->next;
        if (x < xvmin)
            xvmin = x;
        if (x > xvmax)
            xvmax = x;
        if (y < yvmin)
            yvmin = y;
        if (y > yvmax)
            yvmax = y;
    }
    glEnd ();
    glStencilOp (GL_KEEP, GL_KEEP, GL_KEEP);
    /* Second - draw polygons only where stencil is 1 - odd number of intersection */
    glStencilFunc (GL_EQUAL, 1, 1);
    glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glBegin (GL_QUADS);
    glVertex2i (xvmin, yvmin);
    glVertex2i (xvmax, yvmin);
    glVertex2i (xvmax, yvmax);
    glVertex2i (xvmin, yvmax);
    glEnd ();
    glClear (GL_STENCIL_BUFFER_BIT);
    glDisable (GL_STENCIL_TEST);
}

void oglraster (int xpix, int ypix, int xmin, int ymin, int xmax, int ymax, 
                unsigned char **raster_block, int orient, int dither_it)
{
/*< raster >*/
    int i, j, k, l, hnum, wnum, kmax, lmax, width, height;
    GLuint *tex_id;
    float stepy, stepx;
    float xtex, ytex;
    float x_min, x_max, y_min, y_max;

    /* Colormap table */
    if (!dither_it) {
        glPixelTransferi (GL_MAP_COLOR, GL_TRUE);
        glPixelMapfv (GL_PIXEL_MAP_I_TO_R, NCOLOR, color_table);
        glPixelMapfv (GL_PIXEL_MAP_I_TO_G, NCOLOR, &color_table[NCOLOR]);
        glPixelMapfv (GL_PIXEL_MAP_I_TO_B, NCOLOR, &color_table[NCOLOR * 2]);
    }

    height = orient % 2 ? xpix : ypix;
    width = orient % 2 ? ypix : xpix;

    wnum = width / TEX_SIZE + ((width % TEX_SIZE) != 0);
    hnum = height / TEX_SIZE + ((height % TEX_SIZE) != 0);
    stepy = (ymax - ymin + 1) / (height / (float)TEX_SIZE);
    stepx = (xmax - xmin + 1) / (width / (float)TEX_SIZE);

    tex_id = (GLuint*) sf_alloc (hnum*wnum,sizeof (GLuint));
    glGenTextures (hnum*wnum, tex_id);

    /* Loop over tiles and make textures */
    for (i = 0; i < wnum; i++) {
        kmax = i != (wnum - 1) ? TEX_SIZE : width - i * TEX_SIZE;
        for (j = 0; j < hnum; j++) {
            lmax = j != (hnum - 1) ? TEX_SIZE : height - j * TEX_SIZE;
            /* Copy samples to the local buffer */
            switch (orient) {
                case 0:
                    for (l = 0; l < lmax; l++) {
                        for (k = 0; k < kmax; k++) {
                            tex_buf[k*TEX_SIZE + l] =
                            raster_block[height - j*TEX_SIZE - l - 1][i*TEX_SIZE + k];
                        }
                    }
                    break;
                case 1:
                    for (k = 0; k < kmax; k++) {
                        for (l = 0; l < lmax; l++) {
                            tex_buf[k*TEX_SIZE + l] =
                            raster_block[width - i*TEX_SIZE - k - 1][height - j*TEX_SIZE - l - 1];
                        }
                    }
                    break;
                case 2:
                    for (l = 0; l < lmax; l++) {
                        for (k = 0; k < kmax; k++) {
                            tex_buf[k*TEX_SIZE + l] =
                            raster_block[j*TEX_SIZE + l][width - i*TEX_SIZE - k - 1];
                        }
                    }
                    break;
                case 3:
                    for (k = 0; k < kmax; k++) {
                        for (l = 0; l < lmax; l++) {
                            tex_buf[k*TEX_SIZE + l] =
                            raster_block[i*TEX_SIZE + k][j*TEX_SIZE + l];
                        }
                    }
                    break;
            }

            glBindTexture (GL_TEXTURE_2D, tex_id[i*hnum + j]);
            glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
            glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
            glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
            glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            /* Send texture samples to OpenGL */
            glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, TEX_SIZE, TEX_SIZE, 0,
                                         dither_it ? GL_LUMINANCE : 
                                                     GL_COLOR_INDEX,
                                         GL_UNSIGNED_BYTE, tex_buf);
            glBindTexture (GL_TEXTURE_2D, 0);
        }
    }
    /* Loop over tiles and render textures */
    for (i = 0; i < wnum; i++) {
        x_min = xmin + i*stepx;
        x_max = i != (wnum - 1) ? xmin + (i + 1)*stepx : xmax;
        xtex = i != (wnum - 1) ? 1.0 : 1.0*(width - i*TEX_SIZE) / (float)TEX_SIZE;
        for (j = 0; j < hnum; j++) {
            y_min = ymin + j*stepy;
            y_max = j != (hnum - 1) ? ymin + (j + 1)*stepy : ymax;
            ytex = j != (hnum - 1) ? 1.0 : 1.0*(height - j*TEX_SIZE) / (float)TEX_SIZE;
            glEnable (GL_TEXTURE_2D);
            glBindTexture (GL_TEXTURE_2D, tex_id[i*hnum + j]);
            glBegin (GL_QUADS);
            glTexCoord2f (0, 0);
            glVertex2i (x_min, y_min);
            glTexCoord2f (0, xtex);
            glVertex2i (x_max, y_min);
            glTexCoord2f (ytex, xtex);
            glVertex2i (x_max, y_max);
            glTexCoord2f (ytex, 0);
            glVertex2i (x_min, y_max);
            glEnd ();
            glBindTexture (GL_TEXTURE_2D, 0);
            glDisable (GL_TEXTURE_2D);
        }
    }

    if (!dither_it)
        glPixelTransferi (GL_MAP_COLOR, GL_FALSE);

    free (tex_id);
}

void oglredraw (void)
/*< OpenGL redraw >*/
{
    /* Set background color */
    glCallList (ogllists[curr_frame / LIST_CHUNK] +
                (curr_frame % LIST_CHUNK)*2 + 1);
    /* Plotting */
    glCallList (ogllists[curr_frame / LIST_CHUNK] +
                (curr_frame % LIST_CHUNK)*2);
    glutSwapBuffers ();

    if (frames_num > 1 && curr_title_frame != (int) curr_frame) {
        curr_title_frame = curr_frame;
        snprintf (title_strbuf, MAX_TITLE_STR - 1, "%s: frame %d of %d",
                  pen_title, curr_title_frame + 1, frames_num);
        glutSetWindowTitle (title_strbuf);
    }
}

void oglreshape (int width, int height)
/*< OpenGL reshape >*/
{
    GLdouble wh = 1.0;

    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    if (!stretchy) {
        if (VP_ABSOLUTE == size) {
            ny = height*DEV_SCALE;
            nx = ny*width/height;
            glOrtho (-1, nx,
                     -1, ny, -1, 1);
        } else {
            if (width*VP_SCREEN_RATIO <= height) {
                wh = (GLfloat)height / (GLfloat)width / VP_SCREEN_RATIO;
                glOrtho (-1, nx,
                         -1, ny*wh, -1, 1);
            } else {
                wh = ((GLfloat)width / (GLfloat)height)*VP_SCREEN_RATIO;
                glOrtho (-1, nx*wh,
                         -1, ny, -1, 1);
            }
        }
    } else
        glOrtho (-1, nx,
                 -1, ny, -1, 1);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    glViewport (0, 0, width, height);

    glutPostRedisplay ();
}

void oglmenu (int value)
/*< OpenGL menu >*/
{
    switch (value) {
        case MENU_NEXT:
            if (!animate) {
                if (curr_frame == (frames_num - 1))
                    curr_frame = 0;
                else
                    curr_frame++;
                break;
            }
            break;
        case MENU_PREV:
            if (!animate) {
                if (curr_frame == 0)
                    curr_frame = frames_num - 1;
                else
                    curr_frame--;
            }
            break;
        case MENU_RESTART:
            if (!animate) {
                if (direction == 1)
                    curr_frame = frames_num - 1;
                else
                    curr_frame = 0;
            }
            break;
        case MENU_RUN:
            oglstartanimation ();
            oglbuildmenu ();
            break;
        case MENU_STOP:
            animate = false;
            oglbuildmenu ();
            break;
        case MENU_DIRECTION:
            direction++;
            if (direction > 2)
                direction = 0;
            oglbuildmenu ();
            break;
        case MENU_STRETCHY:
            stretchy = (bool) !stretchy;
            oglreshape (glutGet (GLUT_WINDOW_WIDTH), glutGet (GLUT_WINDOW_HEIGHT));
            oglbuildmenu ();
            break;
        case MENU_FULLSCREEN:
            glutFullScreen ();
            break;
        case MENU_QUIT:
            oglquit (0);
            break;
        default:
            break;    
    }
    glutPostRedisplay ();
}

void oglkeyboard (unsigned char key, int x, int y)
/*< OpenGL keyboard callback >*/
{
    switch (key) {
        case 'F': /* Fast */
        case 'f':
            if (frames_num > 1) {
                delay *= 0.5;
                if (delay < MIN_DELAY)
                    delay = MIN_DELAY;
                fprintf (stderr, "Animation delay [ms]: %f\n", delay);
            }
            break;
        case 'S': /* Slow */
        case 's':
            if (frames_num > 1) {
                delay *= 2.0;
                if (delay > MAX_DELAY)
                    delay = MAX_DELAY;
                fprintf (stderr, "Animation delay [ms]: %f\n", delay);
            }
            break;
        case 'N': /* Next */
        case 'n':
	    if (frames_num > 1) {
		if (animate) {
		    animate = false;
		    oglbuildmenu ();
		}
		if (curr_frame == (frames_num - 1))
                    curr_frame = 0;
                else
                    curr_frame++;
                break;
            }
            break;
        case 'm': /* Prev */
        case 'M':
            if (frames_num > 1) {
		if (animate) {
		    animate = false;
		    oglbuildmenu ();
		}
                if (curr_frame == 0)
                    curr_frame = frames_num - 1;
                else
                    curr_frame--;
            }
            break;
        case 'R': /* Run */
        case 'r':
            if (!animate && frames_num > 1) {
                oglstartanimation ();
                oglbuildmenu ();
            }
            break;
        case '.': /* Stop */
            if (animate && frames_num > 1) {
                animate = false;
                oglbuildmenu ();
            }
            break;
        case 'D': /* Direction */
        case 'd':
            if (frames_num > 1) {
                direction++;
                if (direction > 2)
                    direction = 0;
                fprintf (stderr, "Animation direction: %s\n",
                         direction == 2 ? "Both ways" :
                         direction == 1 ? "Backwards" : "Forwards");
                oglbuildmenu ();
            }
            break;
        case 'T': /* Stretchy */
        case 't':
            stretchy = (bool) !stretchy;
            oglreshape (glutGet (GLUT_WINDOW_WIDTH), 
			glutGet (GLUT_WINDOW_HEIGHT));
            oglbuildmenu ();
            break;
        case 'Q': /* Quit */
        case 'q':
            oglquit (0);
            break;
        default:
            break;    
    }
    glutPostRedisplay ();
}

