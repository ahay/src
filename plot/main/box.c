/* Draw a balloon-style label.

Takes: > out.vpl*/
/*
  Copyright (C) 1987 The Board of Trustees of Stanford University
  
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
/* Original author: Joe Dellinger. */

#include <math.h>
#include <rsf.h>
#include <rsfplot.h>

static float longi, lati;

static void proj (float xin, float yin, float zin,
                  float *x,  float *y,  float *z);
static void rot (float ang, float xin, float yin, float *x, float *y);

int main (int argc, char* argv[])
{
    const int oval_pts=400;
    float angle, x0,y0, xt,yt, x_oval,y_oval, scale0;
    float xarray[3], yarray[3], tx_height, xout, yout, zout;
    float xpath, ypath, xup, yup, scalet, length, norm;
    float xdis, ydis, pscale, x,y;
    int color, fat, ii, font;
    bool pointer, rev, boxit;
    const char *string;

    sf_init(argc,argv);

    vp_init();

    if (!sf_getint("lab_color",&color)) color=VP_WHITE;
    /* label color */
    if (!sf_getint("lab_fat",&fat)) fat=0;
    /* label fatness */

    if (!sf_getfloat("pscale",&pscale)) pscale=1.;
    /* scale factor for width of pointer */
    if (!sf_getbool("pointer",&pointer)) pointer=true;
    /* if y, create arrow pointer */
    if (!sf_getbool("reverse",&rev)) rev=false;
    /* if reverse */

    if (!sf_getfloat("lat",&lati)) lati=0.;
    /* latitude of viewpoint in 3-D */
    if (!sf_getfloat("long",&longi)) longi=90.;
    /* longitude of viewpoint in 3-D */

    if (!sf_getfloat("angle",&angle)) angle=0.;
    /* longitude of floating label in 3-D */
    angle *= SF_PI / 180.;

    if (!sf_getfloat("x0",&x0)) x0=0.;
    /* x position of the pointer tip */
    if (!sf_getfloat("y0",&y0)) y0=0.;
    /* y position of the pointer tip */

    if (!sf_getfloat("scale0",&scale0)) scale0=1.;
    /* scale factor for x0 and y0 */
    x0 *= scale0;
    y0 *= scale0;

    if (!sf_getfloat("xt",&xt)) xt=2.;
    /* relative position of text in x */
    if (!sf_getfloat("yt",&yt)) yt=0.;
    /* relative position of text in y */

    if (!sf_getfloat("x_oval",&x_oval)) x_oval=0.;
    /* x size of the oval around pointer */
    if (!sf_getfloat("y_oval",&y_oval)) y_oval=0.;
    /* y size of the oval around pointer */

    if (!sf_getbool("boxit",&boxit)) boxit=true;
    /* if y, create a box around text */

    if (xt != 0. || yt != 0.) {
        if (sf_getfloat ("length",&length)) {
            /* normalization for xt and yt */
            scalet = length / hypotf(xt,yt);
        } else if (!sf_getfloat ("scalet",&scalet)) {
            /*( scalet scale factor for xt and yt (if length is not set) )*/
            scalet=1.;
        }
        xt *= scalet;
        yt *= scalet;
    }

    xt += x0;
    yt += y0;

    if (!sf_getfloat("size",&tx_height)) tx_height = .25;
    /* text height in inches */

    if (NULL == (string=sf_getstring("label"))) string=" ";
    /* text for label */

    if (!sf_getint("font",&font)) font=VP_NO_CHANGE;
    /* text font */

    proj (tx_height * cosf(angle),
          tx_height * sinf(angle),
          0., &xout, &yout, &zout);

    if (((xout < 0.) && rev) || ((xout >= 0.) && !rev)) {
        xpath = xout;
        ypath = yout;
    } else {
        xpath = -xout;
        ypath = -yout;
    }

    proj (0., 0., tx_height, &xout, &yout, &zout);
    xup = xout;
    yup = yout;

    vp_tjust (TH_CENTER, TV_HALF);

    if(x_oval >0. && y_oval >0.){
        vp_color (color);
        vp_penup ();
        vp_move(x0-x_oval,y0);
        vp_fat (fat);
        x=-x_oval;

        for(ii=0; ii < oval_pts; ii++) {
            y=(1.-x*x/x_oval/x_oval);
            if(y<0) y=0;
            y=y_oval*sqrtf(y);
            vp_pendn(x0+x,y+y0);
            x+=2*x_oval/oval_pts;
        }
        for(ii=0; ii < oval_pts+1; ii++) {
            y=(1.-x*x/x_oval/x_oval);
            if(y<0) y=0;
            y=y_oval*sqrtf(y);
            vp_pendn(x0+x,y0-y);
            x-=2*x_oval/oval_pts;
        }
    }

    if (pointer) {
        xdis = (yt - y0);
        ydis = -(xt - x0);
        norm = 3. * hypotf(xdis,ydis);
        if (norm == 0.) {
            xdis = 0.;
            ydis = 0.;
        } else {
            xdis *= pscale * tx_height / norm;
            ydis *= pscale * tx_height / norm;
        }

        xarray[0] = xt + xdis;
        yarray[0] = yt + ydis;
        xarray[1] = x0;
        yarray[1] = y0;
        xarray[2] = xt - xdis;
        yarray[2] = yt - ydis;

        vp_color (VP_BLACK);
        vp_fat (0);
        vp_fill (xarray, yarray, 3);

        /* Now the dark edge to the pointer */
        vp_color (color);
        vp_fat (fat);

        vp_penup ();
        for (ii = 0; ii < 3; ii++) {
            vp_pendn (xarray[ii], yarray[ii]);
        }
    }

    /* Finally the shaded box with the text */
    vp_color (color);
    vp_fat (fat);

    vp_tfont (font, VP_NO_CHANGE, boxit? OVLY_SHADE_BOX: OVLY_NORMAL);
    vp_gtext (xt, yt, xpath, ypath, xup, yup, string);

    vp_tjust (TH_NORMAL, TV_NORMAL);
    vp_tfont (font, VP_NO_CHANGE, OVLY_NORMAL);
    vp_color (VP_WHITE);
    vp_fat (0);

    exit(0);
}

static void proj (float xin, float yin, float zin,
                  float *x, float *y, float *z)
{
    rot (longi, xin, yin, z, x);
    rot (lati, xin, zin, z, y);
    rot (0., *x, *z, x, z);
}

static void rot (float ang, float xin, float yin, float *x, float *y)
{
    float c, s;

    ang *= -SF_PI/180.;
    c = cosf(ang);
    s = sinf(ang);

    *x = xin * c - yin * s;
    *y = xin * s + yin * c;
}
