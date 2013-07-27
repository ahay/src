/* 2-D common scattering-point gathers mapping operator and it's adjoint */
/*
  Copyright (C) 2012 China University of Petroleum (East China)
  
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
/*^*/

#include "csp2d.h"

static int nhe, nxs, nh, nxm, nt;
static float dhe, he0, dxs, xs0, dh, h0, dxm, xm0, dt, t0, v, v2;
static bool weight, linear, verb;

void csp2d_init(int nhe_in, float dhe_in, float he0_in,
                int nxs_in, float dxs_in, float xs0_in,
                int nh_in, float dh_in, float h0_in,
                int nxm_in, float dxm_in, float xm0_in,
                int nt_in, float dt_in, float t0_in,
                float v_in, bool weight_in,
                bool linear_in, bool verb_in)
/*< Initialize >*/
{
    nhe = nhe_in;
    dhe = dhe_in;
    he0 = he0_in;
    nxs = nxs_in;
    dxs = dxs_in;
    xs0 = xs0_in;
 
    nh = nh_in;
    dh = dh_in;
    h0 = h0_in;
    nxm = nxm_in;
    dxm = dxm_in;
    xm0 = xm0_in;

    nt = nt_in;
    dt = dt_in;
    t0 = t0_in;

    v = v_in;
    v2 = v*v;

    weight = weight_in;
    linear=linear_in;
    verb = verb_in;
}

void csp_close(void)
/*< Free allocated storage >*/
{
}

void csp2d_lop(bool adj, bool add, int nm, int nd,
                float *modl, float *data)
/*< Apply >*/
{
    int ihe, im, id;
    float t, t2, xs, xm, x, x2, h, h2, he, fx, gx;
    float w = 1.;

    sf_adjnull(adj,add,nm,nd,modl,data);

    /* Time loop should be parallel*/
    for (int it = 1; it < nt; it++) {
        if(verb) sf_warning("it/nt=%d/%d;",it+1,nt);
        t = t0 + it*dt;

        t2 = t*t;

        /* CSP position loop */
        for (int ixs = 0; ixs < nxs; ixs++) {
            xs = xs0 + ixs*dxs;

            /* CMP position loop */
            for (int ixm = 0; ixm < nxm; ixm++) {
                xm = xm0 + ixm*dxm;
                x = xm - xs;
                x2 = x*x;

                if (t < 2.*x/v) continue;

                /* offset loop */
                for (int ih = 0; ih < nh; ih++) {
                    h = h0 + ih*dh;
                    h2 = h*h;

                    /* EO equation */
                    he = sqrt(x2 + h2 * (1. - 4.*x2/t2/v2));
                    if (h < 0) he = -he;

                    if (weight) w = 1.-x2/(he+1.e-6)/(he+1.e-06);

                    he = (he-he0)/dhe;
                    ihe = floorf(he);

                    id = ih+ixm*nh+it*nxm*nh;
                    im = ihe+ixs*nhe+it*nxs*nhe;

                    if (!linear) {
                    if (im>=0&&im<nm) {
                       if (adj) {
                          modl[im] += w*data[id];
                       } else {
                          data[id] += w*modl[im];
                       }
                    }

                    } else {
                    fx = he - ihe;
                    gx = 1. - fx;

                    if (im<0||im>nm-2) continue;
                    if (adj) {
                       modl[im] += gx * w * data[id];
                       modl[im+1] += fx * w * data[id];
                    } else {
                       data[id] += gx * w * modl[im] + fx * w * modl[im+1];
                    }}
                }
            }
        }
    }
    sf_warning(".");
}
