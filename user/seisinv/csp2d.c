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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "csp2d.h"

static int ompnth, ompchunk;
static int nhe, nxs, nh, nxm, nt;
static float dhe, he0, dxs, xs0, dh, h0, dxm, xm0, dt, t0, v, v2, apt;
static bool weight, linear, verb, half;

void csp2d_init(int ompnth1, int ompchunk1, /* openmp info */
                int nhe_in, float dhe_in, float he0_in,
                int nxs_in, float dxs_in, float xs0_in,
                int nh_in, float dh_in, float h0_in,
                int nxm_in, float dxm_in, float xm0_in,
                int nt_in, float dt_in, float t0_in,
                float apt_in, float v_in, bool weight_in,
                bool linear_in, bool half_in, bool verb_in)
/*< Initialize >*/
{
    ompnth = ompnth1;
    ompchunk = ompchunk1;

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

    apt = apt_in;

    v = v_in;
    v2 = v*v;

    weight = weight_in;
    linear=linear_in;
    half = half_in;
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
    int ihe, im, id, it, ixs, ixm, ih;
    float t, t2, xs, xm, x, x2, h, h2, he, he_max, fx, gx;
    float w;

    sf_adjnull(adj,add,nm,nd,modl,data);

    /* Time loop should be parallel*/
#ifdef _OPENMP
#pragma omp parallel for            \
    private(ihe, im, id, it, ixs, ixm, ih, \
      t, t2, xs, xm, x, x2, h, h2, he, he_max, fx, gx, w)
#endif
    for (it = 1; it < nt; it++) {

#ifdef _OPENMP      
#pragma omp critical
#endif
        if(verb) sf_warning("it/nt=%d/%d;",it+1,nt);
        t = t0 + it*dt;
        t2 = t*t;

        /* CSP position loop */
        for (ixs = 0; ixs < nxs; ixs++) {
            xs = xs0 + ixs*dxs;

            /* CMP position loop */
            for (ixm = 0; ixm < nxm; ixm++) {
                xm = xm0 + ixm*dxm;
                x = fabsf(xm - xs);
                x2 = x*x;

                if (v*t < 2.*x || x > apt) continue;

                /* offset loop */
                for (ih = 0; ih < nh; ih++) {
                    h = h0 + ih*dh;
                    if (!half) h = h/2.;
                    h2 = h*h;
                    he_max = sqrt(x2+h2);

                    /* EO equation */
                    he = sqrt(x2 + h2 * (1. - 4.*x2/t2/v2));
                    if (he > he_max) continue;

                    if (h < 0) he = -he;

                    if (weight) {
                       w = 1.-x2/(he+1.e-12)/(he+1.e-12);
                    } else {
                       w = 1.;
                    }

                    /* if (weight) w = 1.-fabsf(x/(he+1.e-12)); */

                    he = (he-he0)/dhe+0.5;
                    ihe = floorf(he);

                    if (ihe < 0 || ihe >= nhe) continue;

                    id = ixm*nt*nh+ih*nt+it;
                    im = ixs*nt*nhe+ihe*nt+it;

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
#ifdef _OPENMP      
#pragma omp critical
#endif
    sf_warning(".");
}
