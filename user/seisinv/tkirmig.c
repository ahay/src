/* Kirchhoff prestack time modeling/migration with antialiase */
/*
  Copyright (C) 2012 China University of Petroleum
  
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

#include "tkirmig.h"
#include "doubint.h"

static int ompnth, ompchunk;
static int nt, ncmp, ncdp, nh;
static float dt,t0,dcmp,cmp0,cdp0,dcdp,h0,dh,rho,apt,aal;
static float **vrms, *off, **mask;
static float **img, **trace;
static bool verb, amp;

void tkirmig_init(int ompnth1, int ompchunk1, /* openmp info */
                 int n1, float d1, float o1, /* time axis */
                 int n2, float d2, float o2,  /* cmp axis */
                 int n2_im, float d2_im, float o2_im,  /* cdp axis */
                 int n3, float d3, float o3,  /* off axis */
                 float apt1,                   /* aperture */
                 float rho1,
                 float aal1,
                 float **vel1,
                 float *off1,
                 float **mask1,
                 bool amp1,
                 bool verb1)
/*< Initialize >*/
{
    ompnth = ompnth1;
    ompchunk = ompchunk1;
    nt=n1;
    dt=d1;
    t0=o1;

    ncmp=n2;
    dcmp=d2;
    cmp0=o2;

    ncdp=n2_im;
    dcdp=d2_im;
    cdp0=o2_im;

    nh=n3;
    dh=d3;
    h0=o3;

    apt = apt1;
    aal = aal1;
    rho = rho1;

    vrms = vel1;
    off = off1;
    mask = mask1;

    amp = amp1;
    verb = verb1;

    trace = sf_floatalloc2(nt,ompnth);
    img = sf_floatalloc2(nt,ompnth);

    sf_doubint_init(nt);
    sf_halfint_init (true, nt, rho);
}

void tkirmig_close(void)
/*< Free allocated storage >*/
{
    free (*img), free(img);
    free (*trace), free(trace);
    sf_doubint_close();
}

void tkirmig_lop(bool adj, bool add, int nm, int nd,
                float *modl, float *data)
/*< Apply >*/
{
     int ih,it,icmp,icdp,itm,itp,im;
     float cmp,cdp,tau,tau2;
     float disx,dish_plus,dish_plus2,dish_minus,dish_minus2;
     float vrms2,ts,tr,time,slope,tm,tp,h,wt;
     int ompith=0;

     sf_adjnull(adj,add,nm,nd,modl,data);
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static,ompchunk) \
    private(ompith,ih,it,icmp,icdp,itm,itp,im, \
      cmp,cdp,tau,tau2, \
      disx,dish_plus,dish_plus2,dish_minus,dish_minus2, \
      vrms2,ts,tr,time,slope,tm,tp,h,wt) \
    shared(adj, add, nm, nd, modl, data)
#endif
     for (ih=0; ih < nh; ih++) {
#ifdef _OPENMP 
         ompith=omp_get_thread_num();
#pragma omp critical
#endif
         if (verb) sf_warning("ith=%d; offset %d of %d;",ompith, ih+1,nh);

         for (icmp=0; icmp < ncmp; icmp++) {

             if (mask[ih][icmp]==0) continue;

             cmp=cmp0+icmp*dcmp;

             h = fabsf(off[ih*ncmp+icmp]);

             if (adj) {
                for (it=0; it < nt; it++) {
                    trace[ompith][it]=data[ih*ncmp*nt+icmp*nt+it];
                    img[ompith][it] = 0.;
                }
                sf_doubint_lop(adj, false, nt, nt, img[ompith], trace[ompith]);
             } else {
                for (it=0; it < nt; it++) {
                    img[ompith][it] = 0.;
                    trace[ompith][it] = 0.;
                }
             }

           for (icdp=0; icdp < ncdp; icdp++) {
 
              cdp = cdp0+icdp*dcdp;

              disx=cdp-cmp;

              if (fabsf(disx) > dcmp*apt) continue;

              dish_plus=fabsf(disx+h);
              dish_plus2=dish_plus*dish_plus;
              dish_minus=fabsf(disx-h);
              dish_minus2=dish_minus*dish_minus;

              for (it=0; it < nt; it++) {
                  im = ih*ncdp*nt+icdp*nt+it;
                  tau = 0.5*(t0 + it*dt);
                  tau2 = tau*tau;

                  vrms2=vrms[icdp][it]*vrms[icdp][it];
                  ts=sqrtf(tau2+dish_plus2/vrms2);
                  tr=sqrtf(tau2+dish_minus2/vrms2);
                  time=ts+tr;

                  if (time <= t0+(nt-1)*dt && ts != 0. && tr !=0.) {

                     slope=dish_plus/vrms2/(ts+dt)+dish_minus/vrms2/(tr+dt);

                     tm=time-fabsf(slope*dcmp*aal)-dt;
                     tp=time+fabsf(slope*dcmp*aal)+dt;
                     itm=floorf((tm-t0)/dt);
                     itp=floorf((tp-t0)/dt);
                     wt=dt/(dt+tp-tm);
                     wt*=wt;

                     if (amp) wt*=tau*(tr*tr+ts*ts)/ts/tr/sqrtf(ts)/sqrtf(tr)/vrms2;
                     //amp = wt;

                     if (itm>=0&&itp<nt-1) {
                        spotw(adj,-wt,nt,t0,dt,tm,&modl[im],img[ompith]);
                        spotw(adj,2*wt,nt,t0,dt,time,&modl[im],img[ompith]);
                        spotw(adj,-wt,nt,t0,dt,tp,&modl[im],img[ompith]);
                     }
                  }

               } /* time loop */

            } /* cdp loop */

            if (!adj) {
               sf_doubint_lop(adj, false, nt, nt, img[ompith], trace[ompith]);
               for (it=0; it<nt; it++) data[ih*ncmp*nt+icmp*nt+it] = trace[ompith][it];
            }

        } /* cmp loop */
    } /* offset loop */

    if (verb) sf_warning(".");
}

void spotw(bool adj,float w,int nt,float t0,float dt,float t,float *val, float *vec)
/*< scaled linear intepolation >*/
{
     int it;
     float tc,g;

     //sf_adjnull(adj,false,nt,nt,val,vec);

     tc=(t-t0)/dt;
     it=floorf(tc);
     g=tc-it;

     if (it>=0 && it < nt-1) {
        if (!adj) {
           vec[it] += (1.-g)*val[0]*w;
           vec[it+1] += g*val[0]*w;
        }
        else {
           val[0] += ((1-g)*vec[it]+g*vec[it+1])*w;
        }
      }
}
