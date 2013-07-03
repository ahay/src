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

#include "tkirmig.h"

static int nt, ncmp, nh;
static float dt,t0,dcmp,cmpmin,cmpmax,h0,dh,dip;
static float **vrms, **mask;
static float *t_tmp;
static bool half,verb;

void tkirmig_init(int n1, float d1, float o1,
                 int n2, float d2, float o2,
                 int n3, float d3, float o3,
                 float dip1,
                 float **vel1,
                 float **mask1,
                 bool half1, 
                 bool verb1)
/*< Initialize >*/
{
    nt=n1;
    dt=d1;
    t0=o1;
    ncmp=n2;
    dcmp=d2;
    cmpmin=o2;
    cmpmax=cmpmin+(ncmp-1)*dcmp;
    nh=n3;
    dh=d3;
    h0=o3;
    dip = dip1;
    vrms=vel1;
    half=half1;
    verb=verb1;
    t_tmp = sf_floatalloc(nt);
    mask = mask1;
    sf_halfint_init(true,nt,1.-1./nt);
}

void tkirmig_close(void)
/*< Free allocated storage >*/
{
    free (t_tmp);
}

void tkirmig_lop(bool adj, bool add, int nm, int nd,
                float *modl, float *data)
/*< Apply >*/
{
     int ih,it,icmp,icdp,cdp1,cdp2,itime,itm,itp;
     float cmp,cdp,tau,tau2,z,minaper,maxaper;
     float disx,dish_plus,dish_plus2,dish_minus,dish_minus2;
     float vrms2,ts,tr,time,slope,tm,tp,amp,h,wt;

     sf_adjnull(adj,add,nm,nd,modl,data);
    
     for (ih=0; ih < nh; ih++) {
     h = h0+ih*dh;
     if (verb) sf_warning("ih/nh=%d/%d",ih+1,nh);
     if (!half) h = h/2.;
     for (icmp=0; icmp < ncmp; icmp++) {
         if (mask[ih][icmp]==0) continue;
         cmp=cmpmin+icmp*dcmp;
//         if (verb) sf_warning("icmp/ncmp=%d/%d,nt=%d,t0=%f,h=%f,v=%f",icmp+1,ncmp,nt,t0,h,vrms[icmp][10]);
         if (adj) {
            for (it=0; it < nt; it++) t_tmp[it]=data[ih*ncmp*nt+icmp*nt+it];
//            doubint(nt,t_tmp,false);
            sf_halfint (adj,t_tmp);
         }
         else {for (it=0; it < nt; it++) t_tmp[it]=0.;}

         for (it=0; it < nt; it++) {
             tau=(t0+it*dt)/2.0;
             tau2=tau*tau;
             z=vrms[icmp][it]*tau;
             minaper=cmp-z*tanf(dip);
             maxaper=cmp+z*tanf(dip);
             if (minaper <= cmpmin) minaper = cmpmin;
             if (maxaper >= cmpmax) maxaper = cmpmax;
             cdp1=(minaper-cmpmin)/dcmp;
             cdp2=(maxaper-cmpmin)/dcmp;

             for (icdp=cdp1; icdp < cdp2+1; icdp++) {
                 cdp=cmpmin+icdp*dcmp;
                 disx=cdp-cmp;
                 dish_plus=disx+h;
                 dish_plus2=dish_plus*dish_plus;
                 dish_minus=disx-h;
                 dish_minus2=dish_minus*dish_minus;
                 vrms2=vrms[icdp][it]*vrms[icdp][it];
                 ts=sqrtf(tau2+dish_plus2/vrms2);
                 tr=sqrtf(tau2+dish_minus2/vrms2);
                 time=ts+tr;

                 if (time <= t0+(nt-1)*dt && ts != 0. && tr !=0.) {
                    slope=dish_plus/vrms2/ts+dish_minus/vrms2/tr;
                    tm=time-fabsf(slope*dcmp*2.0);
                    tp=time+fabsf(slope*dcmp*2.0);
                    itime=0.5+(time-t0)/dt;
                    itm=0.5+(tm-t0)/dt;
//                    sf_warning("itm=%d,tm=%f,tp=%f,dcmp=%f,time=%f",itm,tm,tp,dcmp,time);
                    itp=0.5+(tp-t0)/dt;
                    wt=dt/(dt+(tp-tm));
                    wt*=wt;
//                    sf_warning("wt=%f",wt);
//                    amp=(dt/(dt+(tp-tm)))*(dt/(dt+(tp-tm)));
//                    amp=(dt/(dt+(tp-tm)))*(dt/(dt+(tp-tm)))*tau*sqrtf(ts/tr)/tr;
//                    amp=(dt/(dt+(tp-tm)))*(dt/(dt+(tp-tm)))*sqrtf(ts)*tau/tr/sqrtf(tr)/2./vrms[icdp][it]/vrms[icdp][it];
                    amp=wt*tau*(tr*tr+ts*ts)/ts/tr/sqrtf(ts)/sqrtf(tr)/vrms2;
                    if (itm>0&&itp<nt-2) {
                       spotw(adj,-amp,nt,t0,dt,tm,&modl[ih*ncmp*nt+icdp*nt+it],t_tmp);
                       spotw(adj,2*amp,nt,t0,dt,time,&modl[ih*ncmp*nt+icdp*nt+it],t_tmp);
                       spotw(adj,-amp,nt,t0,dt,tp,&modl[ih*ncmp*nt+icdp*nt+it],t_tmp);
                    }
                 }
              }
         }
         if (!adj) {
//            doubint2(nt,t_tmp,false);
            sf_halfint (adj,t_tmp);
            for (it=0;it<nt;it++) data[ih*ncmp*nt+icmp*nt+it]=t_tmp[it];
         }
     }
//     for (it=0; it <ncmp*nt; it++) sf_warning("modl[%d]=%f",it,modl[it]);
   }
}

void spotw(bool adj,float w,int nt,float t0,float dt,float t,float *val, float *vec)
/*< scaled linear intepolation >*/
{
     int it,itc;
     float tc,g;

     tc=(t-t0)/dt;
     itc=tc;
     it=1+itc;
     g=tc-itc;

     if (it>=0 && it < nt-1) {
        if (!adj) {
           vec[it]=vec[it]+(1.-g)*val[0]*w;
           vec[it+1]=vec[it+1]+g*val[0]*w;
        }
        else {
           val[0]=val[0]+((1-g)*vec[it]+g*vec[it+1])*w;
        }
      }
      
}

void doubint (int nx, float *xx, bool der)
/*< adjoin double integration >*/
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
        t += xx[i];
        xx[i] = t;
    }

    if (der) return;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
        t += xx[i];
        xx[i] = t;
    }
}

void doubint2 (int nx, float *xx, bool der)
/*< double integration >*/
{
    int i;
    float t;


    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
        t += xx[i];
        xx[i] = t;
    }

    if (der) return;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
        t += xx[i];
        xx[i] = t;
    }
}

