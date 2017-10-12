/* Kirchhoff migration kernel. */
/*
  Copyright (C) 2007 Dough McCowan
   
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
#include <math.h>

#include "kernel.h"

#define KP1MAX 64           /* max half length of triangle aa filter */

void kernel( float *datas           /* input data vector */,  
	     int npts_in            /* number of input data samples */,  
	     int npts_out           /* number of output image samples */,
             float xldist2          /* xline distance from this image position
				       to midpoint squared (m**2) */, 
	     float sldist2          /* sline distance from this image position
				       to midpoint squared (m**2) */, 
             float shot_xdist       /* surface x distance to shot (m) */, 
	     float shot_ydist       /* surface y distance to shot (m) */, 
	     float rcvr_xdist       /* surface x distance to receiver (m) */, 
	     float rcvr_ydist       /* surface y distance to receiver (m) */,
             int   *hits            /* input/output imaged switch 
				       +=1 image contribution made */,  
	     float fsint            /* input data sample interval (msec) */, 
	     float fosint           /* output image sample interval (msec) */,
	     float stretch_mute     /* stretch mute limit (msec) */,
             float *slineApMap      /* [ntimes] aperture^2 in sline */, 
	     float *xlineApMap      /* [ntimes] aperture^2 in xline */,
             float *vrmsGrid        /* [ntimes] rms velocities (m/sec) */, 
	     int ntimes             /* number of table times */, 
	     float dttab            /* vrms table sample interval (msec) */, 
	     float *scratch         /* [11*ntimes] travel time work arrays */,
             int amps               /* amplitude switch
				       =0 for no amplitude correction
				       =1 for spherical divergence correction
				       =2 for obliquity correction
				       =3 for obliq + sphdiv corrections
				       =4 for obliq - sphdiv corrections */, 
	     float *image           /* [npts_in] output image vector */, 
	     float *image_fold      /* [npts_in/8] output image fold vector */,
             int migrated_gathers   /* migrated gathers switch
				       =0 only compute image
				       =1 compute image and gathers */, 
	     float *gathers         /* [npts_in] output gathers vector */, 
	     float *gathers_fold    /* [npts_in/8] output gathers fold */,
             float trace_space      /* average trace spacing (m)
				       =0.001 for no anti-alias filtering */, 
	     float t0first          /* first output time to compute */)
/*< Engine for KTMIG3D migration program. This function utilizes the
  Lumley-Claerbout anti-alias filter >*/
/****************************************************************************
       Written for Oil of Olay  by  E. Coli                      4/12/00
   Changed obliquity from product to sum of cosines              7/23/03
*****************************************************************************/
{
    int   k, l, kk, ahit, mutes;
    register int kp1,ttx;
    
    float travel_time, prev_travel_time;
    float *shot_times, *rcvr_times, *kp1vals, *amplitudes;
    float *delslineaps, *delxlineaps, *delstimes, *delrtimes;
    float *delvrmss, *delscales, *delakp1s;
    float ft2, t0, basis, tnext;
    float anynum, vsq, image_value;
    float test, dt, fsinv;
    float shot_dist2, rcvr_dist2, aa, bb;
    float xlaper2, slaper2;
    float slineap, dslineap, xlineap, dxlineap, stime, dstime;
    float rtime, drtime, vrms, dvrms, scale, dscale, akp1, dakp1;
    float tmp2;

    /*  setup shop */

    ahit        = 0;
    shot_dist2  = shot_xdist*shot_xdist + shot_ydist*shot_ydist;
    rcvr_dist2  = rcvr_xdist*rcvr_xdist + rcvr_ydist*rcvr_ydist;
    fsinv       = 1.0f/fsint;
    dt          = 0.001f*fsint;
    tmp2        = 4000.0f * trace_space/dt;

    shot_times  = scratch;
    rcvr_times  = shot_times  + ntimes;
    kp1vals     = rcvr_times  + ntimes;
    amplitudes  = kp1vals     + ntimes;
    delslineaps = amplitudes  + ntimes;
    delxlineaps = delslineaps + ntimes;
    delstimes   = delxlineaps + ntimes;
    delrtimes   = delstimes   + ntimes;
    delvrmss    = delrtimes   + ntimes;
    delscales   = delvrmss    + ntimes;
    delakp1s    = delscales   + ntimes;

    /*  compute travel time and amplitude on velocity table grid  */

    for(k=0; k<ntimes; k++) {
        t0            = k*dttab;                      /* t0 in msec */
        vsq           = vrmsGrid[k]*vrmsGrid[k];      /* in m^^2/sec^^2 */
        anynum        = t0*0.0005f;                   /* t0/2 in sec */
        ft2           = anynum*anynum;                /* first term squared */
        shot_times[k] = 1000.0f*sqrt(ft2+shot_dist2/vsq); /* time to source */
        rcvr_times[k] = 1000.0f*sqrt(ft2+rcvr_dist2/vsq); /* time to receiv */

        aa = shot_xdist/shot_times[k] + rcvr_xdist/rcvr_times[k];
        bb = shot_ydist/shot_times[k] + rcvr_ydist/rcvr_times[k];
        kk = tmp2*sqrt(aa*aa + bb*bb)/vsq;
        if(kk<1) kk=1;
        if(kk>KP1MAX) kk=KP1MAX;
        kp1vals[k] = kk;

        anynum  = -1.0f/(float)(kk*kk);
        if(amps==1 || amps==3) 
	    anynum *= 1.0e-9f*vsq*(shot_times[k]+rcvr_times[k]);
/*        if(amps>1) anynum *= t0*t0/(shot_times[k]*rcvr_times[k]+0.1f);    */
        if(amps>1) anynum *= t0/(shot_times[k]+0.1f)+t0/(rcvr_times[k]+0.1f);
        if(amps==4) anynum /= 1.0e-9f*vsq*(shot_times[k]+rcvr_times[k]+0.1f);
        amplitudes[k] = anynum;
    }

    basis = fosint/dttab;
    for(k=0; k<ntimes-1; k++) {
        delslineaps[k] = (slineApMap[k+1]-slineApMap[k])*basis;
        delxlineaps[k] = (xlineApMap[k+1]-xlineApMap[k])*basis;
        delstimes[k]   = (shot_times[k+1]-shot_times[k])*basis;
        delrtimes[k]   = (rcvr_times[k+1]-rcvr_times[k])*basis;
        delvrmss[k]    = (vrmsGrid[k+1]-vrmsGrid[k])*basis;
        delscales[k]   = (amplitudes[k+1]-amplitudes[k])*basis;
        delakp1s[k]    = (kp1vals[k+1]-kp1vals[k])*basis;
    }
    delslineaps[ntimes-1] = 0;
    delxlineaps[ntimes-1] = 0;
    delstimes[ntimes-1]   = 0;
    delrtimes[ntimes-1]   = 0;
    delvrmss[ntimes-1]    = 0;
    delscales[ntimes-1]   = 0;
    delakp1s[ntimes-1]    = 0;

    /*  loop over the output time samples  */

    travel_time = prev_travel_time = 0.0f;
    mutes = 0;
    k=0; t0=0.0f; tnext=dttab;
    slineap = slineApMap[k]; dslineap = delslineaps[k];
    xlineap = xlineApMap[k]; dxlineap = delxlineaps[k];
    stime   = shot_times[k]; dstime   = delstimes[k];
    rtime   = rcvr_times[k]; drtime   = delrtimes[k];
    vrms    = vrmsGrid[k];   dvrms    = delvrmss[k];
    scale   = amplitudes[k]; dscale   = delscales[k];
    akp1    = kp1vals[k];    dakp1    = delakp1s[k];
    for(l=0; l<npts_out; l++, t0+=fosint) {

        if(k<ntimes-1 && t0>=tnext) {
           k++; tnext+=dttab;
           slineap = slineApMap[k]; dslineap = delslineaps[k];
           xlineap = xlineApMap[k]; dxlineap = delxlineaps[k];
           stime   = shot_times[k]; dstime   = delstimes[k];
           rtime   = rcvr_times[k]; drtime   = delrtimes[k];
           vrms    = vrmsGrid[k];   dvrms    = delvrmss[k];
           scale   = amplitudes[k]; dscale   = delscales[k];
           akp1    = kp1vals[k];    dakp1    = delakp1s[k];
        }

	/* check the mute */

        if(t0<t0first) goto again;

	/*  check aperture conditions, aperture vectors are squares  */

        slaper2 = slineap;
        xlaper2 = xlineap;
        if(sldist2*slaper2 + xldist2*xlaper2 > xlaper2*slaper2) goto again;

	/*  mute if necessary and set previous travel time  */

        travel_time = stime + rtime;          /* in msec */
        test = travel_time - prev_travel_time;
        if(test < stretch_mute ) goto again;
        if(++mutes<3) goto again;             /* skip first 2 values */

        /*  get the travel time index = ttx  */

        ttx = travel_time*fsinv+0.5f;
        if(ttx > npts_in-1) goto quit;

	/*  migrate the data  */

        kp1 = akp1;
        image_value = scale*(2.0f*datas[ttx]-datas[ttx+kp1]-datas[ttx-kp1]);
        if(image_value == 0.0f) goto again;     /* skip zero values */
        image[l] += image_value;
        image_fold[l>>3] += 1.0f;
        ahit=1;

	/*  get the migrated gathers contribution  */

        if( migrated_gathers ) {
          gathers[l] += image_value;
          gathers_fold[l>>3] +=1.0f;
        }

again:  prev_travel_time = travel_time;

        slineap += dslineap; xlineap += dxlineap;
        stime   += dstime;   rtime   += drtime;
        vrms    += dvrms;    scale   += dscale;
        akp1    += dakp1;
	
    }           /* end the l loop */


quit: if(ahit) (*hits)++;
      return;
}
