/* Test migration kernel. */
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
#include <time.h>

#include <rsf.h>

#include "kernel.h"

#define EXPAND 4            /* the expansion factor for the input data */
#define NPTS 1000
#define NTIMES 42
#define KP1MAX 64           /* max half length of triangle aa filter */

int main(int argc, char **argv)
{
  int   i;
  double cpu;
  clock_t stime, etime;
  float aper;
  float dtdata=4.0, str_mut=0.35, dtout=4.0, dttab=96.0;
  float sldist2=160000.0, xldist2=640000.0, dxbar=50.0, t0first=500.0;
  float shot_xd=-400.0, shot_yd=1000.0, rcvr_xd=2000.0, rcvr_yd=-1800.0;
  int   in_nsamps=NPTS, *hits, *amps;
  int   out_nsamps=NPTS, ncalls;
  int   migrated_gathers=1;
  float temp1[EXPAND*NPTS+2*KP1MAX], temp2[EXPAND*NPTS+2*KP1MAX];
  float *data1, *data2;
  float saper[NTIMES], xaper[NTIMES], vrms[NTIMES], scratch[11*NTIMES];
  float image[NPTS], image_fold[NPTS/8+10], gathers[NPTS], aperture[NTIMES];
  float gathers_fold[NPTS/8+10];

  sf_init(argc,argv);
  if (!sf_getint("ncalls",&ncalls)) ncalls=5000;
  /* number of calls */

  for(i=0; i<EXPAND*NPTS+2*KP1MAX; i++) {
      temp1[i] = temp2[i] = 3.0f;
  }
  data1 = temp1 + KP1MAX; 
  data2 = temp2 + KP1MAX;
  for(i=0; i<EXPAND*NPTS; i++) {
      data1[i] = data2[i] = 3.3f;
  }
  for(i=0; i<NTIMES; i++) {
      if(i<=10) aper = i*500;
      else      aper = 5000.0;
      aperture[i] = aper;
      saper[i] = xaper[i] = aper*aper;
      vrms[i]  = 5000.0 + 38.7*i;
  }
  for(i=0; i<out_nsamps; i++) {
      image[i] = gathers[i] = 0.0;
  }
  for(i=0; i<NPTS/8+10; i++) {
      image_fold[i] = gathers_fold[i] = 0.0;
  }

  hits = sf_intalloc(ncalls);
  amps = sf_intalloc(ncalls);

  for(i=0; i<ncalls; i++) {
      hits[i] = 0;
      amps[i] = 2;
  }

  stime = clock();

  for(i=0; i<ncalls; i++) {
      kernel( data1,  EXPAND*in_nsamps,  out_nsamps,
	      xldist2, sldist2,
	      shot_xd, shot_yd, rcvr_xd, rcvr_yd,
	      &hits[i],  0.25f*dtdata, dtout,
	      str_mut, saper, xaper,
	      vrms, NTIMES, dttab, scratch,
	      amps[i], image, image_fold,
	      migrated_gathers, gathers, gathers_fold,
	      dxbar, t0first);
  }

  etime = clock();

  cpu = ((double) (etime-stime))/ CLOCKS_PER_SEC;
  printf("Kernel CPU time = %f s for %d calls\n", cpu, ncalls);

  exit(0);
}
