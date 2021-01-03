/* Gradient adjoint-state calculation for image-domain WET */
/*
  Copyright (C) 2012 The University of Western Australia
  
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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "adjgradient2d_coupled.h"

int main(int argc, char* argv[])
{
  int nw,nx,nz,nh;
  float ox,oz,dx,dz,ow,dw,oh,dh;
  
  /* I/O files */
  sf_file Fxig1 = NULL; /* Input XIG */
  sf_file Fvel1 = NULL; /* velocity file */
  sf_file Fur1 = NULL; /* receiver wavefield @ nz-1 */
  sf_file Fus1 = NULL; /* source   wavefield @ nz-1 */ 
  sf_file Fgr1 = NULL; /* Output gradient */

  sf_file Fxig2 = NULL; /* Input XIG */
  sf_file Fvel2 = NULL; /* velocity file */
  sf_file Fur2 = NULL; /* receiver wavefield @ nz-1 */
  sf_file Fus2 = NULL; /* source   wavefield @ nz-1 */ 
  sf_file Fgr2 = NULL; /* Output gradient */
  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);
  
  /* cube axes */
  sf_axis ax,aw,az,ah;
  
  /* Set up file objects */
  Fvel1 = sf_input("in");  /* input velocity   */
  Fxig1 = sf_input("xig1"); /* input penalized image */
  Fus1 = sf_input("us1"); sf_settype(Fus1,SF_COMPLEX);/* output wavefield */
  Fur1 = sf_input("ur1"); sf_settype(Fur1,SF_COMPLEX);/* output wavefield */
  Fgr1 = sf_output("out"); /* Output gradient */
  
  Fvel2 = sf_input("vel2"); /* input velocity */
  Fxig2 = sf_input("xig2"); /* input penalized image */
  Fus2 = sf_input("us2"); sf_settype(Fus2,SF_COMPLEX);/* output wavefield */
  Fur2 = sf_input("ur2"); sf_settype(Fur2,SF_COMPLEX);/* output wavefield */
  Fgr2 = sf_output("gr2"); /* Output gradient */

   /* Read in parameter */ 
  int nxtap;
  bool verbose;
  if(! sf_getint("nxtap",&nxtap)) nxtap=40;  /* TAPER size */
  if(! sf_getbool("verbose",&verbose)) verbose=false; /* VERBOSITY flag */

  /* Get axes information */
  ax = sf_iaxa(Fur1, 1); if (verbose) sf_raxa(ax);
  aw = sf_iaxa(Fur1, 2); if (verbose) sf_raxa(aw);
  az = sf_iaxa(Fvel1,2); if (verbose) sf_raxa(az);
  ah = sf_iaxa(Fxig1,1); if (verbose) sf_raxa(ah);
  
  /* Pull out individual values */
  nx = sf_n(ax); ox = sf_o(ax); dx = sf_d(ax);
  nw = sf_n(aw); ow = sf_o(aw); dw = sf_d(aw);
  nz = sf_n(az); oz = sf_o(az); dz = sf_d(az);
  nh = sf_n(ah); oh = sf_o(ah); dh = sf_d(ah);

  /* I/O arrays */
  float **vel1=NULL;  /* v array */
  float **vel2=NULL;  /* v array */
  float ***xig1=NULL; /* Input xig array */
  float ***xig2=NULL; /* Input xig array */
  sf_complex **us1=NULL; /* input SWF @ NZ-1 array */
  sf_complex **us2=NULL; /* input SWF @ NZ-1 array */
  sf_complex **ur1=NULL; /* input RSF @ NZ-1 array */
  sf_complex **ur2=NULL; /* input RSF @ NZ-1 array */
  float **gr1=NULL; /* Output gradient */
  float **gr2=NULL; /* Output gradient */

  /* BASELINE Allocation */
  vel1 = sf_floatalloc2(nx,nz);
  xig1 = sf_floatalloc3(nh,nx,nz);
  us1 =  sf_complexalloc2(nx,nw);
  ur1 =  sf_complexalloc2(nx,nw);
  gr1 =  sf_floatalloc2(nx,nz);

  /* MONITOR Allocation */
  vel2 = sf_floatalloc2(nx,nz);
  xig2 = sf_floatalloc3(nh,nx,nz);
  us2 =  sf_complexalloc2(nx,nw);
  ur2 =  sf_complexalloc2(nx,nw);
  gr2 =  sf_floatalloc2(nx,nz);

  /* Read in velocity */
  sf_warning("Read in Velocity Model");
  sf_floatread( vel1[0],nx*nz,Fvel1 );	
  sf_floatread( vel2[0],nx*nz,Fvel2 );	
 
  /* Read in XIG */
  sf_warning("Read in XIG");
  sf_floatread( xig1[0][0],nh*nx*nz,Fxig1);	
  sf_floatread( xig2[0][0],nh*nx*nz,Fxig2);	

  /* Read in SWF and RWF */
  sf_warning("Reading in SWF and RWF");
  sf_complexread(us1[0],nx*nw,Fus1);
  sf_complexread(ur1[0],nx*nw,Fur1);
  sf_complexread(us2[0],nx*nw,Fus2);
  sf_complexread(ur2[0],nx*nw,Fur2);

  /* Initialize the module */
  adjgradient2d_coupled_init(nx,nz,nw,nh,ox,oz,ow,oh,dx,dz,dw,dh,nxtap,vel1,vel2);
      
  /* IMAGE INDIVIDUAL SHOTS */
  adjgradient2d_coupled_wemig(vel1,xig1,us1,ur1,gr1,vel2,xig2,us2,ur2,gr2);	

  /* Write out SWF and RWF */
  sf_warning("Writing out gradient");
  sf_floatwrite(gr1[0],nx*nz,Fgr1);
  sf_floatwrite(gr2[0],nx*nz,Fgr2);
    
  adjgradient2d_coupled_close();
  
  exit(0);
}
