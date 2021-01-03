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

#include "adjgradient2d.h"

int main(int argc, char* argv[])
{
  int nw,nx,nz,nh;
  float ox,oz,dx,dz,ow,dw,oh,dh;
  
  /* I/O files */
  sf_file Fxig = NULL; /* Input XIG */
  sf_file Fvel = NULL; /* velocity file */
  sf_file Frwf = NULL; /* receiver wavefield @ nz-1 */
  sf_file Fswf = NULL; /* source   wavefield @ nz-1 */ 
  sf_file Fgrd = NULL; /* Output gradient */

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);
  
  /* cube axes */
  sf_axis ax,aw,az,ah;
  
  /* Set up file objects */
  Fvel = sf_input("in");  /* input velocity   */
  Fxig = sf_input("xig"); /* input penalized image */
  Fswf = sf_input("swf"); sf_settype(Fswf,SF_COMPLEX);/* output wavefield */
  Frwf = sf_input("rwf"); sf_settype(Frwf,SF_COMPLEX);/* output wavefield */
  Fgrd = sf_output("out"); /* Output gradient */

   /* Read in parameter */ 
  int nxtap;
  bool verbose;
  if(! sf_getint("nxtap",&nxtap)) nxtap=40;  /* TAPER size */
  if(! sf_getbool("verbose",&verbose)) verbose=false; /* VERBOSITY flag */

  /* Get axes information */
  ax = sf_iaxa(Frwf,1); if (verbose) sf_raxa(ax);
  az = sf_iaxa(Fvel,2); if (verbose) sf_raxa(az);
  aw = sf_iaxa(Frwf,2); if (verbose) sf_raxa(aw);
  ah = sf_iaxa(Fxig,1); if (verbose) sf_raxa(ah);
  
  /* Pull out individual values */
  nx = sf_n(ax); ox = sf_o(ax); dx = sf_d(ax);
  nw = sf_n(aw); ow = sf_o(aw); dw = sf_d(aw);
  nz = sf_n(az); oz = sf_o(az); dz = sf_d(az);
  nh = sf_n(ah); oh = sf_o(ah); dh = sf_d(ah);

  /* I/O arrays */
  float **vel=NULL;  /* v array */
  float ***xig=NULL; /* Input xig array */
  sf_complex **swf=NULL; /* input SWF @ NZ-1 array */
  sf_complex **rwf=NULL; /* input RSF @ NZ-1 array */
  float **grd=NULL; /* Output gradient */

  /* Allocation */
  vel =   sf_floatalloc2(nx,nz);
  xig =   sf_floatalloc3(nh,nx,nz);
  swf = sf_complexalloc2(nx,nw);
  rwf = sf_complexalloc2(nx,nw);
  grd =   sf_floatalloc2(nx,nz);

  /* Read in velocity */
  sf_warning("Read in Velocity Model");
  sf_floatread( vel[0],nx*nz,Fvel );	
 
  /* Read in XIG */
  sf_warning("Read in XIG");
  sf_floatread( xig[0][0],nh*nx*nz,Fxig );	

  /* Read in SWF and RWF */
  sf_warning("Reading in SWF and RWF");
  sf_complexread(swf[0],nx*nw,Fswf);
  sf_complexread(rwf[0],nx*nw,Frwf);

  /* Initialize the module */
  adjgradient2d_init(nx,nz,nw,nh,ox,oz,ow,oh,dx,dz,dw,dh,nxtap,vel);
      
  /* IMAGE INDIVIDUAL SHOTS */
  adjgradient2d_wemig(vel,xig,swf,rwf,grd);	

  /* Write out SWF and RWF */
  sf_warning("Writing out gradient");
  sf_floatwrite(grd[0],nx*nz,Fgrd);
    
  adjgradient2d_close();
  
  exit(0);
}
