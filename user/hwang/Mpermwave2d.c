/* Wavefield Extrapolation for 2D PERM */
/*
   Copyright (C) 2013 King Abduallah University of Science and Technology

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
  
  /*
   * :TODO: 
   *  1.Add option to specify shot position. 
   *    For now the shot position is on every grid point at surface
   * :WARNING: 
   *  1.Refer to parallelized version for large models
   * ========================================================================== 
   * I/O File   :                         Description
   * --------------------------------------------------------------------------
   *  in/out    :    [nh][nx][nz]         reflectivity 
   *  out/in    :    [nt][nh][nx]         data (Time reversed if migration)
   *  left      :    [nzx][m2]            lowrank factorized left matrix
   *  right     :    [m2][nk]             lowrank factorized right matrix
   *  wavelet   :    [nt]                 source wavelet (for modeling)
   * ========================================================================== 
   */

#include <rsf.h>
#include "fft3.h"

int main(int argc, char* argv[])
{
  bool mig;
  int it, nt, ix, nx, ih, nh, iz, nz, nx2, nh2, nz2, nxy, nzx, nzx2;
  int im, i, j, m2, it1, it2, its, ik, n2, nk;
  float dt, dx, dh, dz, x0, h0;
  float *curr, *prev, ***img, **dat, **lft, **rht, *wave;
  sf_complex *cwave, *cwavem;
  sf_file data, image, left, right;
  sf_file wavelet_file=NULL;
  float *wavelet=NULL;
  float *ptrtmp=NULL;
  int imageit;

  sf_init(argc,argv);

  if (!sf_getbool("mig",&mig)) mig=false;/* if n, modeling; if y, migration */


  if (mig) { /* migration */
    data = sf_input("in");
    image = sf_output("out");

    if (!sf_histint(data,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(data,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(data,"o1",&x0)) x0=0.; 

    if (!sf_histint(data,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histfloat(data,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(data,"o2",&h0)) h0=0.; 

    if (!sf_histint(data,"n3",&nt)) sf_error("No n3= in input");
    if (!sf_histfloat(data,"d3",&dt)) sf_error("No d3= in input");

    if (!sf_getint("nz",&nz)) sf_error("Need nz=");/* depth axis */
    if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
    if (!sf_getint("imageit",&imageit)) imageit=0;/* time for extracting image */

    sf_putint(image,"n1",nz);
    sf_putfloat(image,"d1",dz);
    sf_putfloat(image,"o1",0.);
    sf_putstring(image,"label1","Depth");

    sf_putint(image,"n2",nx);
    sf_putfloat(image,"d2",dx);
    sf_putfloat(image,"o2",x0);
    sf_putstring(image,"label2","Midpoint");


    sf_putint(image,"n3",nh);
    sf_putfloat(image,"d3",dh);
    sf_putfloat(image,"o3",h0);
    sf_putstring(image,"label3","Half-Offset");
  } else { /* modeling */
    image = sf_input("in");
    data = sf_output("out");
    wavelet_file = sf_input("wavelet");

    if (!sf_histint(image,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(image,"d1",&dz)) sf_error("No d1= in input");

    if (!sf_histint(image,"n2",&nx))  sf_error("No n2= in input");
    if (!sf_histfloat(image,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(image,"o2",&x0)) x0=0.; 	

    if (!sf_histint(image,"n3",&nh)) sf_error("No n3= in input");
    if (!sf_histfloat(image,"d3",&dh)) sf_error("No d3= in input");
    if (!sf_histfloat(image,"o3",&h0)) h0=0.; 

    if (!sf_histint(wavelet_file,"n1",&nt)) sf_error("No n1= in wavelet");
    if (!sf_histfloat(wavelet_file,"d1",&dt)) sf_error("No d1= in wavelet");

    sf_putint(data,"n1",nx);
    sf_putfloat(data,"d1",dx);
    sf_putfloat(data,"o1",x0);
    sf_putstring(data,"label1","Midpoint");

    sf_putint(data,"n2",nh);
    sf_putfloat(data,"d2",dh);
    sf_putfloat(data,"o2",h0);
    sf_putstring(data,"label2","Half-Offset");

    sf_putint(data,"n3",nt);
    sf_putfloat(data,"d3",dt);
    sf_putfloat(data,"o3",0.);
    sf_putstring(data,"label3","Time");
    sf_putstring(data,"unit3","s");
  }

  nk = fft3_init(false,1,nh,nx,nz,&nh2,&nx2,&nz2);

  nxy = nx*nh;
  nzx = nz*nx*nh;
  nzx2 = nz2*nx2*nh2;

  img = sf_floatalloc3(nz,nx,nh);
  dat = sf_floatalloc2(nx,nh);

  /* propagator matrices */
  left = sf_input("left");
  right = sf_input("right");

  if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
  if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");

  if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
  if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);

  lft = sf_floatalloc2(nzx,m2);
  rht = sf_floatalloc2(m2,nk);

  sf_floatread(lft[0],nzx*m2,left);
  sf_floatread(rht[0],m2*nk,right);

  curr = sf_floatalloc(nzx2);
  prev = sf_floatalloc(nzx2);

  cwave  = sf_complexalloc(nk);
  cwavem = sf_complexalloc(nk);
  wave = sf_floatalloc(nzx2);
  if (!mig)
    wavelet = sf_floatalloc(nt);

  ifft3_allocate(cwavem);

  memset(prev,0,sizeof(float)*nzx2);
  memset(curr,0,sizeof(float)*nzx2);

  if (mig) { /* migration */
    /* step backward in time */
    it1 = nt-1;
    it2 = -1;
    its = -1;	
  } else { /* modeling */
    sf_floatread(wavelet,nt,wavelet_file);
    sf_floatread(img[0][0],nzx,image);

    /* step forward in time */
    it1 = 0;
    it2 = nt;
    its = +1;
  }

  /* time stepping */
  for (it=it1; it!=it2; it+=its) {
    sf_warning("it=%d;",it);

    if (mig) /* migration <- read data */
      sf_floatread(dat[0],nxy,data);
    else 
      memset(dat[0],0,sizeof(float)*nxy);

    if (!mig) { /* modeling -> inject source */
      for (ih=0; ih<nh; ih++) 
        for (ix=0; ix<nx; ix++) 
          for (iz=0; iz<nz; iz++) 
            curr[ih+nh2*(ix+iz*nx2)] += img[ih][ix][iz]*wavelet[it];
    }

    for (ih=0; ih < nh; ih++) {
      for (ix=0; ix < nx; ix++) {
        if (mig) 
          curr[ih+nh2*ix] += dat[ih][ix];
        else 
          dat[ih][ix] = curr[ih+nh2*ix];
      }
    }

    fft3(curr,cwave);

    for (j=0; j<nzx2; j++)
      prev[j] *= -1;

    /* matrix multiplication */
    for (im=0; im<m2; im++) {
      for (ik=0; ik<nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
        cwavem[ik] = cwave[ik]*rht[ik][im];
#else
        cwavem[ik] = sf_crmul(cwave[ik],rht[ik][im]);
#endif
      }

      ifft3(wave,cwavem);

      for (ih=0; ih<nh; ih++) {
        for (ix=0; ix<nx; ix++) {
          for (iz=0; iz<nz; iz++) {
            i = ih+nh*(ix+iz*nx);   /* original grid */
            j = ih+nh2*(ix+iz*nx2); /* padded grid */
            prev[j] += lft[im][i]*wave[j];
          }
        }
      }
    } /* End of m2 loop */

    /* loop over pointers */
    ptrtmp = prev;
    prev = curr;
    curr = ptrtmp;

    if (!mig) { /* modeling -> write out data */
      sf_floatwrite(dat[0],nxy,data);
    } else {
      if (it==imageit ) {  /* migration -> write out image */
        /* transpose */
        for (ih=0; ih<nh; ih++) 
          for (ix=0; ix<nx; ix++) 
            for (iz=0; iz<nz; iz++) 
              img[ih][ix][iz] = curr[ih+nh2*(ix+iz*nx2)];
        sf_floatwrite(img[0][0],nzx,image);
      }
    }
  } /* End of time loop */
  sf_warning(".");

  exit(0);
}
