// Lowrank decomposition for 2D PERM in DTI media
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
   * :WARNING: 
   *  1.nx*nz*nh and nk should be smaller than the maximum value of 4-byte int.
   *  2.Large models need parallelized version.
   * ========================================================================== 
   * I/O File                             Description
   * --------------------------------------------------------------------------
   *  in        :    [nz][nx]             vp 
   *  vnmo      :    [nz][nx]             vnmo 
   *  eta       :    [nz][nx]             eta
   *  theta     :    [nz][nx]             theta
   *  left      :    [nzx][n2]            lowrank factorized left matrix
   *  right     :    [n2][nk]             lowrank factorized right matrix
   * ========================================================================== 
   */

#include <time.h>
#include <rsf.hh>
#include "vecmatop.hh"
#include "serialize.hh"
using namespace std;

static valarray<float> v,vnmo,eta,theta;
static int nh,nx,nz;
static int nkh,nkx,nkz;
static float dt,kh,kx,kz;
static float dkh,dkx,dkz;
static float kh0,kx0,kz0;
static float dx,dz,dh,ox,oz,oh;
static double minv,maxv,maxvs2pvr2;
static float freqmax;
static int equation;


static int 
sample(vector<int>& rs, vector<int>& cs, DblNumMat &res)
{ // phase operator evaluation in DTI media 
  int ix, ih, iz;
  double phi;
  int nr = rs.size();
  int nc = cs.size();
  res.resize(nr,nc);
  setvalue(res,0.0);

  double hnow;
  sf_warning("nr=%d nc=%d",nr,nc);

  double alpha,sin2_alpha,cos2_alpha;
  double coef1,coef2,coef3;
  double vt_s,vt_r,eta_s,eta_r,vnmo_s,vnmo_r;

  for (int b=0; b<nc; b++) { // column loop

    int i = cs[b];
    ih = i%nkh; i = static_cast<int>(floor(static_cast<double>(i)/nkh));
    ix = i%nkx; i = static_cast<int>(floor(static_cast<double>(i)/nkx)); 
    iz = i%nkz; 
    kh = kh0+ih*dkh;
    kx = kx0+ix*dkx;
    kz = kz0+iz*dkz;

    double kp = kx+kh;
    double km = kx-kh;
    double kxh = kx*kh;

    alpha = atan2(kh,kz); 
    if (alpha > 0.5*SF_PI)
      alpha -= SF_PI;
    else if (alpha < -0.5*SF_PI)
      alpha += SF_PI;
    sin2_alpha = sin(alpha);   sin2_alpha *= sin2_alpha;
    cos2_alpha = cos(alpha);    cos2_alpha *= cos2_alpha;

    kp *= kp;
    km *= km;
    kz *= kz; 
    kh *= kh; 
    kx *= kx;
    double x,y;

    for (int a=0; a<nr; a++) { // row loop

      int i = rs[a];
      ih = i%nh; i = static_cast<int>(floor(static_cast<double>(i)/nh));
      ix = i%nx; i = static_cast<int>(floor(static_cast<double>(i)/nx));
      iz = i%nz; 
      hnow=(oh+ih*dh)/dx;
      int is = SF_MIN(SF_MAX(0,static_cast<int>(ix-hnow)),nx-1);
      int ir = SF_MIN(SF_MAX(0,static_cast<int>(ix+hnow)),nx-1);
      double vs = v[is+iz*nx];
      double vr = v[ir+iz*nx];
      vnmo_s = vnmo[is+iz*nx];    vnmo_s *= vnmo_s;
      vt_s = v[is+iz*nx];         vt_s *= vt_s;
      eta_s = eta[is+iz*nx];
      coef1 = vnmo_s*(1.+2.*eta_s);
      coef2 = 2.*vnmo_s*vt_s*(1. - 2*eta_s);
      coef3 = 4.*vt_s*vt_s;
      vs = 0.5*(coef1*sin2_alpha + vt_s*cos2_alpha)
        + 0.25*sqrt(4.*coef1*coef1*sin2_alpha*sin2_alpha
            + 4.*coef2*sin2_alpha*cos2_alpha + coef3*cos2_alpha*cos2_alpha );
      vnmo_r = vnmo[ir+iz*nx];    vnmo_r *= vnmo_r;
      vt_r = v[ir+iz*nx];         vt_r *= vt_r;
      eta_r = eta[ir+iz*nx];
      coef1 = vnmo_r*(1.+2.*eta_r);
      coef2 = 2.*vnmo_r*vt_r*(1. - 2*eta_r);
      coef3 = 4.*vt_r*vt_r;
      vr = 0.5*(coef1*sin2_alpha + vt_r*cos2_alpha)
        + 0.25*sqrt(4.*coef1*coef1*sin2_alpha*sin2_alpha
            + 4.*coef2*sin2_alpha*cos2_alpha + coef3*cos2_alpha*cos2_alpha );
      double vm = vs-vr;

      x=(kh + kz)*(kx + kz)*vr*vs;
      y=kxh*vm +2*kz*(vs+vr);
      if(y<1e-13)
        y=1e-13;
      phi=x/y;

      if(x<0 || y<0)
        cerr<<"x="<<x<<" y="<<y<<" phi="<<phi<<endl;

      double kzmin = sqrt(kh*kx);

      if (kz>=kzmin &&
          ((kx+kh+kz)<=(2*freqmax*freqmax*maxvs2pvr2)) &&
          (kp<=(freqmax*freqmax*4/minv/minv)) &&
          (km<=(freqmax*freqmax*4/minv/minv)))
      {
          res(a,b) = 2*cos(2*SF_PI*sqrt(phi)*dt); 
      } else {
          res(a,b) = 0.0;
      }
    } // End of column(nc) loop
  } // End of row(nr) loop
  return 0;
}



int
main(int argc, char** argv)
{
  sf_init(argc,argv);

  int seed;
  int npk;
  int nxz;
  float tol;

  iRSF par(0);
  iRSF vp_file;  // vz 
  iRSF vnmo_file("vnmo");  // vnmo 
  iRSF eta_file("eta");  // eta 
  iRSF theta_file("theta");  // theta  
  oRSF left("left");  // left maxtrix 
  oRSF right;  // right matrix 


  par.get("seed",seed,time(NULL)); // seed for random number generator
  srand48(seed);
  par.get("tol",tol,1.e-4); // tolerance
  par.get("npk",npk,20); // maximum rank
  par.get("freqmax",freqmax,60.); // if subtract one
  par.get("dt",dt); // time step
  par.get("nh",nh); // half-offset 
  par.get("dh",dh);
  par.get("oh",oh); 

  vp_file.get("n1",nx);
  vp_file.get("n2",nz);
  vp_file.get("d1",dx);
  vp_file.get("d2",dz);
  vp_file.get("o1",ox);
  vp_file.get("o2",oz);
  nxz = nx*nz;

  v.resize(nxz);
  vnmo.resize(nxz);
  eta.resize(nxz);
  theta.resize(nxz);

  vp_file >> v;    
  vnmo_file >> vnmo;
  eta_file >> eta;
  theta_file >> theta;
  for (int k=0;k<nxz;k++)
    theta[k] *= SF_PI/180.0;

  nkh = kiss_fft_next_fast_size((nh+1)/2)+1;                                             
  nkx = kiss_fft_next_fast_size(nx);                                             
  nkz = kiss_fft_next_fast_size(nz);                                             
  kh0 = 0.;          kx0 = -0.5/dx;          kz0 = -0.5/dz;                 
  dkh = 1./(2*(nkh-1)*dh);      dkx = 1./(nkx*dx);      dkz = 1./(nkz*dz); 

  int m = nxz*nh;
  int n = nkh*nkx*nkz;
  sf_warning("nx=%d nz=%d nh=%d m=%d n=%d",nx,nxz,nh,m,n);

  // set up coefficients for stability 
  int ih,ix,iz,is,ir;
  float vp_s,vp_r;
  float hh;
  float tmpcoef;
  float vx_pt,eta_pt;

  minv = v.min();
  maxv = v.max();
  maxvs2pvr2 = 0.;

  for (ih=0;ih<nh;ih++) {
    hh = (oh+ih*dh)/dx;
    for (ix=0;ix<nx;ix++) {
      is = static_cast<int>(SF_MIN(SF_MAX(0,ix-hh),nx-1));
      ir = static_cast<int>(SF_MIN(SF_MAX(0,ix+hh),nx-1));
      for (iz=0;iz<nz;iz++) {
        vp_s = v[is+iz*nx];
        vp_r = v[ir+iz*nx];
        eta_pt = eta[ix+iz*nx];
        vx_pt = sqrt(1.+2.*eta_pt)*vnmo[ix+iz*nx];
        if (maxv<vx_pt)
          maxv = vx_pt;
        vp_s *= vp_s;
        vp_r *= vp_r;
        tmpcoef = 1./vp_s + 1./vp_r;
        if (maxvs2pvr2<tmpcoef)
          maxvs2pvr2 = tmpcoef;
      }
    }
  }

  // lowrank decomposition
  vector<int> lidx, ridx;
  DblNumMat mid;

  iC( ddlowrank(m,n,sample,tol,npk,lidx,ridx,mid) );

  int m2=mid.m();
  int n2=mid.n();

  vector<int> midx(m), nidx(n);
  for (int k=0; k < m; k++) 
    midx[k] = k;
  for (int k=0; k < n; k++) 
    nidx[k] = k;    

  // collect left matrix left=left x mid 
  DblNumMat lmat(m,m2);
  iC ( sample(midx,lidx,lmat) );

  DblNumMat lmat2(m,n2);
  iC( ddgemm(1.0, lmat, mid, 0.0, lmat2) );

  double* ldat = lmat2.data();
  valarray<float> ldata(m*n2);
  for (int k=0; k < m*n2; k++) 
    ldata[k] = ldat[k];

  left.put("n1",m);
  left.put("n2",n2);
  left << ldata;

  // collect right matrix
  DblNumMat rmat(n2,n);
  iC ( sample(ridx,nidx,rmat) );
  double* rdat = rmat.data();

  valarray<float> rdata(n2*n);    
  for (int k=0; k < n2*n; k++) 
    rdata[k] = rdat[k];

  right.put("n1",n2);
  right.put("n2",n);
  right << rdata;

  return 0;

}
