// Lowrank decomposition for 3-D orthorhombic wave propagation with linearization. 
//   Copyright (C) 2010 University of Texas at Austin
//  
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//  
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//  
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include <time.h>

#include <rsf.hh>

#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

static std::valarray<float>  vx, dv, vz, yi1, dyi, dmu;
static std::valarray<float> kx, ky, kz;
static float dt;

int sample(vector<int>& rs, vector<int>& cs, FltNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0f);
    for(int a=0; a<nr; a++) {
        int i=rs[a];
        float wx = vx[i]*vx[i];
        float wy = wx;
        float wz = vz[i]*vz[i];
        float dw = dv[i]*vx[i];
        float e1 = yi1[i];
        float e2 = dyi[i];
        float e3 = dmu[i];

	for(int b=0; b<nc; b++) {
           int j=cs[b];
           float x=kx[j]*kx[j];
           float y=ky[j]*ky[j];
           float z=kz[j]*kz[j];
          
           float wxx = wx*x;
           float wyy = wy*y;
           float wzz = wz*z;
           float wxy = dw*y;
           float tmp1 = (wxx+wyy)*(1+2*e1);
	   float phi0 = tmp1*tmp1+wzz*wzz-2.0*(wxx+wyy)*wzz*(2*e1-1);
           phi0 = sqrt(phi0) + tmp1 + wzz; 
           float phi02 = phi0/2.0;
           phi0 = sqrt(phi02);
           float tmp2 = -2.0*tmp1*phi0*phi02+3.0*phi0*phi02*phi02+2*wzz*phi0*((wxx+wyy)*e1-phi02);
           float phiv = (abs(tmp2)>1.0e-7)?-wxy*(wxx*(1+2.0*e1)-phi02)*(-2.0*wzz*e1+(1+2.0*e1)*phi02)/tmp2:0;
           float phie = (abs(tmp2)>1.0e-7)?-wyy*((wxx*(1+2.0*e1)-phi02)*phi02+wzz*(-2.0*wxx*e1+phi02))/tmp2:0;
           float phig = (abs(tmp2)>1.0e-7)?wxx*wyy*(1+2.0*e1)*(-2*wzz*e1+(1+2.0*e1)*phi02)/tmp2:0;
           //float phig = (abs(tmp2)>1.0e-7)?wxx*wyy*(1+2.0*e1)*(-2*wzz*e1+(1+2.0*e1)*phi02)/(-tmp2):0;
           float phi = phi0 + phiv + phie*e2 + phig*e3;
	   res(a,b) = 2*(cos(phi*dt)-1); 
	}
    }
    return 0;
}

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0);
    int seed;

    par.get("seed",seed,time(NULL)); // seed for random number generator
    srand48(seed);

    float eps;
    par.get("eps",eps,1.e-4); // tolerance

    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step

    iRSF velz, velx("velx"), vely("vely"), etax("etax"), etay("etay"), muz("mu");

    int nz,nx,ny;
    velz.get("n1",nz);
    velz.get("n2",nx);
    velz.get("n3",ny);
    int m = nx*ny*nz;

    vx.resize(m);
    vz.resize(m);
    yi1.resize(m);
    dmu.resize(m);
    dv.resize(m);
    dyi.resize(m);

    velx >> vx;
    vely >> dv;
    velz >> vz;
    etax >> yi1;
    etay >> dyi;
    muz  >> dmu;
     

    iRSF fft("fft");

    int nkz,nkx,nky;
    fft.get("n1",nkz);
    fft.get("n2",nkx);
    fft.get("n3",nky);

    float dkz,dkx,dky;
    fft.get("d1",dkz);
    fft.get("d2",dkx);
    fft.get("d3",dky);
    
    float kz0,kx0,ky0;
    fft.get("o1",kz0);
    fft.get("o2",kx0);
    fft.get("o3",ky0);


    int n = nkx*nky*nkz;
    kx.resize(n);
    ky.resize(n);
    kz.resize(n);
    int i = 0; 
    for (int iy=0; iy < nky; iy++) {
        for (int ix=0; ix < nkx; ix++) {
            for (int iz=0; iz < nkz; iz++) {
                kx[i] = 2*SF_PI*(kx0+ix*dkx);
                ky[i] = 2*SF_PI*(ky0+iy*dky);
                kz[i] = 2*SF_PI*(kz0+iz*dkz);
                i++;
            }
        }
    }

    for (i=0; i < m; i++){
         dv[i]  -=  vx[i];
         dyi[i] -=  yi1[i];
         dmu[i] -=  1.0;
    }

        

    vector<int> lidx, ridx;
    FltNumMat mid;

    iC( lowrank(m,n,sample,(float) eps,npk,lidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();
    float *dmid = mid.data();

    std::valarray<float> fmid(m2*n2);
    for (int k=0; k < m2*n2; k++) {
	fmid[k] = dmid[k];
    }

    oRSF middle;
    middle.put("n1",m2);
    middle.put("n2",n2);
    middle.put("n3",1);
    middle << fmid;

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    FltNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );
    float *ldat = lmat.data();

    std::valarray<float> ldata(m*m2);
    for (int k=0; k < m*m2; k++) 
	ldata[k] = ldat[k];
    oRSF left("left");
    left.put("n1",m);
    left.put("n2",m2);
    left.put("n3",1);
    left << ldata;

    FltNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );
    float *rdat = rmat.data();

    std::valarray<float> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = rdat[k];
    oRSF right("right");
    right.put("n1",n2);
    right.put("n2",n);
    right.put("n3",1);
    right << rdata;

    return 0;
}
