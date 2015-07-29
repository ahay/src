// 2D high-order TTI Lowrank FD coefficient

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

//FltNumVec vs; //c
//FltNumVec ks; //k
static std::valarray<float> vz,vx,q,t;
static std::valarray<double> kx, kz;
static float dt;

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
        int i=rs[a];
        float wx = vx[i]*vx[i];
        float wz = vz[i]*vz[i];
        float qq = q[i];
        float tt = t[i];
        double c = cos(tt);
        double s = sin(tt);

	for (int b=0; b<nc; b++) {
            double x0 = kx[cs[b]];
            double z0 = kz[cs[b]];
            // rotation of coordinates
            double x = x0*c+z0*s;
            double z = z0*c-x0*s;

            z = wz*z*z;
            x = wx*x*x;
            double r = x+z;
            r = r+sqrt(r*r-qq*x*z);
            r = sqrt(0.5*r);
            res(a,b) = 2.0*cos(2.0*SF_PI*r*dt);
           
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
    par.get("eps",eps,1.e-6); // tolerance

    int npk;
    par.get("npk",npk,50); // maximum rank

    par.get("dt",dt); // time step

    int SIZE, DE;
    par.get("size",SIZE,17); // stencil length 
    par.get("de",DE,1); // stencil length 
    iRSF velz, velx("velx"), eta("eta"), seta("seta");
    oRSF outm, s1f("s1"), s2f("s2");
    float dx, dz;


    int nz, nx;
    velz.get("n1",nz);
    velz.get("d1",dz);
    velz.get("n2",nx);
    velz.get("d2",dx);
    float dkz, dkx;
    dkx = 1.0/(dx*nx);
    dkz = 1.0/(dz*nz);

    int nxz = nx*nz;
    vx.resize(nxz);
    vz.resize(nxz);
    t.resize(nxz);
    q.resize(nxz);
    velz >> vz;
    velx >> vx;
    eta  >> q;
    seta >> t;


    
    int m = nxz;
    int n = nxz;

    /* from eta to q */
    for (int im=0; im < m; im++) {
        q[im] = 8*q[im]/(1.0+2*q[im]);
    }

    /* fram degrees to radians */
    if (DE)
    {
       for (int im=0; im < m; im++) {
           t[im] *= SF_PI/180.;
       }
    }
  

    int COUNT= 0;
    kx.resize(nxz);
    kz.resize(nxz);
    float kx1, kz1;
    float kx0 = -dkx*nx/2.0; 
    float kz0 = -dkz*nz/2.0; 
    float a = nx/3.0*dkx; 
    float b = nz/3.0*dkz; 
    int i=0;
    float dkxz=dkx+dkz;
    int SMK=0;
    std::valarray<double> ks(nxz);
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
            kz1 = kz0+iz*dkz; 
            ks[iz+nz*ix] = sqrtf(kx1*kx1+kz1*kz1);
            if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) COUNT++;
            if (ks[i] < (dkxz+0.00001)) SMK++;
            kx[i] = kx1;
            kz[i] = kz1;
            i++;
        }
    }
    vector<int> ksc(COUNT), smallk(SMK);
    int nk=0, mk=0; 
    i=0;
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
            kz1 = kz0+iz*dkz; 
            if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) {
               ksc[nk] = i;
               nk++;
            }
            if (ks[i] < (dkxz+0.00001)){ 
               smallk[mk] = i;
               mk++;
            }
            i++;
               
        }
    }
    
    
    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m,n,sample,eps,npk,cidx,ridx,mid) );

    DblNumMat M1(m,cidx.size());
    vector<int> rs(m);
    for(int k=0; k<m; k++) rs[k]=k;
    iC( sample(rs,cidx,M1) );
    DblNumMat M2(ridx.size(),n);
    vector<int> cs(n);
    for(int k=0; k<n; k++) cs[k]=k;
    iC( sample(ridx,cs,M2) );


/*  Next */

    //float stmp[] = {0,1,2,3,4,5};
    std::valarray<float> stmp(SIZE);
    int band = (SIZE-1)/2;
    
    for (int ix=0; ix<SIZE; ix++) {stmp[ix]= (float) (ix-band);}
    int gdc=0;
    for (int ix=0; ix<SIZE; ix++){
        for (int iz=0; iz<SIZE; iz++){
            if((stmp[iz] > 0 || (stmp[iz] == 0 && stmp[ix] >=0)) && ((stmp[iz]*stmp[iz]+stmp[ix]*stmp[ix])< (float) (band*band)+0.00001)) gdc++;
        }
    }
    nk = 0;
    DblNumMat s1(gdc,1), s2(gdc,1); 
    for (int ix=0; ix<SIZE; ix++){
        for (int iz=0; iz<SIZE; iz++){
            if((stmp[iz] > 0 || (stmp[iz] == 0 && stmp[ix] >=0)) && ((stmp[iz]*stmp[iz]+stmp[ix]*stmp[ix])<(float) (band*band)+0.00001)){
              s1._data[nk]=stmp[iz];
              s2._data[nk]=stmp[ix];
              nk++;
            }
        }
    }
    

    DblNumMat kxtmp(1,nxz); for(int k=0; k<nxz; k++) kxtmp._data[k]=kx[k];
    DblNumMat kztmp(1,nxz); for(int k=0; k<nxz; k++) kztmp._data[k]=kz[k];
 //   DblNumMat kxtmpc(1,COUNT); for(int k=0; k<COUNT; k++) kxtmpc._data[k]=kx[ksc[k]];
 //   DblNumMat kztmpc(1,COUNT); for(int k=0; k<COUNT; k++) kztmpc._data[k]=kz[ksc[k]];
    int LEN = s1._m;
    DblNumMat Bc(LEN,nxz), B(LEN,nxz);
//    DblNumMat Bxc(LEN,COUNT), Bc(LEN,COUNT);
    iC(ddgemm(dz,s1,kztmp,0.0,B));
    iC(ddgemm(dx,s2,kxtmp,0.0,Bc));
//    iC(dgemm(dz,s1,kztmpc,0.0,Bc));
//    iC(dgemm(dx,s2,kxtmpc,0.0,Bxc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=cos(2.0*SF_PI*(B._data[k]+Bc._data[k]));
//    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=cos(2.0*pi*(Bc._data[k]+Bxc._data[k]));
    float ACCU=1e-4;
    DblNumMat WGT(nxz,1);
    for(int k=0; k<nxz; k++) WGT._data[k]=1.0;
    for(int k=0; k<COUNT; k++) WGT._data[ksc[k]]=1.0/ACCU;
    DblNumMat M2c;
    iC( sample(ridx,ksc,M2c) );
    for (int x=0; x<M2._m; x++){
        for (int k=0; k<nxz; k++){
            M2(x,k) = 0.0;
        }
    } 
    for (int x=0; x<M2._m; x++){
        for (int k=0; k<COUNT; k++){
            M2(x,ksc[k]) = M2c(x,k)*WGT._data[ksc[k]];
        }
    } 
    for (int x=0; x<B._m; x++){
        for (int k=0; k<nxz; k++){
            Bc(x,k) = B(x,k)*WGT._data[k];
        }
    } 

    DblNumMat IBc(nxz,LEN);    iC( ddpinv(Bc, 1e-16, IBc) );
    DblNumMat coef(ridx.size(),LEN);
//    DblNumMat M2c;
//    iC( sample(ridx,ksc,M2c) );
    
    iC(ddgemm(1.0,M2,IBc,0.0,coef));

    DblNumMat G(nxz,LEN), tmpG(mid._m,LEN);
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,G));

    Bc.resize(LEN,SMK);

    for(int k=0; k<LEN; k++) {
       for (int j=0; j<SMK; j++) {
           Bc(k,j) =B(k,smallk[j]);
       }
    }
    DblNumMat tmpB(nxz,SMK), maxB(nxz,1);
    iC(ddgemm(1.0,G,Bc,0.0,tmpB));
    float tmpmax;
    for (int k=0; k<nxz; k++) {
        tmpmax=-9999999.0;
        for (int j=0; j<SMK; j++) {
            if (fabs(tmpB(k,j)) > tmpmax) tmpmax=fabs(tmpB(k,j));
        }
        maxB._data[k]=tmpmax;
    }
    i=0;
    for(int k=0; k<nxz; k++) if (maxB._data[k]>2.0) i++;
    if (i>0) {
       for(int k=0; k<nxz; k++) { maxB._data[k]=2.0/fabs(maxB._data[k]); }
       for (int x=0; x<nxz; x++){
           for (int k=0; k<LEN; k++){
                  G(x,k) = G(x,k)*maxB._data[x];
           }
       } 
    }

    std::valarray<float> fMlr(nxz*LEN);
    double *ldat = G.data();
    for (int k=0; k < nxz*LEN; k++) {
        fMlr[k] = ldat[k];
    } 
    outm.put("n1",nxz);
    outm.put("n2",LEN);
    outm << fMlr;
    //s1f.type(SF_INT);
    //s2f.type(SF_INT);
    s1f.put("n1",LEN);
    s1f.put("n2",1);
    s2f.put("n1",LEN);
    s2f.put("n2",1);
    //std::valarray<int> fs(LEN);
    std::valarray<float> fs(LEN);
    ldat = s1.data();
    for (int k=0; k < LEN; k++) {
        //fs[k] = (int) ldat[k];
        fs[k] = ldat[k];
    } 
    s1f << fs;
    ldat = s2.data();
    for (int k=0; k < LEN; k++) {
        fs[k] =  ldat[k];
    } 
    s2f << fs;

  //  */
    return 0;
}





