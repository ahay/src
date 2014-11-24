// 2D 10th-order Lowrank FD coefficient

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
static std::valarray<float> vs;
static std::valarray<double> ks;
static std::valarray<double> kx, kz;
static float dt;

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
        res(a,b) = 2.0*cos(2.0*SF_PI*vs[rs[a]]*ks[cs[b]]*dt);
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

    int SIZE;
    par.get("size",SIZE,9); // stencil length 
    iRSF velf;
    oRSF outm;
    float dx, dz;


    int nz, nx;
    velf.get("n1",nz);
    velf.get("d1",dz);
    velf.get("n2",nx);
    velf.get("d2",dx);
    float dkz, dkx;
    dkx = 1.0/(dx*nx);
    dkz = 1.0/(dz*nz);

    int nxz = nx*nz;
    vs.resize(nxz);
    ks.resize(nxz);
    velf >> vs;
    
    int m = nxz;
    int n = nxz;

    int COUNT= 0;
    kx.resize(nxz);
    kz.resize(nxz);
    float kx1, kz1;
    float kx0 = -dkx*nx/2.0; 
    float kz0 = -dkz*nz/2.0; 
    float a = nx/4.0*dkx; 
    float b = nz/4.0*dkz; 
    int i=0;
    float dkxz=dkx+dkz;
    int SMK=0;
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
            kz1 = kz0+iz*dkz; 
            ks[i] = sqrtf(kx1*kx1+kz1*kz1);
            //ks[iz+nz*ix] = sqrtf(kx1*kx1+kz1*kz1);
            if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) COUNT++;
            if (ks[i] < (dkxz+0.00001)) SMK++;
            kx[i] = kx1;
            kz[i] = kz1;
            i++;
        }
    }
    vector<int> ksc(COUNT), smallk(SMK);
    sf_warning("COUNT=%d",COUNT);
    sf_warning("SMK=%d",SMK);
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
    sf_warning("nk=%d",nk);
    
    
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
    sf_warning("BAND=%d",band);
    
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
    DblNumMat kxtmpc(1,COUNT); for(int k=0; k<COUNT; k++) kxtmpc._data[k]=kx[ksc[k]];
    DblNumMat kztmpc(1,COUNT); for(int k=0; k<COUNT; k++) kztmpc._data[k]=kz[ksc[k]];
    int LEN = s1._m;
    DblNumMat Bx(LEN,nxz), B(LEN,nxz);
    DblNumMat Bxc(LEN,COUNT), Bc(LEN,COUNT);
    iC(ddgemm(dz,s1,kztmp,0.0,B));
    iC(ddgemm(dx,s2,kxtmp,0.0,Bx));
    iC(ddgemm(dz,s1,kztmpc,0.0,Bc));
    iC(ddgemm(dx,s2,kxtmpc,0.0,Bxc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=cos(2.0*SF_PI*(B._data[k]+Bx._data[k]));
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=cos(2.0*SF_PI*(Bc._data[k]+Bxc._data[k]));
    DblNumMat IBc(COUNT,LEN);    iC( ddpinv(Bc, 1e-16, IBc) );
    DblNumMat coef(ridx.size(),LEN);
    DblNumMat M2c;
    iC( sample(ridx,ksc,M2c) );
    
    iC(ddgemm(1.0,M2c,IBc,0.0,coef));

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
    sf_warning("i=%d",i);
    if (i>0) {
       for(int k=0; k<nxz; k++) { maxB._data[k]=2.0/fabs(maxB._data[k]); }
       for (int x=0; x<nxz; x++){
           for (int k=0; k<LEN; k++){
                  G(x,k) = G(x,k)*maxB._data[x];
           }
       } 
    }

    float dk2=sqrt(dkx*dkx+dkz*dkz);
    std::valarray<float> fMlr(nz*3),ktmp(nz);
    DblNumMat vratio(nz,3);
    for (int k=0; k < nz; k++) {
        ktmp[k] = k*dk2*2.0*SF_PI;
    } 
    float gt, seta, dangle=SF_PI/8;
    for (int aj=0; aj<3; aj++) {
        seta = dangle*aj;
        vratio(0,aj) = 1.0;
        for (int k=1; k < nz; k++) {
            gt = 0;
            for (int j=0; j<LEN; j++) {
                gt += G(0,j)*cos(s1(j,0)*ktmp[k]*cos(seta)*dz+s2(j,0)*ktmp[k]*sin(seta)*dx);
            }
            gt = acos(0.5*gt);
            vratio(k,aj) = gt/(vs[0]*ktmp[k]*dt);
        }
    }
    
    double *ldat = vratio.data();
    for (int k=0; k < nz*3; k++) {
        fMlr[k] = ldat[k];
    } 
    outm.put("n1",nz);
    outm.put("d1",dk2);
    outm.put("n2",3);
    outm.put("d2",dangle);
    outm << fMlr;
    return 0;
}





