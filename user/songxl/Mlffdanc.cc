// 2D high-order TTI Lowrank FFD coefficient

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
static float pi=SF_PI;
static float dt;
static double vx0, vz0, t0, q00;

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    //setvalue(res,0.0f);
    double c0 = cos(t0);
    double s0 = sin(t0);
    double wx0 = vx0*vx0;
    double wz0 = vz0*vz0;
    for(int a=0; a<nr; a++) {
        int i=rs[a];
        double wx = vx[i]*vx[i];
        double wz = vz[i]*vz[i];
        double qq = q[i];
        double tt = t[i];
        double c = cos(tt);
        double s = sin(tt);

	for (int b=0; b<nc; b++) {
            double x0 = kx[cs[b]];
            double z0 = kz[cs[b]];
            // rotation of coordinates
            double x = x0*c+z0*s;
            double z = z0*c-x0*s;
            double x00 = x0*c0+z0*s0;
            double z00 = z0*c0-x0*s0;

            z = wz*z*z;
            x = wx*x*x;
            double r = x+z;
            r = r+sqrt(r*r-qq*x*z);
            r = sqrt(0.5*r);
            z00 = wz0*z00*z00;
            x00 = wx0*x00*x00;
            double r0 = x00+z00;
            r0 = r0+sqrt(r0*r0-q00*x00*z00);
            r0 = sqrt(0.5*r0);
            res(a,b) = cos(2.0*pi*r*dt)/cos(2.0*pi*r0*dt);
           
	}
/*
	for (int b=0; b<nc; b++) {
            double x0 = kx[cs[b]];
            double z0 = kz[cs[b]];
            // rotation of coordinates
            double x = x0*c+z0*s;
            double z = z0*c-x0*s;
            double x00 = x0*c0+z0*s0;
            double z00 = z0*c0-x0*s0;

            res(a,b) = (x*x+z*z);

            z = wz*z*z;
            x = wx*x*x;
            double r = x+z;
            r = r+sqrt(r*r-qq*x*z);
            r = sqrt(0.5*r);
            z00 = wz0*z00*z00;
            x00 = wx0*x00*x00;
            double r0 = x00+z00;
            r0 = r0+sqrt(r0*r0-q00*x00*z00);
            r0 = sqrt(0.5*r0);
            if(r0<0.0000001) res(a,b)=0;
            else res(a,b) *= (cos(2.0*pi*r*dt)-1.0)/(cos(2.0*pi*r0*dt)-1.0);
           
	}
*/
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
    float pr;
    par.get("pr",pr,0.15); // time step

    int SIZE, DE;
    par.get("size",SIZE,9); // stencil length 
    par.get("de",DE,1); // stencil length 
    iRSF velz, velx("velx"), eta("eta"), seta("seta");
    oRSF outm, s1f("s1"), s2f("s2"), paraf("paras");
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

    double q0=0;
    vx0=0; vz0=0; t0=0; 
    for (int im=0; im < m; im++) {
        vx0 += vx[im]*vx[im];
        vz0 += vz[im]*vz[im];
         q0 += q[im];
         t0 += t[im];
    }
    vx0 = sqrt(vx0/((double) m));
    vz0 = sqrt(vz0/((double) m));
     q0 = q0/((double) m);
     t0 = t0/((double) m);

    /* from eta to q */
    for (int im=0; im < m; im++) {
        q[im] = 8*q[im]/(1.0+2*q[im]);
    }
    q00 = 8*q0/(1.0+2*q0);

    /* fram degrees to radians */
    if (DE)
    {
       for (int im=0; im < m; im++) {
           t[im] *= SF_PI/180.;
       }
       t0  *= SF_PI/180.;
    }
    std::valarray<float> parms(4);
    parms[0]=(float) vz0;
    parms[1]=(float) vx0;
    parms[2]=(float) q0;
    parms[3]=(float) t0;
    paraf.put("n1",4);
    paraf.put("n2",1);
    paraf << parms;

    
  

    int COUNT= 0;
    kx.resize(nxz);
    kz.resize(nxz);
    double kx1, kz1;
    float kx0 = -dkx*nx/2.0; 
    float kz0 = -dkz*nz/2.0; 
    float a = nx*pr*dkx; 
    float b = nz*pr*dkz; 
    int i=0;
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
            kz1 = kz0+iz*dkz; 
            if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) COUNT++;
            kx[i] = kx1;
            kz[i] = kz1;
            i++;
        }
    }
    vector<int> ksc(COUNT);
    sf_warning("COUNT=%d",COUNT);
    int nk=0; 
    i=0;
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
            kz1 = kz0+iz*dkz; 
            if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) {
               ksc[nk] = i;
               nk++;
            }
            i++;
               
        }
    }
    
    
    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m,n,sample,(double) eps,npk,cidx,ridx,mid) );

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
            //if((stmp[iz] > 0 || (stmp[iz] == 0 && stmp[ix] >=0)) && ((stmp[iz]*stmp[iz]+stmp[ix]*stmp[ix])< (float) (band*band)+0.00001)) gdc++;
            if(( (stmp[iz] > 0 && stmp[ix] == 0) || (stmp[iz] == 0 && stmp[ix] >=0))) gdc++;
        }
    }
    nk = 0;
    DblNumMat s1(gdc,1), s2(gdc,1); 
    for (int ix=0; ix<SIZE; ix++){
        for (int iz=0; iz<SIZE; iz++){
            //if((stmp[iz] > 0 || (stmp[iz] == 0 && stmp[ix] >=0)) && ((stmp[iz]*stmp[iz]+stmp[ix]*stmp[ix])<(float) (band*band)+0.00001)){
            if((  (stmp[iz] > 0 && stmp[ix] == 0) || (stmp[iz] == 0 && stmp[ix] >=0))){
              s1._data[nk]=stmp[iz];
              s2._data[nk]=stmp[ix];
              nk++;
            }
        }
    }
    

    DblNumMat kxtmp(1,nxz); for(int k=0; k<nxz; k++) kxtmp._data[k]=kx[k];
    DblNumMat kztmp(1,nxz); for(int k=0; k<nxz; k++) kztmp._data[k]=kz[k];
    int LENt = s1._m;
    int LEN = LENt+1;
    DblNumMat Bc(LENt,nxz), Bb(LENt,nxz),B(LEN,nxz);
    iC(ddgemm(dz,s1,kztmp,0.0,Bb));
    iC(ddgemm(dx,s2,kxtmp,0.0,Bc));
    for (int x=0; x<LENt; x++){
        for (int k=0; k<nxz; k++){
            B(x,k) = cos(2.0*pi*(Bc(x,k)+Bb(x,k)));
        }
    }
    for (int k=0; k<nxz; k++) B(LENt,k) = cos(2.0*pi*(kz[k]*dz+kx[k]*dx))+cos(2.0*pi*(kz[k]*dz-kx[k]*dx));
    //for(int k=0; k<B._m*B._n; k++) B._data[k]=cos(2.0*pi*(B._data[k]+Bc._data[k]));
    float ACCU=1e-4;
    DblNumMat WGT(nxz,1);
    for(int k=0; k<nxz; k++) WGT._data[k]=1.0;
    for(int k=0; k<COUNT; k++) WGT._data[ksc[k]]=1.0/ACCU;
/*
    DblNumMat M2c;
    iC( sample(ridx,ksc,M2c) );
    for (int x=0; x<M2._m; x++){
        for (int k=0; k<nxz; k++){
            M2(x,k) = 0.0;
        }
    } 
*/
    DblNumMat M2c(M2._m,nxz);
    for (int x=0; x<M2._m; x++){
        for (int k=0; k<nxz; k++){
            M2c(x,k) = M2(x,k)*WGT._data[k];
        }
    } 
    Bc.resize(LEN,nxz);  
    for (int x=0; x<B._m; x++){
        for (int k=0; k<nxz; k++){
            Bc(x,k) = B(x,k)*WGT._data[k];
        }
    } 

    DblNumMat IBc(nxz,LEN);    iC( ddpinv(Bc, 1e-16, IBc) );
    DblNumMat coef(ridx.size(),LEN);
//    DblNumMat M2c;
//    iC( sample(ridx,ksc,M2c) );
    
    iC(ddgemm(1.0,M2c,IBc,0.0,coef));

    DblNumMat G(nxz,LEN), tmpG(mid._m,LEN);
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,G));


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
    s1f.put("n1",LENt);
    s1f.put("n2",1);
    s2f.put("n1",LENt);
    s2f.put("n2",1);
    //std::valarray<int> fs(LEN);
    std::valarray<float> fs(LENt);
    ldat = s1.data();
    for (int k=0; k < LENt; k++) {
        //fs[k] = (int) ldat[k];
        fs[k] = ldat[k];
    } 
    s1f << fs;
    ldat = s2.data();
    for (int k=0; k < LENt; k++) {
        fs[k] =  ldat[k];
    } 
    s2f << fs;

  //  */
    return 0;
}





