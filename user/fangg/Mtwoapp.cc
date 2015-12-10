// 2D weighted least square Lowrank FD coefficient

//   Copyright (C) 2015 China Geological Survey
//   Author: Gang Fang
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

static std::valarray<float> vs;
static std::valarray<double> ks, kss;
static std::valarray<double> kx, kz;
static float dt;
static float dfrq, a0;
static float twopi = 2.0*SF_PI;

static double sinc(double x)
{
    if (fabs(x)<=SF_EPS) return 1.0;
    return sinf(x)/(x+SF_EPS);
}

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{ 
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr, nc);
    setvalue(res, 0.0);
    for (int a=0; a<nr; a++) {
	for (int b=0; b<nc; b++) {
	    res(a,b) = 2.0*cos(2.0*SF_PI*vs[rs[a]]*ks[cs[b]]*dt);
	}
    }
    return 0;
}

float wfun(float kx, float kz, float k0, float a0)
{
    float s1, s2, s3, s4;
    s1 = -1.0*((kx-k0)*(kx-k0)+(kz-k0)*(kz-k0))/a0;
    s2 = -1.0*((kx+k0)*(kx+k0)+(kz+k0)*(kz+k0))/a0;
    s3 = -1.0*((kx-k0)*(kx+k0)+(kz-k0)*(kz+k0))/a0;
    s4 = -1.0*((kx+k0)*(kx-k0)+(kz+k0)*(kz-k0))/a0;
    return (exp(s1)+exp(s2)+exp(s3)+exp(s4))/4.0;

}

int main(int argc, char** argv)
{
    sf_init(argc, argv); // Initialize RSF
    
    iRSF par(0);
    int  seed;

    par.get("seed",seed,time(NULL)); // seed for random number generator
    srand48(seed);

    float eps;
    par.get("eps", eps, 1.e-6); // tolerance
    int npk;
    par.get("npl", npk, 50); //maxium rank

    par.get("dt",dt); // time step

    int size;
    par.get("size", size, 6); //stencil length
    iRSF velf;
    oRSF outm;
    oRSF fsx("sx");
    oRSF fsz("sz");
    
    oRSF Mexact("Mexact"); // Exact operator: W(x,k)
    oRSF MGB("MGB");       // lowrank FD operator: GB
    //oRSF MWAW("Mwaw");     // lowrank operator: WAW
    //oRSF MGsB("MGsB");     // LFD with Sigm filter: GsB
    
    float dx, dz;

    float wavnumcut;
    par.get("wavnumcut",wavnumcut,1.0); // wavenumber cut percentile

    bool weight;
    par.get("weight", weight, true);
    /* using weighted least square */
    par.get("dfrq", dfrq, 15);
    /* dominant frequency */
    par.get("a0", a0, 1000); 
    /* weight parameters */
    int filter;
    par.get("filter", filter, 0);
    /*0- close; 1- sigma filter; 2- cos filter; 3- expontial filter*/
    
    float p;
    par.get("p", p, 0.1);
    /* power of sigma factor, 0<p<=1 */

    int nz, nx;
    velf.get("n1", nz);
    velf.get("d1", dz);
    velf.get("n2", nx);
    velf.get("d2", dx);
    float dkz, dkx;
    dkx = 1.0/(dx*nx);
    dkz = 1.0/(dz*nz);

    int nxz = nx*nz;
    vs.resize(nxz);
    ks.resize(nxz);
    velf >> vs;

    int m = nxz;
    int n = nxz;

    sf_warning("=======================");
    sf_warning(" 2D Lowrank FD coefficients with weighted least square and sigma filter");
    sf_warning("weight=%d",weight);
    sf_warning("wavnumcut=%f", wavnumcut);
    sf_warning("n=%d",n);
    sf_warning("m=%d",m);
    sf_warning("dfrq=%f", dfrq);
    switch (filter) {
	case 1 : // sigma filter
	    sf_warning("filter=%d, sinc filter, p=%f", filter, p);
	    break;
	case 2 : // cose filter
	    sf_warning("filter=%d, cose filter, p=%f", filter, p);
	    break;
	case 3 : // exponential filter
	    sf_warning("filter=%d, exp filter, p=%f",  filter, p);
	default :
	    sf_warning("filter=%d, sinc filter, p=%f", filter, p);
    } 
    sf_warning("sigm_p=%f", p);
    sf_warning("=======================");

    int COUNT = 0;
    kx.resize(nxz);
    kz.resize(nxz);
    float kx1, kz1;
    float kx0 = -dkx*nx/2.0;
    float kz0 = -dkz*nz/2.0;
    float a = nx/3.0*dkx*wavnumcut;
    float b = nx/3.0*dkz*wavnumcut;
    
    int i = 0;
    float dkxz = dkx+dkz;
    int SMK = 0;
    for (int ix=0; ix<nx; ix++) {
	kx1 = kx0+ix*dkx;
	for (int iz=0; iz<nz; iz++) {
	    kz1 = kz0+iz*dkz;
	    ks[i] = sqrtf(kx1*kx1+kz1*kz1);
	    if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) COUNT++;
	    if (ks[i] < (dkxz+0.00001)) SMK++;
	    kx[i] = kx1;
	    kz[i] = kz1;
	    i++;
	}
    }

    vector<int> ksc(COUNT), smallk(SMK);
    sf_warning("COUNT=%d", COUNT);
    sf_warning("SMK=%d", SMK);

    int nk=0, mk=0;
    i=0;
    for (int ix=0; ix<nx; ix++) {
	kx1 = kx0+ix*dkx;
	for (int iz=0; iz<nz; iz++) {
	    kz1 = kz0+iz*dkz;
	    if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) {
		ksc[nk] = i;
		nk++;
	    }
	    if (ks[i] < (dkxz+0.00001)) {
		smallk[mk] = i;
		mk++;
	    }
	    i++;
	}
    }
    sf_warning("nk=%d", nk);
    
    /* --- Exact Mat --- */
    vector<int> rs(m);
    for (int k=0; k<m; k++) rs[k]=k;
    vector<int> cs(n);
    for(int k=0; k<n; k++) cs[k]=k;

    int smlm = 250; 
    DblNumMat W(smlm, n);
    vector<int> smallx(smlm);
    i = 0;
    for (int k=0; k<smlm; k++) { smallx[k]=rs[i]; i+=(int)m/smlm; }
    iC( sample(smallx, cs, W) );
    valarray<float> fW(smlm*n);
    double *ddat=W.data();
    for (int k=0; k<n*smlm; k++) {
	fW[k] = ddat[k];
    }
    float maxexact, minexact;
    maxexact = fW.max();
    minexact = fW.min();
    sf_warning("Exact matrix: max value=%f; min value=%f",maxexact, minexact);

    const char *label1 = "Distance";
    const char *label2 = "Wavenumber";
    const char *unit1 = "km" ; 
    const char *unit2 = "1/km" ;
    
    Mexact.put("n1", smlm);
    Mexact.put("n2", n);
    Mexact.put("d1", dx);
    Mexact.put("d2", dkx*twopi);
    Mexact.put("o1", 0);
    float o2 = -1*dkx*nx*SF_PI ;  
    Mexact.put("o2", o2);
    Mexact.put("label1", label1);
    Mexact.put("unit1", unit1);
    Mexact.put("label2", label2);
    Mexact.put("unit2", unit2);
    Mexact << fW;


    /* --- Lowrank decomposition --- */
   
    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m, n, sample, eps, npk, cidx, ridx, mid) );

    DblNumMat M1(m, cidx.size());
    //vector<int> rs(m);
    //for (int k=0; k<m; k++) rs[k]=k;
    iC( sample(rs, cidx, M1) );
    DblNumMat M2(ridx.size(),n);
    //vector<int> cs(n);
    //for(int k=0; k<n; k++) cs[k]=k;
    iC( sample(ridx,cs,M2) );

    /* --- FD coefficients --- */

    /* Define Grids */
    std::valarray<float> stmp(size);
    for (int ix=0; ix<size; ix++) {stmp[ix] = (float)ix;}
    int gdc = 0;
    for (int ix=0; ix<size; ix++) {
	for (int iz=0; iz<size; iz++) {
	    if( (stmp[iz]>=0 && stmp[ix]==0) || (stmp[iz]==0 && stmp[ix]>=0) ) gdc++;
	}
    }
    nk = 0;
    gdc += 2;
   
    DblNumMat sx(gdc, 1), sz(gdc, 1);
    for (int ix=0; ix<size; ix++) {
	for (int iz=0; iz<size; iz++) {
	    if( (stmp[iz]>=0 && stmp[ix]==0) || (stmp[iz]==0 && stmp[ix]>=0) ) {
		sx._data[nk] = stmp[ix];
		sz._data[nk] = stmp[iz];
		nk++;
	    }
	}
    }
    sz._data[nk] = 1; sz._data[nk+1] = -1;
    sx._data[nk] = 1; sx._data[nk+1] = 1;

    cerr<<"[x,z]="; for(int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<","<<sz._data[k]<<"] ";    
    cerr<<endl;
    cerr<<"gdc "<<gdc<<" "<<endl;
    cerr<<"nk "<<nk<<" "<<endl;

    /* G */
    DblNumMat kxtmp(1, nxz); for (int k=0; k<nxz; k++) kxtmp._data[k]=kx[k];
    DblNumMat kztmp(1, nxz); for (int k=0; k<nxz; k++) kztmp._data[k]=kz[k];
    DblNumMat kxtmpc(1, COUNT); for (int k=0; k<COUNT; k++) kxtmpc._data[k] = kx[ksc[k]];
    DblNumMat kztmpc(1, COUNT); for (int k=0; k<COUNT; k++) kztmpc._data[k] = kz[ksc[k]];
    
    int LEN = sz._m;
    DblNumMat Bx(LEN, nxz), B(LEN, nxz);
    DblNumMat Bxc(LEN, COUNT), Bc(LEN, COUNT);
    iC(ddgemm(dz, sz, kztmp,  0.0,  B));
    iC(ddgemm(dx, sx, kxtmp,  0.0,  Bx));
    iC(ddgemm(dz, sz, kztmpc, 0.0,  Bc));
    iC(ddgemm(dx, sx, kxtmpc, 0.0,  Bxc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=cos(2.0*SF_PI*(B._data[k]+Bx._data[k]));
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=cos(2.0*SF_PI*(Bc._data[k]+Bxc._data[k]));
    DblNumMat IBc(COUNT, LEN); iC( ddpinv(Bc, 1e-16, IBc) );
    DblNumMat coef(ridx.size(),LEN);
    DblNumMat M2c;
    iC( sample(ridx,ksc,M2c) );
    //

    DblNumMat G(nxz,LEN), tmpG(mid._m,LEN);
    //iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    //iC(ddgemm(1.0,M1,tmpG,0.0,G));
    
    /* weight */
    int itmp;
    float w=0.0, k0=0.0, kkx=0.0, kkz=0.0;
    valarray<float> wf(Bc._n*ridx.size());
    DblNumMat wBc(LEN, COUNT), wM2c(1, COUNT);
    DblNumMat tmpcoef(1, LEN);

    setvalue(coef, 0.0);
    setvalue(tmpG, 0.0);
    if ( weight==true ) {
	for (int ixm=0; ixm<(int)ridx.size(); ixm++) {
	    k0 = 2.0*SF_PI*dfrq/vs[ridx[ixm]];
	    for (int ik=0; ik<COUNT; ik++) {
		kkx = 2*SF_PI*kxtmpc._data[ik];
		kkz = 2*SF_PI*kztmpc._data[ik];
		w = wfun(kkx, kkz, k0, a0);
		wM2c(0,ik) = w*M2c(ixm, ik);
		wf[ik+ixm*Bc._n] = w;
		for (int il=0; il<Bc._m; il++) {
		    itmp = ik*B._m + il;
		    wBc._data[itmp] = w*Bc._data[itmp];
		}
	    }
	    setvalue(IBc, 0.0);
	    iC( ddpinv(wBc, 1e-16, IBc) );
	    setvalue(tmpcoef, 0.0);
	    iC( ddgemm(1.0, wM2c, IBc, 0.0, tmpcoef) );
	    for (int il=0; il<LEN; il++) {
		coef(ixm,il) = tmpcoef(0,il);
	    }
	}
    } else {
	iC(ddgemm(1.0,M2c,IBc,0.0,coef));
    }
    
    iC(ddgemm(1.0, mid, coef, 0.0, tmpG));
    iC(ddgemm(1.0,M1,tmpG, 0.0,G));

    DblNumMat sigm(LEN, 1);
    sigm._data[0]=1.0 ;
    int sizep = size -1;

    switch (filter) {
	case 1 : // sigma filter
	    for (int k=1; k<LEN; k++) {
		sigm._data[k] = fabs(sinc(sx._data[k]*SF_PI/sizep)*sinc(sz._data[k]*SF_PI/sizep));
		//sf_warning("k=%i input=%g sincx=%g sincz=%g sinc=%g sigmk=%g",k, sx._data[k]*SF_PI/sizep, sinc(sx._data[k]*SF_PI/sizep), sinc(sz._data[k]*SF_PI/sizep), sinc(sx._data[k]*SF_PI/sizep)*sinc(sz._data[k]*SF_PI/sizep), sigm._data[k]);
	    }
	    break;
	case 2 : // cose filter
	    for (int k=1; k<LEN; k++) {
		sigm._data[k] = 0.5*(1+cos(sx._data[k]*SF_PI/sizep))*0.5*(1+cos(sz._data[k]*SF_PI/sizep));
	    }
	    break;
	case 3 : // exponential filter
	    double alpha;
	    alpha = -1.0*log(1.0e-15);
	    sf_warning("ALPHA=%g", alpha);
	    for (int k=1; k<LEN; k++) {
		sigm._data[k] = exp(-1.0*alpha*pow(sx._data[k]/sizep,2.0))*exp(-1.0*alpha*pow(sz._data[k]/sizep,2.0));
	    }
	default :
	   for (int k=0; k<LEN; k++) {
		sigm._data[k] = 1.0; 
	    }
    } 

    for(int k=0; k<LEN; k++) {
	sigm._data[k] = pow(sigm._data[k], p);
    }

    if(filter!=0) {
	for(int k=0; k<LEN; k++) {
	    for(int ix=0; ix<nxz; ix++) {
		G._data[k*nxz+ix] = G._data[k*nxz+ix]*sigm._data[k];
	    }
	}
    }
    
    cerr<<"sigm="; for(int k=0; k<LEN; k++) cerr<<sigm._data[k]<<",";
    cerr<<endl;

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
    
    std::valarray<float> fMlr(nxz*LEN);
    double *ldat = G.data();
    for (int k=0; k < nxz*LEN; k++) {
        fMlr[k] = ldat[k];
    } 
    outm.put("n1",nxz);
    outm.put("n2",LEN);
    outm << fMlr;

    fsx.put("n1", LEN);
    fsx.put("n2", 1);
    fsz.put("n1", LEN);
    fsz.put("n2", 1);

    std::valarray<float> fs(LEN);
    ldat = sx.data();
    for (int k=0; k < LEN; k++) {
         fs[k] = (ldat[k]);
    } 
    fsx << fs;
    ldat = sz.data();
    for (int k=0; k < LEN; k++) {
        fs[k] = ldat[k];
    } 
    fsz << fs;

    /* ----Output Operator Matrix ---*/
    // GB
    DblNumMat smlG(smlm, LEN);
    for (int k=0; k<smlm; k++) {
	for (int j=0; j<LEN; j++) {
	    smlG(k, j) = G(smallx[k],j);
	}
    }
    DblNumMat GB(smlm, n);
    iC( ddgemm(1.0, smlG, B, 0.0, GB));
    ddat=GB.data();
    for (int k=0; k<n*smlm; k++) {
	fW[k] = ddat[k];
    }
    maxexact = fW.max();
    minexact = fW.min();
    sf_warning("GB matrix: max value=%f; min value=%f",maxexact, minexact);

    MGB.put("n1", smlm);
    MGB.put("n2", n);
    MGB.put("d1", dx);
    MGB.put("d2", dkx*twopi);
    MGB.put("o1", 0);
    MGB.put("o2", o2);
    MGB.put("label1", label1);
    MGB.put("unit1", unit1);
    MGB.put("label2", label2);
    MGB.put("unit2", unit2);
    MGB << fW;

    return 0;
}	    
















