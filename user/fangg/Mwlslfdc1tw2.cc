/*  Weighted least square Lowrank FD coefficient on staggered grid (optimized)
*/

/* Copyright (C) 2010 University of Texas at Austin
  
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

#include <rsf.hh>

#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

static std::valarray<float> vs;
static std::valarray<double> ks;

static float pi=SF_PI;
static float twopi = 2.0*SF_PI;
static float dt, dx;
static float tpa, tpb;
static float dfrq, a0; //dominant frequency

float wghtfun(float k, float k0)
{
    return 0.5*(exp(-1*a0*(k-k0)*(k-k0))+exp(-1*a0*(k+k0)*(k+k0)));
    //return exp(-1*a0*(fabs(k)-k0)*(fabs(k)-k0));
}

float tpfun(float a, float b, float k, float k0) 
{
    float res   = 0.0;
    float halfa = 0.5*a; 
    float kkh, kkh2;
    if (fabs(k-k0) < halfa) {
	kkh  = (k-k0)/halfa;
	kkh2 = kkh*kkh;
	res  = 1.0 - b*exp(1.0-1.0/(1.0-kkh2));
    } 
    else if (fabs(k+k0) < halfa ) {
	kkh  = (k+k0)/halfa;
	kkh2 = kkh*kkh;
	res  = 1.0 - b*exp(1.0-1.0/(1.0-kkh2));
    } else {
	res  = 1.0;
    }
    return res;
}	


/*int sampletp(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    float tmp=0.0, k=0.0, k0=0.0;
    float tp=0.0;
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) { // x 
	k0 = 0.5*pi/(vs[rs[a]]*dt*0.5);
	for(int b=0; b<nc; b++) { // k
	    k  = twopi*ks[cs[b]];
	    tp = tpfun(tpa, tpb, k, k0);
	    //sf_warning("tp=%f, k0=%f, k=%f", tp, k0, k);
	    tmp = tp*sin(pi*vs[rs[a]]*(ks[cs[b]])*dt);
	    res(a,b) =2.0*tmp/vs[rs[a]]/dt;
	}
    }
    return 0;
    }*/

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    float tmp = 0.0;
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) { // x 
	for(int b=0; b<nc; b++) { // k
	    tmp = sin(pi*vs[rs[a]]*(ks[cs[b]])*dt);
	    res(a,b) =2.0*tmp/vs[rs[a]]/dt;
	}
    }
    return 0;
}



int main(int argc, char** argv)
{
    sf_init(argc, argv);
    
    iRSF par(0);
    int seed;

    par.get("seed", seed, time(NULL));// seed for random number generator
    srand48(seed);

    par.get("f0", dfrq, 15); //dominant frequency
    
    par.get("a0", a0, 0.0001); // weight parameters

    bool taper;
    par.get("taper",taper, true);
    par.get("tpa", tpa, 0.0); //taper for stability
    par.get("tpb", tpb, 0.0); 
    
    float eps;
    par.get("eps",eps,1.e-4); // tolerance
    
    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step
    
    bool weight;
    par.get("weight", weight, true);
        
    iRSF velf;
    oRSF outm;  
    oRSF fsx("sx");//, Mexactfile("Mexact"),Mlrfile("Mlr"), Mappfile("Mapp"); 
    
    oRSF Mwatpwfile("Mwatpw");
    oRSF Mwfun("wfun");
    oRSF Mappfile("Mapp");
    oRSF Mtpfile("tp");
    oRSF Mwawfile("waw");

    int nx;
    velf.get("n1",nx);
    float dk;
    velf.get("d1",dx);
    dk = 1.0/(dx*nx);

    vs.resize(nx);
    ks.resize(nx);
    velf >> vs;
    
    int size;
    par.get("size",size,6); // stencil length 
    
    int m = nx;
    int n = nx;

    for (int ix=0; ix < nx; ix++) {
        ks[ix] = -dk*nx/2.0+ix*dk; 
    }

    sf_warning("==========================");
    sf_warning("Weighted LS LFD");
    //sf_warning("LEN=%d", sx._m);
    sf_warning("weight=%d",weight);
    sf_warning("taper=%d",taper);
    sf_warning("n=%d",n);
    sf_warning("m=%d",m);
    sf_warning("tpa=%f", tpa);
    sf_warning("tpb=%f", tpb);
    sf_warning("dfrq=%f", dfrq);
    sf_warning("==========================");

    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m, n, sample, (double)eps, npk, cidx, ridx, mid) );
    
    sf_warning("cidx.size = %d", cidx.size());
    sf_warning("ridx.size = %d", ridx.size());
    
    DblNumMat M1(m, cidx.size());
    vector<int> rs(m);
    for (int k=0; k<m; k++)  rs[k]=k;
    iC( sample(rs, cidx, M1) );

    DblNumMat M2(ridx.size(), n);
    vector<int> cs(n);
    for (int k=0; k<n; k++) cs[k] = k;
    iC( sample(ridx, cs, M2) );
    
    /* FD coefficient*/

    /* d/dx */
    int len = 0;
    len = (size+1)/2; //x: -len ,... -1, 1, ... len;
    std::valarray<float> xtmp(len);
    for (int ix = 1; ix<=len; ix++) {
	xtmp[ix-1] = (2.0*ix-1.0)/2.0;     //   (2l-1)/2
    }

    cerr<<"len  = " <<len <<endl;
    cerr<<"xtmp = "; for (int ix=0; ix<len; ix++) cerr<<xtmp[ix]<<", "; cerr<<endl;
    
    //nk =0;
    DblNumMat sx(len,1);
    for (int ix=0; ix<len; ix++) {
	sx._data[ix] = xtmp[ix];
    }
    
    cerr<<"[x]="; for(int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<"] ";    
    cerr<<endl;
    
    DblNumMat ktmp(1, nx); for(int k=0; k<nx; k++) ktmp._data[k]=ks[k];

    int LEN = sx._m;
    DblNumMat B(LEN,nx);
    
    iC(ddgemm(2*pi*dx,sx,ktmp,0.0,B));

    for(int k=0; k<B._m*B._n; k++) B._data[k]=sin(B._data[k]);
    
    int itmp;
    float w=0.0, k0=0.0, kk=0.0;
    float tp=0.0, tpk0=0.0, tpk=0.0;
    std::valarray<float> wfun(B._n*ridx.size());
    
    DblNumMat wtB(LEN, nx);
    DblNumMat tmpcoef(1, LEN);
    DblNumMat wtM2(1, nx), tpM2(ridx.size(), n);
    DblNumMat IB(nx,LEN);
    DblNumMat coef(ridx.size(),LEN);

    for (int ixm=0; ixm<(int)ridx.size(); ixm++) {
	k0 = twopi*dfrq/vs[ridx[ixm]];
	tpk0 = 0.5*pi/(vs[ridx[ixm]]*dt*0.5);
	sf_warning("ixm=%d, ridx[ixm]=%d, vs=%f", ixm, ridx[ixm], vs[ridx[ixm]]);
	for (int ik=0; ik<nx; ik++) {
	    kk = twopi*ktmp._data[ik];
	    if (weight == true ) {
		w  = wghtfun(kk, k0);
		/*weight*/
	    } else {
		w = 1.0;
	    }
	    tpk = twopi*ks[ik];

	    if (taper == true ) {
		tp  = tpfun(tpa, tpb, tpk, tpk0);
	    } else {
		tp = 1.0;
	    }

	    wtM2(0,ik) = tp*w*M2(ixm, ik);
	    tpM2(ixm,ik) = tp*M2(ixm, ik);
	    /*weight*M2*/
	    wfun[ik+ixm*B._n] = w;
	    for (int il=0; il<B._m; il++) {
		itmp = ik*B._m+il;
		wtB._data[itmp] = w*B._data[itmp];
		/* Weight*B */
	    }
	}
	setvalue(IB, 0.0);
	iC( ddpinv(wtB, 1e-16, IB) );
	
	setvalue(tmpcoef, 0.0);
	iC( ddgemm(1.0,wtM2,IB,0.0,tmpcoef) );
	
	for (int il=0; il<LEN; il++) {
	    coef(ixm, il) = tmpcoef(0, il);
	}
    }
    
    DblNumMat G(nx,LEN), tmpG(mid._m,LEN);
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,G));

    // reconstruct propagation metrix
    DblNumMat M1A(nx, mid._n);
    iC( ddgemm(1.0, M1, mid, 0.0, M1A));
    
    // Exact matre
    DblNumMat Mwatpw(m,n);
    iC( ddgemm(1.0, M1A, tpM2, 0.0, Mwatpw));

    Mwatpwfile.put("n1",nx);
    Mwatpwfile.put("n2",nx);
    Mwatpwfile.put("d2",dk);
    Mwatpwfile.put("o2",0);

    std::valarray<float> fMex(nx*nx);
    double *tmpdat;
    tmpdat = Mwatpw.data();
    for (int k=0; k < nx*nx; k++) {
        fMex[k] = tmpdat[k];
    } 
    Mwatpwfile << fMex;

    DblNumMat Mwaw(m,n);
    iC( ddgemm(1.0, M1A, M2, 0.0, Mwaw));
    Mwawfile.put("n1",nx);
    Mwawfile.put("n2",nx);
    Mwawfile.put("d2",dk);
    Mwawfile.put("o2",0);
    
    tmpdat = Mwaw.data();
    for (int k=0; k < nx*nx; k++) {
        fMex[k] = tmpdat[k];
    } 
    Mwawfile << fMex;

    
    std::valarray<float> fMlr(nx*LEN);
    double *ldat = G.data();
    for (int k=0; k < nx*LEN; k++) {
        fMlr[k] = ldat[k];
    } 
    outm.put("n1",nx);
    outm.put("n2",LEN);
    outm << fMlr;
    
    fsx.put("n1", LEN);
    fsx.put("n2", 1);

    std::valarray<float> fs(LEN);
    ldat = sx.data();
    for (int k=0; k < LEN; k++) {
	fs[k] = (ldat[k]+0.5);
    } 
    fsx << fs;

    
    DblNumMat Mapp(nx, nx);
    std::valarray<float> fwtfull(nx*nx);
    iC( ddgemm(1.0,G,B,0.0,Mapp) );

    Mappfile.put("n1",nx);
    Mappfile.put("n2",nx);
    Mappfile.put("d2",dk*twopi);
    Mappfile.put("o1",0);
    float o2 = -1*dk*B._n/2.0;
    Mappfile.put("o2",o2*twopi);
    
    ldat = Mapp.data();
    for (int k=0; k < nx*nx; k++) {
        fMex[k] = ldat[k];
    } 
    Mappfile << fMex;
    
    //wfun3
    Mwfun.put("n1",B._n);
    Mwfun.put("n2",(int)ridx.size());
    Mwfun.put("d1",dk*twopi);
    Mwfun.put("d2",dx);
    float o1 = -1*dk*B._n/2.0;
    Mwfun.put("o1",o1*twopi);
    Mwfun << wfun ;


    //taper function
    std::valarray<float> Mtp(nx*nx);
    for (int ix=0; ix<nx; ix++) {
	k0 = 0.5*pi/(vs[ix]*dt*0.5);
	for (int ik=0; ik<nx; ik++) {
	    kk  = twopi*ks[ik];
	    tp = tpfun(tpa, tpb, kk, k0);
	    Mtp[ix*nx+ik] =tp;
	}
    }

    Mtpfile.put("n1",nx);
    Mtpfile.put("n2",nx);
    Mtpfile.put("d1",twopi*dk);
    Mtpfile.put("d2",dx);
    Mtpfile.put("o1",o1*twopi);
    Mtpfile << Mtp ;

    return 0;
}

    
    
    
    

    
    



























