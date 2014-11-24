/*  Weighted least square Lowrank FD coefficient on staggered grid (optimized)

   Copyright (C) 2010 University of Texas at Austin
  
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

static float twopi = 2.0*SF_PI;
static float dt, dx;
static float taper;
static float dfrq, a0; //dominant frequency

/*static float sinc(float x)
{
    if (fabs(x)<=SF_EPS) return 1.0;
    return sinf(x)/(x+SF_EPS);
    }*/

static float wghtfun(float k, float k0)
{
    return exp(-1*a0*(k-k0)*(k-k0))+exp(-1*a0*(k+k0)*(k+k0));
}

int samplex(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    float tmp = 0.0;
    float w = 0.0, k0=0.0, k=0.0;
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    tmp = sin(SF_PI*vs[rs[a]]*(ks[cs[b]])*dt);
	    if (tmp > 1-taper) {
		tmp = 1-taper;
	    }
	    else if (tmp < -1+taper) {
		tmp = -1+taper;
	    }
	    k  = twopi*ks[cs[b]];
	    k0 = twopi*dfrq/vs[rs[a]];
	    w  = wghtfun(k,k0);
	    res(a,b) = w*2.0*tmp/vs[rs[a]]/dt;
	}
    }
    return 0;
}


int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    float tmp = 0.0;
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    tmp = sin(SF_PI*vs[rs[a]]*(ks[cs[b]])*dt);
	    if (tmp > 1-taper) {
		tmp = 1-taper;
	    }
	    else if (tmp < -1+taper) {
		tmp = -1+taper;
	    }
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

    par.get("taper", taper, 0.0); //taper for stability
    
    float eps;
    par.get("eps",eps,1.e-4); // tolerance
    
    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step

    iRSF velf;
    oRSF outm;  
    oRSF fsx("sx");//, Mexactfile("Mexact"),Mlrfile("Mlr"), Mappfile("Mapp"); 

    oRSF Mexactfile("Mexact");
    oRSF Mwfun("wfun");

    float wavnumcut;
    par.get("wavnumcut",wavnumcut,1.0); // wavenumber cut percentile
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

    int COUNT = 0;
    float CUT = 0.5*nx*dk*wavnumcut;
    
    for (int ix=0; ix < nx; ix++) {
        ks[ix] = -dk*nx/2.0+ix*dk; 
	if (fabs(ks[ix]) < CUT) {
	    COUNT++;
	}
    }

    vector<int> ksc(COUNT);
    
    int nk=0, mk=0;
    // ks: -fn --- fn
    for (int ix=0; ix < nx; ix++) {
	ks[ix] = -dk*nx/2.0 + ix*dk;
	if (fabs(ks[ix])<CUT) {
	    ksc[nk] = ix;
	    nk++;
	}
    }
   
    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m, n, samplex, (double)eps, npk, cidx, ridx, mid) );
    
    sf_warning("cidx.size = %d", cidx.size());
    sf_warning("ridx.size = %d", ridx.size());

    DblNumMat M1(m, cidx.size());
    vector<int> rs(m);
    for (int k=0; k<m; k++)  rs[k]=k;
    iC( samplex(rs, cidx, M1) );
    
    DblNumMat M2(ridx.size(), n);
    vector<int> cs(n);
    for (int k=0; k<n; k++) cs[k] = k;
    iC( samplex(ridx, cs, M2) );

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
	
    nk =0;
    DblNumMat sx(len,1);
    for (int ix=0; ix<len; ix++) {
	sx._data[ix] = xtmp[ix];
    }

    cerr<<"[x]="; for(int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<"] ";    
    cerr<<endl;

    DblNumMat ktmp(1, nx); for(int k=0; k<nx; k++) ktmp._data[k]=ks[k];
    DblNumMat ktmpc(1,COUNT); for(int k=0; k<COUNT; k++) ktmpc._data[k]=ks[ksc[k]];
    
    int LEN = sx._m;
    DblNumMat B(LEN,nx);
    DblNumMat Bc(LEN,COUNT);
    iC(ddgemm(2*SF_PI*dx,sx,ktmp,0.0,B));
    iC(ddgemm(2*SF_PI*dx,sx,ktmpc,0.0,Bc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=sin(B._data[k]);
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=sin(Bc._data[k]);
    
    int itmp;
    float w = 0.0, k0=0.0, kk=0.0;
    std::valarray<float> wfun(B._n);
    
    for (int ik=0; ik<B._n; ik++) {
	kk = twopi*ktmp._data[ik];
	k0 = twopi*dfrq/vs[ik];
	w  = wghtfun(kk, k0);
	//sf_warning("kk=%f", kk);
	wfun[ik] = w;
	for (int il=0; il<B._m; il++) {
	    itmp = ik*B._m+il;
	    B._data[itmp] = w*B._data[itmp];
	}
    }
	    
    for (int ik=0; ik<Bc._n; ik++) {
	kk = twopi*ktmpc._data[ik];
	k0 = twopi*dfrq/vs[ik];
	w  = wghtfun(kk, k0);
	for (int il=0; il<Bc._m; il++) {
	    itmp = ik*Bc._m+il;
	    Bc._data[itmp] = w*Bc._data[itmp];
	}
    }

    DblNumMat IB(nx,LEN);    iC( ddpinv(B, 1e-16, IB) );
    DblNumMat IBc(COUNT,LEN);    iC( ddpinv(Bc, 1e-16, IBc) );
          
    DblNumMat coef(ridx.size(),LEN);
    DblNumMat M2c;
    iC( samplex(ridx,ksc,M2c) );
    iC(ddgemm(1.0,M2c,IBc,0.0,coef));
 
    DblNumMat G(nx,LEN), tmpG(mid._m,LEN);
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,G));
        
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
  
    // Exact matre
    DblNumMat Mexact(nx,nx);
    iC( sample(rs,cs,Mexact) );
       
    //float dk2=dk/2.0;
    Mexactfile.put("n1",nx);
    Mexactfile.put("n2",nx);
    Mexactfile.put("d2",dk);
    Mexactfile.put("o2",0);
    
    std::valarray<float> fMex(nx*nx);
    ldat = Mexact.data();
    for (int k=0; k < nx*nx; k++) {
        fMex[k] = ldat[k];
    } 
    Mexactfile << fMex;

    //wfun
    Mwfun.put("n1",B._n);
    Mwfun.put("d1",dk);
    float o1 = -1*dk*B._n/2.0;
    Mwfun.put("o1",o1);
    Mwfun << wfun ;

    sf_warning("==========================");
    sf_warning("Weighted LS LFD");
    sf_warning("COUNT=%d", COUNT);
    sf_warning("nk=%d",nk);
    sf_warning("mk=%d",mk);
    sf_warning("taper=%f", taper);
    sf_warning("k0=%f", k0);
    sf_warning("==========================");
    
    return 0;
}
    
    
    
    
    
    
    
