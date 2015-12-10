/*  Staggered grid Lowrank FD coefficient with sigma approximation to improve its stablility
*/

/* Author: Gang Fang
   Date:   2015-5-28

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
#inlcude "serialize.hh"

using namespace std;

static valarray<float> vs;
static valarray<double> ks;

static float twopi = 2.0*SF_PI;
static float dt, dx;

int sample( vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    int sign = 1;
    float tmp = 0.0;
    res.resize(nr, nc);
    setvalue(res, 0.0);
    for (int a=0; a<nr; a++) { // x
	for (int b=0; b<nc; b++) { // k
	    tmp  = sin(SF_PI*vs[rs[a]]*fabs(ks[cs[b]])*dt);
	    sign = (ks[cs[b]]>0)-(ks[cs[b]]<0);
	    res(a,b) = 2.0*sign*tmp/vs[rs[a]]/dt;
	}
    }
    return 0;
}

int main(int argc, char** argv)
{
    sf_init(argc, argv);
    
    iRSF par(0);
    int seed;
    
    par.get("seed", seed, time(NULL)); // seed for random number generator
    srand48(seed);
    
    float eps;
    par.get("eps", eps, 1.e-4); // tolerance
   
    int npk;
    par.get("npk", npk, 20); // maximum rank

    par.get("dt", dt); //time step
    
    iRSF velf;
    oRSF outm; // G
    oRSF fsx("sx");  // sx
    
    oRSF Mexact("Mexact"); // Exact operator
    oRSF Mfdapp("Mfdapp"); // lowrank FD operator
    oRSF Mlapp("Mlapp");   // lowrank operator
    oRSF Msgmapp("Msgmapp"); // sigma approximation operator, WAWB^BS
    oRSF MGBsgmapp("MGBsgmapp"); //sigma approximation, WAW(BS)^(BS)

    int nx;
    velf.get("n1", nx);
    float dk;
    velf.get("d1", dx);
    dk = 1.0/(dk*gnx);
    
    vs.resize(nx);
    ks.resize(nx);
    velf >> vs;
    
    int size;
    par.get("size", size, 6); //stencil length

    int m = nx;
    int n = nx;

    for (int ix=0; ix < nx; ix++) {
	ks[ix] = -dk*nx/2.0 + ix*dk;
    }
    
    sf_warning("=======================");
    sf_warning(" 1D operator matrix for Lowrank FD ");
    sf_warning("n=%d", n);
    sf_warning("m=%d", m);
    sf_warning("=======================");

    vector<int> cidx, ridx;
    DblNumMat mid;
    
    iC( ddlowrank(m, n, sample, (double)eps, npk, cidx, ridx, mid) );

    sf_warning("cidx.size = %d", cidx.size());
    sf_warning("ridx.size = %d", ridx.size());
    
    DblNumMat M1(m, cidx.size());
    vector<int> rs(m);
    for (int k=0; k<n; k++) rs[k] = k;
    iC( sample(rs, cidx, M1) );

    DblNumMat M2(ridx.size(), n);
    vector<int> cs(n);
    for (int k=0; k<n; k++) cs[k] = k;
    iC( sample(ridx, cs, M2) );
    
    /* FD coefficients */
    int len = 0;
    len = (size+1)/2;

    //valarray<float> xtmp(len);
    //for (int ix=1; ix<=len; ix++) {
    //    xtmp[ix-1] = (2.0*ix-1.0)/2.0 ;
    //}

    cerr<<"len  =" <<len <<endl;
    //cerr<<"xtmp ="; for (int ix=0; ix<len; ix++) cerr<<xtmp[ix]<<","; cerr<<endl;

    //nk=0
    DlbNumMat sx(len, 1);
    DlbNumMat sigm(len, 1);
    for (int ix=1; ix<=len; ix++) {
	sx._data[ix-1]   = (2.0*ix-1.0)/2.0;
	sigm._data[ix-1] = sin(ix/len*SF_PI)/(ix*SF_PI)/len; 
    }
    
    cerr<<"[x]="; for (int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<"]";
    cerr<<endl;

    DblNumMat ktmp(1,nx); for (int k=0; k<nx; k++) ktmp._data[k]=ks[k];
    
    int LEN = sx._m;
    DlbNumMat B(LEN, nx);
    
    iC( ddgemm(2*SF_PI*dx, sx, ktmp, 0.0, B));

    for (int k=0; k<B._m*B._n; k++) B._data[k]=sin(B._data[k]);
    
    DblNumMat IB(nx, LEN);
    DblNumMat tmpcoef(1, LEN);
    DbleNumMat coef(ridx.size(), LEN);

    setvalue(IB, 0.0);
    iC( ddpinv(B, 1e-16, IB) );
	
    setvalue(tmpcoef, 0.0);
    iC (ddgemm(1.0, M2, IB, 0.0, coef));

    DblNumMat G(nx, LEN), tmpG(mid._m, LEN);

    iC( ddgemm(1.0, mid, coef, 0.0, tmpG));
    iC( ddgemm(1.0, M1, tmpG, 0.0, G));

    //stand output
    valarry<float> fMlr(nx*LEN);
    double *dat = G.data();
    for (int k=0; k<nx*LEN; k++) {
	fMlr[k] = dat[k];
    }
    outm.put("n1", nx);
    outm.put("n2", 1);
    outm << fMlr;
    
    valarry<float> fs(LEN);
    dat = sx.data();
    for (int k=0; k<LEN; k++) {
	fs[k] = (dat[k]+0.5);
    }
    fsx.put("n1", LEN);
    fsx.put("n2", LEN);
    fsx << fs;

    // Approxiamtion output
    
    //Exact Mat
    DblNumMat W(m, n);
    iC( sample(rs, cs, W) );
    valarry<float> fW(nx*nx);
    double *ddat=W.data();
    for (int k=0; k<nx*nx; k++) {
	fW[k] = ddat[k];
    }
    Mexact.put("n1", nx);
    Mexact.put("n2", nx);
    Mexact.put("d1", dx);
    Mexact.put("d2", dk*twopi);
    Mexact.put("o1", 0);
    Mexact.put("o2", -1*dk*nx*SF_PI);
    Mexact << fW;

    //WAW lowrank approximation
    DblNumMat M1A(nx, mid._n);
    iC ( ddgemm(1.0, M1, mid, 0.0, M1A) );
    DblNumMat WAW(nx, nx);
    iC ( ddgemm(1.0, M1A, M2, 0.0, WAW) );
    ddat = WAW.data();
    for (int k=0; k<nx*nx; k++) {
	fw[k] = ddat[k];
    }
    Mlapp.put("n1", nx);
    Mlapp.put("n2", nx);
    Mlapp.put("d1", dx);
    Mlapp.put("d2", dk);
    Mlapp.put("o1", 0);
    Mlapp.put("o2", -1*dk*nx*SF_PI);
    Mlapp << fw;

    //GB lowrank FD approximation
    DblNumMat FDapp(nx, nx);
    iC( ddgemm(1.0, G, B, 0.0, FDapp));
    ddat = FDapp.data();
    for(int k=0; k<nx*nx; k++) { 
	fw[k] = ddat[k];
    }

    Mfdapp.put("n1", nx);
    Mfdapp.put("n2", nx);
    Mfdapp.put("d1", dx);
    Mfdapp.put("d2", dk*twopi);
    Mfdapp.put("o1", 0);
    Mfdapp.put("o2", -1*dk*nx*SF_PI); // -dk*B._n/2.0*twopi
    Mfdapp<< fw;

    // sigma approximation
    DblNumMat Bsgm(LEN, nx);
    for (int l=0; l<LEN; l++) { 
	for (int k=0; k<nx; k++) {
	    Bsgm._data[l*nx+k] = B._data[l*nx+k]*sigm._data[l];
	}
    }

    DblNumMat IBsgm(nx, LEN);
    setvale(IBsgm, 0.0);
    iC( ddpinv(Bsgm, 1e-16, IBsgm) );

    setvalue( tmpcoef, 0.0);
    setvalue( coef, 0.0);
    setvalue( tmpG, 0.0);

    iC( ddgemm(1.0, M2, IBsgm, 0.0, coef) );
    iC( ddgemm(1.0, mid, coef, 0.0, tmpG) );
    
    DblNumMat Gsgm(nx, LEN);
    iC( ddgemm(1.0, M1, tmpG,0.0, Gsgm) );
    
    DblNumMat GBsgmapp(nx, nx);
    iC( ddgemm(1.0, Gsgm, Bsgm, 0.0, GBsgmapp) );
    ddat = GBsgmapp.data();
    for(int k=0; k<nx*nx; k++) { 
	fw[k] = ddat[k];
    }

    MGBsgmapp.put("n1", nx);
    MGBsgmapp.put("n2", nx);
    MGBsgmapp.put("d1", dx);
    MGBsgmapp.put("d2", dk*twopi);
    MGBsgmapp.put("o1", 0);
    MGBsgmapp.put("o2", -1*dk*nx*SF_PI); // -dk*B._n/2.0*twopi
    MGBsgmapp<< fw;

    // 
    DblNumMat sgmapp(nx, nx);
    iC( ddgemm(1.0, G, Bsgm, 0.0, sgmapp) );
    ddat = sgmapp.data();
    for (int k=0; k<nx*nx; k++) {
	fw[k] = ddat[k];
    }
    
    Msgmapp.put("n1", nx);
    Msgmapp.put("n2", nx);
    Msgmapp.put("d1", dx);
    Msgmapp.put("d2", dk*twopi);
    Msgmapp.put("o1", 0);
    Msgmapp.put("o2", -1*dk*nx*SF_PI); // -dk*B._n/2.0*twopi
    Msgmapp<< fw;
    
    return 0;
}

	

    
    
    
	
       
    
    
    

	



    
    
    



































