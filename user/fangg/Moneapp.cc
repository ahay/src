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
#include "serialize.hh"

using namespace std;

static valarray<float> vs;
static valarray<double> ks, kss;

static float twopi = 2.0*SF_PI;
static float dt, dx;
static float tpa, tpb;
static float dfrq, a0; //dominant frequency

static float sinc(float x)
{
    if (fabs(x)<=SF_EPS) return 1.0;
    return sinf(x)/(x+SF_EPS);
}


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
	    //tmp  = sin(SF_PI*vs[rs[a]]*fabs(ks[cs[b]])*dt);
	    //sign = (ks[cs[b]]>0)-(ks[cs[b]]<0);
	    //res(a,b) = 2.0*sign*tmp/vs[rs[a]]/dt;
	    tmp  = sin(SF_PI*vs[rs[a]]*fabs(kss[cs[b]])*dt/dx);
	    sign = (kss[cs[b]]>0)-(kss[cs[b]]<0);
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

    bool weight;
    par.get("weight", weight, true);
    par.get("f0", dfrq, 15); //dominant frequency
    par.get("a0", a0, 0.0001); // weight parameters

    bool taper;
    par.get("taper",taper, true);
    par.get("tpa", tpa, 0.0); //taper for stability
    par.get("tpb", tpb, 0.0);
    
    float eps;
    par.get("eps", eps, 1.e-4); // tolerance
   
    int npk;
    par.get("npk", npk, 20); // maximum rank

    par.get("dt", dt); //time step

    //bool rescale;
    //par.get("rescale", rescale, false); // rescale the sigma approximation
    int filter;
    par.get("filter", filter, 0);
    /*0- sigma filter; 1- cos filter; 2-expontial filter*/
    
    float p;
    par.get("p", p, 0.1);
    /* power of sigma factor, 0<p<=1 */

    bool disp;  
    par.get("dispersion", disp, false);
    /* Dispersion */

    float wavnumcut;
    par.get("wavnumcut", wavnumcut, 0.75); // wavenumber cut percentile
    
    iRSF velf;
    oRSF outm; // G
    oRSF fsx("sx");  // sx
    oRSF fsigm("sigm");
    
    oRSF Mexact("Mexact"); // Exact operator
    oRSF Mfdapp("Mfdapp"); // lowrank FD operator
    oRSF Mlapp("Mlapp");   // lowrank operator
    oRSF Msgmapp("Msgmapp"); // sigma approximation operator, WAWB^BS
    oRSF MGBsgmapp("MGBsgmapp"); //sigma approximation, WAW(BS)^(BS)
    oRSF MsgmG("MsgmG"); // FD coefficients with sigma approximation
    oRSF MwtG("MwtG"); // FD coefficients get from weight and taper 
    oRSF Mwtapp("Mwtapp"); // WLS and taper approximation
    oRSF MwsG("MwsG");     // WLS and sigma approximation
    oRSF Mwsgmapp("Mwsgmapp"); //WLS with sigma approximation


    int nx;
    velf.get("n1", nx);
    float dk, dkk;
    velf.get("d1", dx);
    dk = 1.0/(dx*nx);
    dkk = dk*dx;
    vs.resize(nx);
    ks.resize(nx);
    kss.resize(nx);
    velf >> vs;

    //for (int k=0; k<nx; k++) {
    //sf_warning("v[%d]=%f", k, vs[k]);
    //}
    
    int size;
    par.get("size", size, 6); //stencil length

    bool norm;
    par.get("norm", norm, true);

    int m = nx;
    int n = nx;

    sf_warning("=======================");
    sf_warning(" 1D operator matrix for Lowrank FD ");
    sf_warning("weight=%d",weight);
    sf_warning("taper=%d",taper);
    sf_warning("wavnumcut=%f", wavnumcut);
    sf_warning("n=%d",n);
    sf_warning("m=%d",m);
    sf_warning("tpa=%f", tpa);
    sf_warning("tpb=%f", tpb);
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
	    break;
	default :
	    sf_warning("filter=%d, close filter", filter);
    } 
    sf_warning("sigm_p=%f", p);
    sf_warning("=======================");

    for (int ix=0; ix < nx; ix++) {
	ks[ix] = -dk*nx/2.0 + ix*dk;
	kss[ix] = ks[ix]*dx;
    }

    vector<int> cidx, ridx;
    DblNumMat mid;

    //count
    int COUNT = 0;
    float CUT = 0.5*nx*dk*wavnumcut;
    int SMK = 0;
    
    for (int ix=0; ix<nx; ix++) {
	if (fabs(ks[ix]) < CUT ) {
	    COUNT++;
	}
	if ( fabs(ks[ix] )<( dk+0.00001 )) SMK++;
    }
    
    vector<int> ksc(COUNT), smallk(SMK);
    int nk=0, mk=0;
    
    for (int ix=0; ix<nx; ix++) {
	if ( fabs(ks[ix]) < CUT ) {
	    ksc[nk] = ix ;
	    nk++;
	}
	if ( fabs(ks[ix]) < (dk+0.00001)) {
		smallk[mk] = ix;
		mk++;
	}
    }

    
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
    //len = (size+1)/2;
    len = size;

    //valarray<float> xtmp(len);
    //for (int ix=1; ix<=len; ix++) {
    //    xtmp[ix-1] = (2.0*ix-1.0)/2.0 ;
    //}

    cerr<<"len  =" <<len <<endl;
    //cerr<<"xtmp ="; for (int ix=0; ix<len; ix++) cerr<<xtmp[ix]<<","; cerr<<endl;

    //nk=0
    DblNumMat sx(len, 1);
    DblNumMat sigm(len, 1);
    DblNumMat sigmp(len, 1);
    for (int ix=1; ix<=len; ix++) {
	sx._data[ix-1]   = (2.0*ix-1.0)/2.0;
    }
    
    switch (filter) {
	case 1 : // sigma filter
	    for (int ix=0; ix<len; ix++) {
		//sigm._data[ix-1] = sin(ix*SF_PI/len)/(ix*SF_PI/len); 
		sigm._data[ix] = fabs(sinc((fabs(sx._data[ix])+0.5)*SF_PI/len));
	    }
	    break;
	case 2 : // cose filter
	    for (int ix=0; ix<len; ix++) {
		//sigm._data[ix-1] = 0.5*(1+cos(ix*SF_PI/len));
		sigm._data[ix] = 0.5*(1+cos((fabs(sx._data[ix])+0.5)*SF_PI/len));
	    }
	    break;
	case 3 : // exponential filter
	    double alpha;
	    alpha = -1.0*log(1.0e-15);
	    sf_warning("ALPHA=%g", alpha);
	    for (int ix=0; ix<len; ix++) {
		//sigm._data[ix-1] = exp(-1.0*alpha*pow(ix/len,2.0));
		sigm._data[ix] = exp(-1.0*alpha*pow((fabs(sx._data[ix])+0.5)/len,4.0));
	    }
	    break;
	default :
	   for (int ix=0; ix<len; ix++) {
		sigm._data[ix] = 1.0; 
	    }
    } 
    
    for(int ix=0; ix<len; ix++) {
	sigmp._data[ix] = pow(sigm._data[ix], p);
    }

    cerr<<"[x]="; for (int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<"]";
    cerr<<endl;

    DblNumMat ktmp(1,nx);      for (int k=0; k<nx; k++) ktmp._data[k]=ks[k];
    DblNumMat ktmpc(1, COUNT); for (int k=0; k<COUNT; k++) ktmpc._data[k]=ks[ksc[k]];

    int LEN = sx._m;
    DblNumMat B(LEN, nx);
    DblNumMat Bc(LEN, COUNT);
    
    //iC( ddgemm(2*SF_PI*dx, sx, ktmp, 0.0, B));
    iC( ddgemm(2.0*SF_PI*dx, sx, ktmp, 0.0, B));
    iC( ddgemm(2.0*SF_PI*dx, sx, ktmpc, 0.0, Bc));

    for (int k=0; k<B._m*B._n; k++) B._data[k]=sin(B._data[k]);
    for (int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=sin(Bc._data[k]);

    DblNumMat IB(nx, LEN);     iC( ddpinv(B,  1e-16, IB));
    DblNumMat IBc(COUNT, LEN); iC( ddpinv(Bc, 1e-16, IBc));
    DblNumMat coef(ridx.size(), LEN);
    DblNumMat M2c(ridx.size(), COUNT);
    iC(sample(ridx, ksc, M2c) );
    
    iC (ddgemm(1.0, M2c, IBc, 0.0, coef));

    DblNumMat G(nx, LEN), tmpG(mid._m, LEN);

    iC( ddgemm(1.0, mid, coef, 0.0, tmpG));
    iC( ddgemm(1.0, M1, tmpG, 0.0, G));

    //stand output G
    valarray<float> fMlr(nx*LEN);
    double *dat = G.data();
    for (int k=0; k<nx*LEN; k++) {
	fMlr[k] = dat[k];
    }
    outm.put("n1", nx);
    outm.put("n2", LEN);
    outm << fMlr;
    
    valarray<float> fs(LEN);
    dat = sx.data();
    for (int k=0; k<LEN; k++) {
	fs[k] = (dat[k]+0.5);
    }
    fsx.put("n1", LEN);
    fsx.put("n2", 1);
    fsx << fs;
    
    dat = sigmp.data();
    for (int k=0; k<LEN; k++) {
	fs[k] = dat[k];
    }
    fsigm.put("n1", LEN);
    fsigm.put("n2", 1);
    fsigm << fs;
    
    // Approxiamtion output
    const char *label1 = "Distance";
    const char *label2 = "Wavenumber";
    const char *unit1 = "km" ; 
    const char *unit2 = "1/km" ;
    
    //Exact Mat
    DblNumMat W(m, n);
    iC( sample(rs, cs, W) );
    valarray<float> fW(nx*nx);
    double *ddat=W.data();
    int sign = 1;
    if (norm) {
	for (int k=0; k<nx; k++) {
	    sign = (kss[k]>0)-(kss[k]<0);
	    for (int ix=0; ix<nx; ix++) {
		fW[k*nx+ix] = ddat[k*nx+ix]*vs[ix]*dt*0.5*sign;
	    } 
	}
    } else {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = ddat[k];
	}
    }

    float maxexact, minexact;
    maxexact = fW.max();
    minexact = fW.min();

    sf_warning("Exact matrix: max value=%f; min value=%f",maxexact, minexact);

    Mexact.put("n1", nx);
    Mexact.put("n2", nx);
    Mexact.put("d1", dx);
    Mexact.put("d2", dk*twopi);
    Mexact.put("o1", 0);
    float o2 = -1*dk*nx*SF_PI ;  
    Mexact.put("o2", o2);
    Mexact.put("label1", label1);
    Mexact.put("unit1", unit1);
    Mexact.put("label2", label2);
    Mexact.put("unit2", unit2);
    Mexact << fW;

    //WAW lowrank approximation
    DblNumMat M1A(nx, mid._n);
    iC ( ddgemm(1.0, M1, mid, 0.0, M1A) );
    DblNumMat WAW(nx, nx);
    iC ( ddgemm(1.0, M1A, M2, 0.0, WAW) );
    ddat = WAW.data();
    if (norm) {
	for (int k=0; k<nx; k++) {
	    sign = (kss[k]>0)-(kss[k]<0);
	    for (int ix=0; ix<nx; ix++) {
		fW[k*nx+ix] = ddat[k*nx+ix]*vs[ix]*dt*0.5*sign;
	    } 
	}
    } else {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = ddat[k];
	}
    }
    
    Mlapp.put("n1", nx);
    Mlapp.put("n2", nx);
    Mlapp.put("d1", dx);
    Mlapp.put("d2", dk*twopi);
    Mlapp.put("o1", 0);
    Mlapp.put("o2", o2);
    Mlapp.put("label1", label1);
    Mlapp.put("unit1", unit1);
    Mlapp.put("label2", label2);
    Mlapp.put("unit2", unit2);
    Mlapp << fW;

    //GB lowrank FD approximation
    DblNumMat FDapp(nx, nx);
    iC( ddgemm(1.0, G, B, 0.0, FDapp));
    ddat = FDapp.data();
    if (norm) {
	for (int k=0; k<nx; k++) {
	    sign = (kss[k]>0)-(kss[k]<0);
	    for (int ix=0; ix<nx; ix++) {
		fW[k*nx+ix] = ddat[k*nx+ix]*vs[ix]*dt*0.5*sign;
	    } 
	}
    } else {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = ddat[k];
	}
    }
   
    Mfdapp.put("n1", nx);
    Mfdapp.put("n2", nx);
    Mfdapp.put("d1", dx);
    Mfdapp.put("d2", dk*twopi);
    Mfdapp.put("o1", 0);
    Mfdapp.put("o2", o2); // -dk*B._n/2.0*twopi
    Mfdapp.put("label1", label1);
    Mfdapp.put("unit1", unit1);
    Mfdapp.put("label2", label2);
    Mfdapp.put("unit2", unit2);
    Mfdapp<< fW;

    //GB sigma approximation, WAW(BS)^(BS)
    DblNumMat Bcsgm(LEN, COUNT);
    DblNumMat Bsgm(LEN, nx);
   
    for (int k=0; k<COUNT; k++) {
	for (int l=0; l<LEN; l++) { 
	    Bcsgm._data[k*LEN+l] = Bc._data[k*LEN+l]*sigmp._data[l];
	}
    }

    for (int k=0; k<nx; k++) {
	for (int l=0; l<LEN; l++) { 
	    Bsgm._data[k*LEN+l] = B._data[k*LEN+l]*sigmp._data[l];
	}
    }

    DblNumMat IBcsgm(COUNT, LEN);
    iC( ddpinv(Bcsgm, 1e-16, IBcsgm) );

    setvalue( coef, 0.0);
    setvalue( tmpG, 0.0);

    iC( ddgemm(1.0, M2c, IBcsgm, 0.0, coef) );
    iC( ddgemm(1.0, mid, coef, 0.0, tmpG) );
        
    DblNumMat Gsgm(nx, LEN);
    iC( ddgemm(1.0, M1, tmpG,0.0, Gsgm) );

    DblNumMat GBsgmapp(nx, nx);
    iC( ddgemm(1.0, Gsgm, Bsgm, 0.0, GBsgmapp) );

    ddat = GBsgmapp.data();
    if (norm) {
	for (int k=0; k<nx; k++) {
	    sign = (kss[k]>0)-(kss[k]<0);
	    for (int ix=0; ix<nx; ix++) {
		fW[k*nx+ix] = ddat[k*nx+ix]*vs[ix]*dt*0.5*sign;
	    } 
	}
    } else {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = ddat[k];
	}
    }
    
    MGBsgmapp.put("n1", nx);
    MGBsgmapp.put("n2", nx);
    MGBsgmapp.put("d1", dx);
    MGBsgmapp.put("d2", dk*twopi);
    MGBsgmapp.put("o1", 0);
    MGBsgmapp.put("o2", o2); // -dk*B._n/2.0*twopi
    MGBsgmapp.put("label1", label1);
    MGBsgmapp.put("unit1", unit1);
    MGBsgmapp.put("label2", label2);
    MGBsgmapp.put("unit2", unit2);
    MGBsgmapp<< fW;

    // sigma approximation WAWB^sB, GsB

    DblNumMat sgmG(nx, LEN);
    for (int l=0; l<LEN; l++) {
	for (int ix=0; ix<nx; ix++) {
	    sgmG._data[l*nx+ix] = G._data[l*nx+ix]*sigmp._data[l];
	}
    }
    dat = sgmG.data();
    for (int k=0; k<nx*LEN; k++) {
	fMlr[k] = dat[k];
    }
    MsgmG.put("n1", nx);
    MsgmG.put("n2", LEN);
    MsgmG << fMlr;

    DblNumMat sgmapp(nx, nx);
    iC( ddgemm(1.0, sgmG, B, 0.0, sgmapp) );
    ddat = sgmapp.data();
    if (norm) {
	for (int k=0; k<nx; k++) {
	    sign = (kss[k]>0)-(kss[k]<0);
	    for (int ix=0; ix<nx; ix++) {
		fW[k*nx+ix] = ddat[k*nx+ix]*vs[ix]*dt*0.5*sign;
	    } 
	}
    } else {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = ddat[k];
	}
    }

    float maxsgm, minsgm;
    maxsgm = fW.max();
    minsgm = fW.min();
    sf_warning("Sigma matrix: max value=%f; min value=%f",maxsgm, minsgm);

    /*if (rescale) {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = fW[k]*maxexact/maxsgm;
	}
	}*/
	
    Msgmapp.put("n1", nx);
    Msgmapp.put("n2", nx);
    Msgmapp.put("d1", dx);
    Msgmapp.put("d2", dk*twopi);
    Msgmapp.put("o1", 0);
    Msgmapp.put("o2", o2); // -dk*B._n/2.0*twopi
    Msgmapp.put("label1", label1);
    Msgmapp.put("unit1", unit1);
    Msgmapp.put("label2", label2);
    Msgmapp.put("unit2", unit2);
    Msgmapp<< fW;
    
    // weight and taper 
    int itmp;
    float w=0.0, k0=0.0, kk=0.0;
    float tp=0.0, tpk0=0.0, tpk=0.0;
    valarray<float> wfun(Bc._n*ridx.size());
    DblNumMat wtG(nx, LEN), wBc(LEN, COUNT), wM2c(1, COUNT);
    DblNumMat tmpcoef(1, LEN);

    setvalue(coef, 0.0);
    setvalue(tmpG, 0.0);
    for (int ixm=0; ixm<(int)ridx.size(); ixm++) {
	k0 = twopi*dfrq/vs[ridx[ixm]];
	tpk0 = 0.5*SF_PI/(vs[ridx[ixm]]*dt*0.5);
	//for (int ik=0; ik<nx; ik++) {
	for (int ik=0; ik<COUNT; ik++) {
	    kk = twopi*ktmpc._data[ik];
	    if (weight == true ) {
		w = wghtfun(kk, k0);
		/* weight */
	    } else {
		w = 1.0;
	    }
	    tpk = twopi*ks[ksc[ik]];

	    if (taper == true ) {
		tp = tpfun(tpa, tpb, tpk, tpk0);
	    } else {
		tp = 1.0;
	    }

	    wM2c(0,ik) = tp*w*M2c(ixm, ik);
	    wfun[ik+ixm*Bc._n] = w;
	    for (int il=0; il<Bc._m; il++) {
		itmp = ik*B._m+il;
		wBc._data[itmp] = w*Bc._data[itmp];
		/* Weight*B */
	    }
	}
	setvalue(IBc, 0.0);
	iC( ddpinv(wBc, 1e-16, IBc) );
	
	setvalue(tmpcoef, 0.0);
	iC( ddgemm(1.0,wM2c,IBc,0.0,tmpcoef) );
	
	for (int il=0; il<LEN; il++) {
	    coef(ixm, il) = tmpcoef(0, il);
	}
    }
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,wtG));

    //DblNumMat wsG(nx, LEN), Msgm(LEN, LEN);
    //setvalue(Msgm, 0.0);
    //for (int ix=1; ix<=len; ix++) {
    //Msgm._data[(ix-1)*len+(ix-1)] = sin(ix*SF_PI/len)/(ix*SF_PI/len); 
    //}
    //iC(ddgemm(1.0,wG,Msgm,0.0,wsG));

    dat = wtG.data();
    for (int k=0; k<nx*LEN; k++) {
	fMlr[k] = dat[k];
    }
    MwtG.put("n1", nx);
    MwtG.put("n2", LEN);
    MwtG << fMlr;


    DblNumMat wtapp(nx, nx);
    //iC( ddgemm(1.0, wsG, B, 0.0, wsgmapp) );
    iC( ddgemm(1.0, wtG, B, 0.0, wtapp) );
    ddat = wtapp.data();
    if (norm) {
	for (int k=0; k<nx; k++) {
	    sign = (kss[k]>0)-(kss[k]<0);
	    for (int ix=0; ix<nx; ix++) {
		fW[k*nx+ix] = ddat[k*nx+ix]*vs[ix]*dt*0.5*sign;
	    } 
	}
    } else {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = ddat[k];
	}
    }
        
    Mwtapp.put("n1", nx);
    Mwtapp.put("n2", nx);
    Mwtapp.put("d1", dx);
    Mwtapp.put("d2", dk*twopi);
    Mwtapp.put("o1", 0);
    Mwtapp.put("o2", o2); // -dk*B._n/2.0*twopi
    Mwtapp.put("label1", label1);
    Mwtapp.put("unit1", unit1);
    Mwtapp.put("label2", label2);
    Mwtapp.put("unit2", unit2);
    Mwtapp<< fW;

    // weight and sigma approximation
    DblNumMat wsG(nx, LEN);
    setvalue(wM2c, 0.0);
    setvalue(wBc,0.0);
    setvalue( coef, 0.0);
    for (int ixm=0; ixm<(int)ridx.size(); ixm++) {
	k0 = twopi*dfrq/vs[ridx[ixm]];
	for (int ik=0; ik<COUNT; ik++) {
	    kk = twopi*ktmpc._data[ik];
	    if (weight == true ) {
		w = wghtfun(kk, k0);
		/* weight */
	    } else {
		w = 1.0;
	    }

	    wM2c(0,ik) = w*M2c(ixm, ik);
	    wfun[ik+ixm*Bc._n] = w;
	    for (int il=0; il<Bc._m; il++) {
		itmp = ik*Bc._m+il;
		wBc._data[itmp] = w*Bc._data[itmp];
		/* Weight*B */
	    }
	}
	setvalue(IBc, 0.0);
	iC( ddpinv(wBc, 1e-16, IBc) );
	
	setvalue(tmpcoef, 0.0);
	iC( ddgemm(1.0,wM2c,IBc,0.0,tmpcoef) );
	
	for (int il=0; il<LEN; il++) {
	    coef(ixm, il) = tmpcoef(0, il);
	}
    }
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,wsG));

    for (int l=0; l<LEN; l++) {
	for (int ix=0; ix<nx; ix++) {
	    wsG._data[l*nx+ix] = wsG._data[l*nx+ix]*sigmp._data[l];
	}
    }
    
    dat = wsG.data();
    for (int k=0; k<nx*LEN; k++) {
	fMlr[k] = dat[k];
    }
    MwsG.put("n1", nx);
    MwsG.put("n2", LEN);
    MwsG << fMlr;

    DblNumMat wsgmapp(nx, nx);
    iC( ddgemm(1.0, wsG, B, 0.0, wsgmapp) );
    ddat = wsgmapp.data();
    if (norm) {
	for (int k=0; k<nx; k++) {
	    sign = (kss[k]>0)-(kss[k]<0);
	    for (int ix=0; ix<nx; ix++) {
		fW[k*nx+ix] = ddat[k*nx+ix]*vs[ix]*dt*0.5*sign;
	    }
	}
    } else {
	for (int k=0; k<nx*nx; k++) {
	    fW[k] = ddat[k];
	}
    }
    
    Mwsgmapp.put("n1", nx);
    Mwsgmapp.put("n2", nx);
    Mwsgmapp.put("d1", dx);
    Mwsgmapp.put("d2", dk*twopi);
    Mwsgmapp.put("o1", 0);
    Mwsgmapp.put("o2", o2); // -dk*B._n/2.0*twopi
    Mwsgmapp.put("label1", label1);
    Mwsgmapp.put("unit1", unit1);
    Mwsgmapp.put("label2", label2);
    Mwsgmapp.put("unit2", unit2);
    Mwsgmapp<< fW;

    /* Dispersion */
    std::valarray<float> fDisp(nx);
    double tg;
    if (disp == true ) {
	/* Dispersion */
	oRSF Dlfd("Dlfd"); 
	oRSF Dwt("Dwt");
	oRSF Dsgm("Dsgm");
    
	// LFD 
	fDisp[0] = 1.0; 
	for (int k=0; k < nx; k++) {
	    tg = 0.0;
	    for (int l=1; l <= len; l++) {
		tg += G(0,l-1)*sin((2*l-1)*dx*(ks[k])*SF_PI);
	    }
	    fDisp[k] = asin( dt*vs[0]*tg/2.0 )/((ks[k])*SF_PI*dt*vs[0])-1.0; 
	}
	Dlfd.put("n1", nx);
	Dlfd.put("n2", 1);
	Dlfd.put("d2", dk);
	Dlfd.put("o2",0);
 	Dlfd << fDisp;
	
	// LFD with weight and taper
	fDisp[0] = 1.0; 
	for (int k=0; k < nx; k++) {
	    tg = 0.0;
	    for (int l=1; l <= len; l++) {
		tg += wtG(0,l-1)*sin((2*l-1)*dx*(ks[k])*SF_PI);
	    }
	    fDisp[k] = asin( dt*vs[0]*tg/2.0 )/((ks[k])*SF_PI*dt*vs[0])-1.0; 
	}
	Dwt.put("n1", nx);
	Dwt.put("n2", 1);
	Dwt.put("d2", dk);
	Dwt.put("o2",0);
 	Dwt << fDisp;

	// LFD with weight and taper
	fDisp[0] = 1.0; 
	for (int k=0; k < nx; k++) {
	    tg = 0.0;
	    for (int l=1; l <= len; l++) {
		tg += sgmG(0,l-1)*sin((2*l-1)*dx*(ks[k])*SF_PI);
	    }
	    fDisp[k] = asin( dt*vs[0]*tg/2.0 )/((ks[k])*SF_PI*dt*vs[0])-1.0; 
	}
	Dsgm.put("n1", nx);
	Dsgm.put("n2", 1);
	Dsgm.put("d2", dk);
	Dsgm.put("o2",0);
 	Dsgm << fDisp;
    }
	
    

       


    // test
    /*DblNumMat tA(4,1);
    DblNumMat tB(1,5);
    DblNumMat tAB(4,5);
    for (int k=0; k<4; k++) { 
	tA._data[k]=1;
    }
    for (int k=0; k<5; k++) { 
	tB._data[k]=k;
    }

    iC( ddgemm(1.0, tA, tB, 0.0, tAB));

    
    for (int kk=0; kk<5; kk++){
	for (int ik=0; ik<4; ik++) {    
	    sf_warning("test[%d, %d]=%f",kk, ik, tAB._data[ik*5+kk]);
	}
    }
    for (int ik=0; ik<20; ik++) {    
	sf_warning("test[%d, %d]=%f",ik, ik, tAB._data[ik]);
	}*/
    
    

    
    return 0;
}

	

    
