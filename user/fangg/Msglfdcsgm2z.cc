// 2D Lowrank FD coefficient of d/dx on staggered grid (optimized)

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

static std::valarray<float> vs;
static std::valarray<double> ks;
static std::valarray<double> kx, kz;
static float dt;
static float dfrq, a0;

static float sinc(float x)
{
    if (fabs(x)<=SF_EPS) return 1.0;
    return sinf(x)/(x+SF_EPS);
}

int samplez(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    res(a,b) = 2.0*SF_PI*kz[cs[b]]*sinc(SF_PI*vs[rs[a]]*ks[cs[b]]*dt);
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

float tpfun(float a, float b, float kx, float kz, float k0)
{  
    float res   = 0.0;
    float halfa = a/2.0;
    float ks    = sqrt(kx*kx+kz*kz);
    float kh, kh2;
    
    if ( fabs(ks-k0) < halfa ) {
	kh  = (ks - k0)/halfa;
	kh2 = kh*kh; 
	res = 1.0 - b*exp(1.0-1.0/(1.0-kh2));
    }
    else if ( fabs(ks+k0) < halfa ) {
	kh  = (ks + k0)/halfa;
	kh2 = kh*kh; 
	res = 1.0 - b*exp(1.0-1.0/(1.0-kh2));
    }
    else {
	res = 1.0;
    }
    return res;
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

    int size;
    par.get("size",size,6); // stencil length 
    iRSF velf;
    oRSF outm;  // FD coefficient of d/dx
    oRSF fsx("sx");
    oRSF fsz("sz"); 
    float dx, dz;

    float wavnumcut;
    par.get("wavnumcut",wavnumcut,1.0); // wavenumber cut percentile

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

    bool taper;
    float tpa, tpb;
    par.get("taper",taper, true);
    par.get("tpa", tpa, 0.5); //taper for stability
    par.get("tpb", tpb, 0.8);

    sf_warning("=======================");
    sf_warning(" 2D Lowrank FD coefficients with weighted least square and sigma filter");
    sf_warning("nx=%i, nz=%i, dx=%f dz=%f", nx, nz, dx, dz);
    sf_warning("size=%i", size);
    sf_warning("weight=%d",weight);
    sf_warning("taper=%d", taper);
    sf_warning("taper_a=%f, taper_b=%f", tpa, tpb);
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

    
    int COUNT= 0;
    kx.resize(nxz);
    kz.resize(nxz);
    float kx1, kz1;
    float kx0 = -dkx*nx/2.0; 
    float kz0 = -dkz*nz/2.0;
    
    float a = nx/3.0*dkx*wavnumcut; 
    float b = nz/3.0*dkz*wavnumcut;

    int i=0;
    float dkxz=dkx+dkz;
    int SMK=0;
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
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

 
    /*Low rank decomposition*/   
    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m,n,samplez,eps,npk,cidx,ridx,mid) );

    sf_warning("cidx.size=%d", cidx.size());
    sf_warning("ridx.size=%d", ridx.size());

    DblNumMat M1(m,cidx.size());
    vector<int> rs(m);
    for(int k=0; k<m; k++) rs[k]=k;
    iC( samplez(rs,cidx,M1) );
    DblNumMat M2(ridx.size(),n);
    vector<int> cs(n);
    for(int k=0; k<n; k++) cs[k]=k;
    iC( samplez(ridx,cs,M2) );

    /*FD coefficient*/

    /* d/dz */
    int len=0;
    len = (size+1)/2; //x: -len ,... -1, 1, ... len;
    std::valarray<float> ztmp(2*len);
    std::valarray<float> xtmp(len);
    
    for (int ix = 0; ix<len; ix++) {
	ztmp[ix] = -(len-ix)+0.5;  // - (2l-1)/2
	ztmp[ix+len] = ix+0.5;     //   (2l-1)/2
	xtmp[ix] = ix;
    }
    cerr<<"len  = " <<len <<endl;
    cerr<<"xtmp = "; for (int ix=0; ix<len; ix++) cerr<<xtmp[ix]<<", "; cerr<<endl;
    cerr<<"ztmp = "; for (int ix=0; ix<2*len; ix++) cerr<<ztmp[ix]<<", ";   cerr<<endl;

    int gdc = 0;
    for (int ix=0; ix<len; ix++) {
	for (int iz=0; iz<2*len; iz++) { 
	    if ((xtmp[ix]>0 && ztmp[iz]<1.5 && ztmp[iz]>-1.5) || (ztmp[iz]>0 && xtmp[ix] == 0)) {
		if ( (xtmp[ix]*dx*xtmp[ix]*dx+ztmp[iz]*dz*ztmp[iz]*dz) < dz*dz*(ztmp[0]*ztmp[0]+0.0001)) {gdc++;
		}
	    }
	    if ( (ztmp[iz]==1.5||ztmp[iz]==-1.5) && xtmp[ix]==1) {gdc++;}
	}
    }

    nk =0;
    DblNumMat sx(gdc,1), sz(gdc,1);
    for (int ix=0; ix<len; ix++) {
	for (int iz=0; iz<2*len; iz++) { 
	    if ((xtmp[ix]>0 && ztmp[iz]<1.5 && ztmp[iz]>-1.5) || (ztmp[iz]>0 && xtmp[ix] == 0)) {
		if ( (xtmp[ix]*dx*xtmp[ix]*dx+ztmp[iz]*dz*ztmp[iz]*dz) < dz*dz*(ztmp[0]*ztmp[0]+0.0001)) {
		    sz._data[nk] = ztmp[iz];
		    sx._data[nk] = xtmp[ix];
		    nk++;
		} 
	    }
	    if ( (ztmp[iz]==1.5||ztmp[iz]==-1.5) && xtmp[ix]==1) {
		sz._data[nk] = ztmp[iz];
		sx._data[nk] = xtmp[ix];
		nk++;
	    } 	
	}
    }
      
    cerr<<"[x,z]="; for(int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<","<<sz._data[k]<<"] ";    
    cerr<<endl;
    cerr<<"gdc "<<gdc<<" "<<endl;
    cerr<<"nk "<<nk<<" "<<endl;
    
    /* G */
    DblNumMat kxtmp(1,nxz); for(int k=0; k<nxz; k++) kxtmp._data[k]=kx[k];
    DblNumMat kztmp(1,nxz); for(int k=0; k<nxz; k++) kztmp._data[k]=kz[k];
    DblNumMat kxtmpc(1,COUNT); for(int k=0; k<COUNT; k++) kxtmpc._data[k]=kx[ksc[k]];
    DblNumMat kztmpc(1,COUNT); for(int k=0; k<COUNT; k++) kztmpc._data[k]=kz[ksc[k]];
    int LEN = sx._m;
    DblNumMat Bx(LEN,nxz), B(LEN,nxz);
    DblNumMat Bxc(LEN,COUNT), Bc(LEN,COUNT);
    iC(ddgemm(dz,sz,kztmp,0.0,B));
    iC(ddgemm(dx,sx,kxtmp,0.0,Bx));
    iC(ddgemm(dz,sz,kztmpc,0.0,Bc));
    iC(ddgemm(dx,sx,kxtmpc,0.0,Bxc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=sin(2.0*SF_PI*(B._data[k]+Bx._data[k]));
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=sin(2.0*SF_PI*(Bc._data[k]+Bxc._data[k]));
    DblNumMat IBc(COUNT,LEN);    iC( ddpinv(Bc, 1e-16, IBc) );
      
    DblNumMat coef(ridx.size(),LEN);
    DblNumMat M2c;
    iC( samplez(ridx,ksc,M2c) );
    //iC(ddgemm(1.0,M2c,IBc,0.0,coef));
 
    DblNumMat G(nxz,LEN), tmpG(mid._m,LEN);
    //iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    //iC(ddgemm(1.0,M1,tmpG,0.0,G));
    
    /* weight */
    int itmp;
    float w=0.0, k0=0.0, kkx=0.0, kkz=0.0;
    float tp=0.0, tpk0=0.0;
    valarray<float> wf(Bc._n*ridx.size());
    DblNumMat wBc(LEN, COUNT), wM2c(1, COUNT);
    DblNumMat tmpcoef(1, LEN); 
    
    setvalue(coef, 0.0);
    setvalue(tmpG, 0.0);
    if ( weight==true || taper == true ) {
	for (int ixm=0; ixm<(int)ridx.size(); ixm++) {
	    k0   = 2.0*SF_PI*dfrq/vs[ridx[ixm]];
	    tpk0 = SF_PI/(vs[ridx[ixm]]*dt);
	    for (int ik=0; ik<COUNT; ik++) {
		kkx = 2*SF_PI*kxtmpc._data[ik];
		kkz = 2*SF_PI*kztmpc._data[ik];
		if (weight == true) {
		    w = wfun(kkx, kkz, k0, a0);
		}else {
		    w = 1.0;
		}
		if (taper == true ) {
		    tp = tpfun(tpa, tpb, kkx, kkz, tpk0);
		}else {
		    tp = 1.0;
		}
		wM2c(0,ik) = tp*w*M2c(ixm, ik);
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
    
    switch (filter) {
	case 1 : // sigma filter
	    for (int k=0; k<LEN; k++) {
		sigm._data[k] = fabs(sinc(fabs(sx._data[k]+1)*SF_PI/len)*sinc((fabs(sz._data[k])+0.5)*SF_PI/len));
		//sf_warning("k=%i input=%g sincx=%g sincz=%g sinc=%g sigmk=%g",k, sx._data[k]*SF_PI/sizep, sinc(sx._data[k]*SF_PI/sizep), sinc(sz._data[k]*SF_PI/sizep), sinc(sx._data[k]*SF_PI/sizep)*sinc(sz._data[k]*SF_PI/sizep), sigm._data[k]);
	    }
	    break;
	case 2 : // cose filter
	    for (int k=0; k<LEN; k++) {
		sigm._data[k] = 0.5*(1+cos(fabs(sx._data[k]+1)*SF_PI/len))*0.5*(1+cos((fabs(sz._data[k])+0.5)*SF_PI/len));
	    }
	    break;
	case 3 : // exponential filter
	    double alpha;
	    alpha = -1.0*log(1.0e-15);
	    sf_warning("ALPHA=%g", alpha);
	    for (int k=0; k<LEN; k++) {
		sigm._data[k] = exp(-1.0*alpha*pow(fabs(sx._data[k]+1)/len,2.0))*exp(-1.0*alpha*pow((fabs(sz._data[k])+0.5)/len,4.0));
		//sigm._data[k] = exp(-1.0*alpha*pow((fabs(sx._data[k])-0.5)/len,2.0))*exp(-1.0*alpha*pow(fabs(sz._data[k])/len,2.0));
	    }
	    break;
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
    cerr<<"[x,z]="; for(int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<","<<sz._data[k]<<"] ";
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

    /* -------- */
    // max and min
    sf_warning("Begin caculate max and min value...");
    float maxval=0.0, minval=0.0, tmpsum=0.0;
    for (int ix=0; ix<nxz; ix+=100){
	for (int ik=0; ik<nxz; ik+=100) {
	    tmpsum=0.0;
	    for (int il=0; il<LEN; il++) tmpsum+=G(ix, il)*B(il, ik);
	    tmpsum = tmpsum*ks[ik]*vs[ix]*dt/(kz[ik]*2);
	    if (maxval<tmpsum) maxval=tmpsum;
	    if (minval>tmpsum) minval=tmpsum;
	}
    }
    sf_warning("Max value=%f; Min value=%f",maxval, minval);

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
           fs[k] = ldat[k];
    } 
    fsx << fs;
    ldat = sz.data();
    for (int k=0; k < LEN; k++) {
        fs[k] =  ldat[k]+0.5;
    } 
    fsz << fs;
    
    //  */
    return 0;
}

    
    
    
    
    
    
    
    
    






