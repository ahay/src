/* 2-D two-components elastic wavefield extrapolation and vector decomposition
   simultaneously using low-rank approximation on the base of 
   displacement wave equation and polarization projection in VTI media.

   Copyright (C) 2014 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
     
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

#include <rsf.hh>
#include <assert.h>

/* low rank decomposition  */
#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

/* prepared head files by myself */
#include "_cjb.h"

/* head files aumatically produced from C programs */
extern "C"{
#include "zero.h"
#include "ricker.h"
#include "kykxkztaper.h"
#include "fwpvtielowrank.double.h"
#include "eigen2x2.h"
#include "decomplowrank.h"
}

static std::valarray<float> vp, vs, ep, de, th;
static double dt1;

static std::valarray<double> sinx, cosx, rk2;

/* dual-domain operators based on low-rank decomp. */
int sampleopx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleopz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);

int sampleosx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx);
int sampleosz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx);

static void map2d1d(double *d, DblNumMat mat, int m, int n);
/*****************************************************************************************/
int main(int argc, char* argv[])
{
    sf_init(argc,argv);

    clock_t t1, t2, t3;
    float   timespent;

    t1=clock();

    iRSF par(0);
    int seed;
    par.get("seed",seed,time(NULL)); // seed for random number generator
    srand48(seed);

    float eps;
    par.get("eps",eps,1.e-8); // tolerance
       
    int npk;
    par.get("npk",npk,60); // maximum rank

    int   ns;
    float dt;
    par.get("ns",ns);
    par.get("dt",dt);
    dt1 = (double)dt;

    sf_warning("ns=%d dt=%f",ns,dt);
    sf_warning("npk=%d ",npk);
    sf_warning("eps=%16.12f",eps);
    sf_warning("read velocity model parameters"); 
    /* setup I files */
    iRSF vp0, vs0("vs0"), epsi("epsi"), del("del"), the("the");

    /* Read/Write axes */
    int nxv, nzv;
    vp0.get("n1",nzv);
    vp0.get("n2",nxv);

    float az, ax;
    vp0.get("o1",az);
    vp0.get("o2",ax);

    float fx, fz;
    fx=ax*1000.0;
    fz=az*1000.0;

    float dx, dz;
    vp0.get("d1",az);
    vp0.get("d2",ax);
    dz = az*1000.0;
    dx = ax*1000.0;

    sf_warning("dx=%f dz=%f",dx,dz);

    /* wave modeling space */
    int nx, nz, nxz;
    nx=nxv;
    nz=nzv;
    nxz=nx*nz;

    sf_warning("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    sf_warning("Warning: 2nd-order spectral need odd-based FFT");
    sf_warning("Warning: 2nd-order spectral need odd-based FFT");
    sf_warning("Warning: 2nd-order spectral need odd-based FFT");
    if(nx%2==0||nz%2==0) exit(0);
    sf_warning("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");

    vp.resize(nxz);
    vs.resize(nxz);
    ep.resize(nxz);
    de.resize(nxz);
    th.resize(nxz);
 
    vp0>>vp;
    vs0>>vs;
    epsi>>ep;
    del>>de;
    the>>th;

    for(int i=0;i<nxz;i++)
	th[i] *= SF_PI/180.0;

    /* Fourier spectra demension */
    int nkz,nkx,nk;
    nkx=nx;
    nkz=nz;
    nk = nkx*nkz;

    double dkz,dkx,kz0,kx0;

    dkx=2*SF_PI/dx/nx;
    dkz=2*SF_PI/dz/nz;

    kx0=-SF_PI/dx;
    kz0=-SF_PI/dz;

    sf_warning("dkx=%f dkz=%f",dkx,dkz);

    sinx.resize(nk);
    cosx.resize(nk);
    rk2.resize(nk);

    double *akx = (double*)malloc(sizeof(double)*nk);
    double *akz = (double*)malloc(sizeof(double)*nk);

    double kx, kz, rk, k2;
    int    i=0, j=0, k=0, ix, iz;
   
    for(ix=0; ix < nkx; ix++)
    {
	kx = kx0+ix*dkx;

	for (iz=0; iz < nkz; iz++)
	{
            kz = kz0+iz*dkz;

	    akx[i] = kx;
	    akz[i] = kz;

            k2 = kx*kx+kz*kz;
            rk = sqrt(k2);

            sinx[i] = kx/rk;
            cosx[i] = kz/rk;
            rk2[i] = k2;

            i++;
	}
    }

    /*****************************************************************************
     *  Calculating polarization deviation operator for wave-mode separation
     * ***************************************************************************/
    vector<int> md(nxz), nd(nk);
    for (k=0; k < nxz; k++)  md[k] = k;
    for (k=0; k < nk; k++)  nd[k] = k;

    vector<int> lid, rid;
    DblNumMat mid, mat;

    /********* qP-wave low rank decomposition of operator BpxAx + BpxzAxz applying to ux **********/
    int   m2opx1, n2opx1;
    double *ldataopx1, *fmidopx1, *rdataopx1;

    iC( ddlowrank(nxz,nk,sampleopx1,eps,npk,lid,rid,mid) );
    m2opx1=mid.m();
    n2opx1=mid.n();
    sf_warning("m2opx1=%d n2opx1=%d",m2opx1, n2opx1);

    fmidopx1  = (double*)malloc(sizeof(double)*m2opx1*n2opx1);
    ldataopx1 = (double*)malloc(sizeof(double)*nxz*m2opx1);
    rdataopx1 = (double*)malloc(sizeof(double)*n2opx1*nk);

    map2d1d(fmidopx1, mid, m2opx1, n2opx1);

    iC ( sampleopx1(md,lid,mat) );
    map2d1d(ldataopx1, mat, nxz, m2opx1);

    iC ( sampleopx1(rid,nd,mat) );
    map2d1d(rdataopx1, mat, n2opx1, nk);
    /*  FILE *fmid,*fleft,*fright;
	fmid = fopen("mid","wb");
	fleft = fopen("left","wb");
	fright = fopen("right","wb");
	float temp1;
	for(i=0;i<m2opx1*n2opx1;i++)
	{
	temp1 = fmidopx1[i];
	fwrite(&temp1,sizeof(float),1,fmid);
	}
	for(i=0;i<nxz*m2opx1;i++)
	{
	temp1 = ldataopx1[i];
	fwrite(&temp1,sizeof(float),1,fleft);
	}
	for(i=0;i<n2opx1*nk;i++)
	{
	temp1 = rdataopx1[i];
	fwrite(&temp1,sizeof(float),1,fright);
	}
	fclose(fmid);
	fclose(fleft);
	fclose(fright);*/

    sf_warning("rdata=%ef",rdataopx1[0]);
    sf_warning("ldata=%ef",ldataopx1[0]);
    sf_warning("fmid=%ef",fmidopx1[0]);

    /********* qP-wave low rank decomposition of operator BpxAxz + BpxzAz applying to uz **********/
    int   m2opx2, n2opx2;
    double *ldataopx2, *fmidopx2, *rdataopx2;

    iC( ddlowrank(nxz,nk,sampleopx2,eps,npk,lid,rid,mid) );
    m2opx2=mid.m();
    n2opx2=mid.n();
    sf_warning("m2opx2=%d n2opx2=%d",m2opx2, n2opx2);

    fmidopx2  = (double*)malloc(sizeof(double)*m2opx2*n2opx2);
    ldataopx2 = (double*)malloc(sizeof(double)*nxz*m2opx2);
    rdataopx2 = (double*)malloc(sizeof(double)*n2opx2*nk);

    map2d1d(fmidopx2, mid, m2opx2, n2opx2);

    iC ( sampleopx2(md,lid,mat) );
    map2d1d(ldataopx2, mat, nxz, m2opx2);

    iC ( sampleopx2(rid,nd,mat) );
    map2d1d(rdataopx2, mat, n2opx2, nk);

    /********* qP-wave low rank decomposition of operator BpxzAx + BpzAxz applying to ux **********/
    int   m2opz1, n2opz1;
    double *ldataopz1, *fmidopz1, *rdataopz1;

    iC( ddlowrank(nxz,nk,sampleopz1,eps,npk,lid,rid,mid) );
    m2opz1=mid.m();
    n2opz1=mid.n();
    sf_warning("m2opz1=%d n2opz1=%d",m2opz1, n2opz1);

    fmidopz1  = (double*)malloc(sizeof(double)*m2opz1*n2opz1);
    ldataopz1 = (double*)malloc(sizeof(double)*nxz*m2opz1);
    rdataopz1 = (double*)malloc(sizeof(double)*n2opz1*nk);

    map2d1d(fmidopz1, mid, m2opz1, n2opz1);

    iC ( sampleopz1(md,lid,mat) );
    map2d1d(ldataopz1, mat, nxz, m2opz1);

    iC ( sampleopz1(rid,nd,mat) );
    map2d1d(rdataopz1, mat, n2opz1, nk);

    /********* qP-wave low rank decomposition of operator BpxAxz + BpxzAz applying to uz **********/
    int   m2opz2, n2opz2;
    double *ldataopz2, *fmidopz2, *rdataopz2;

    iC( ddlowrank(nxz,nk,sampleopz2,eps,npk,lid,rid,mid) );
    m2opz2=mid.m();
    n2opz2=mid.n();
    sf_warning("m2opz2=%d n2opz2=%d",m2opz2, n2opz2);

    fmidopz2  = (double*)malloc(sizeof(double)*m2opz2*n2opz2);
    ldataopz2 = (double*)malloc(sizeof(double)*nxz*m2opz2);
    rdataopz2 = (double*)malloc(sizeof(double)*n2opz2*nk);

    map2d1d(fmidopz2, mid, m2opz2, n2opz2);

    iC ( sampleopz2(md,lid,mat) );
    map2d1d(ldataopz2, mat, nxz, m2opz2);

    iC ( sampleopz2(rid,nd,mat) );
    map2d1d(rdataopz2, mat, n2opz2, nk);

    /********* qSV-wave low rank decomposition of operator BsxAx + BsxzAxz applying to ux **********/
    int   m2osx1, n2osx1;
    double *ldataosx1, *fmidosx1, *rdataosx1;

    iC( ddlowrank(nxz,nk,sampleosx1,eps,npk,lid,rid,mid) );
    m2osx1=mid.m();
    n2osx1=mid.n();
    sf_warning("m2osx1=%d n2osx1=%d",m2osx1, n2osx1);

    fmidosx1  = (double*)malloc(sizeof(double)*m2osx1*n2osx1);
    ldataosx1 = (double*)malloc(sizeof(double)*nxz*m2osx1);
    rdataosx1 = (double*)malloc(sizeof(double)*n2osx1*nk);

    map2d1d(fmidosx1, mid, m2osx1, n2osx1);

    iC ( sampleosx1(md,lid,mat) );
    map2d1d(ldataosx1, mat, nxz, m2osx1);

    iC ( sampleosx1(rid,nd,mat) );
    map2d1d(rdataosx1, mat, n2osx1, nk);

    /********* qP-wave low rank decomposition of operator BpxAxz + BpxzAz applying to uz **********/
    int   m2osx2, n2osx2;
    double *ldataosx2, *fmidosx2, *rdataosx2;

    iC( ddlowrank(nxz,nk,sampleosx2,eps,npk,lid,rid,mid) );
    m2osx2=mid.m();
    n2osx2=mid.n();
    sf_warning("m2osx2=%d n2osx2=%d",m2osx2, n2osx2);

    fmidosx2  = (double*)malloc(sizeof(double)*m2osx2*n2osx2);
    ldataosx2 = (double*)malloc(sizeof(double)*nxz*m2osx2);
    rdataosx2 = (double*)malloc(sizeof(double)*n2osx2*nk);

    map2d1d(fmidosx2, mid, m2osx2, n2osx2);

    iC ( sampleosx2(md,lid,mat) );
    map2d1d(ldataosx2, mat, nxz, m2osx2);

    iC ( sampleosx2(rid,nd,mat) );
    map2d1d(rdataosx2, mat, n2osx2, nk);

    /********* qP-wave low rank decomposition of operator BpxzAx + BpzAxz applying to ux **********/
    int   m2osz1, n2osz1;
    double *ldataosz1, *fmidosz1, *rdataosz1;

    iC( ddlowrank(nxz,nk,sampleosz1,eps,npk,lid,rid,mid) );
    m2osz1=mid.m();
    n2osz1=mid.n();
    sf_warning("m2osz1=%d n2osz1=%d",m2osz1, n2osz1);

    fmidosz1  = (double*)malloc(sizeof(double)*m2osz1*n2osz1);
    ldataosz1 = (double*)malloc(sizeof(double)*nxz*m2osz1);
    rdataosz1 = (double*)malloc(sizeof(double)*n2osz1*nk);

    map2d1d(fmidosz1, mid, m2osz1, n2osz1);

    iC ( sampleosz1(md,lid,mat) );
    map2d1d(ldataosz1, mat, nxz, m2osz1);

    iC ( sampleosz1(rid,nd,mat) );
    map2d1d(rdataosz1, mat, n2osz1, nk);

    /********* qP-wave low rank decomposition of operator BpxAxz + BpxzAz applying to uz **********/
    int   m2osz2, n2osz2;
    double *ldataosz2, *fmidosz2, *rdataosz2;

    iC( ddlowrank(nxz,nk,sampleosz2,eps,npk,lid,rid,mid) );
    m2osz2=mid.m();
    n2osz2=mid.n();
    sf_warning("m2osz2=%d n2osz2=%d",m2osz2, n2osz2);

    fmidosz2  = (double*)malloc(sizeof(double)*m2osz2*n2osz2);
    ldataosz2 = (double*)malloc(sizeof(double)*nxz*m2osz2);
    rdataosz2 = (double*)malloc(sizeof(double)*n2osz2*nk);

    map2d1d(fmidosz2, mid, m2osz2, n2osz2);

    iC ( sampleosz2(md,lid,mat) );
    map2d1d(ldataosz2, mat, nxz, m2osz2);

    iC ( sampleosz2(rid,nd,mat) );
    map2d1d(rdataosz2, mat, n2osz2, nk);

    /****************End of Calculating Projection Deviation Operator****************/
    t2=clock();
    timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
    sf_warning("CPU time for low-rank decomp: %f(second)",timespent);

    /****************begin to calculate wavefield****************/
    /****************begin to calculate wavefield****************/
    /*  wavelet parameter for source definition */
    float A, f0, t0;
    f0=30.0;
    t0=0.04;
    A=1.0;

    sf_warning("fx=%f fz=%f dx=%f dz=%f",fx,fz,dx,dz);
    sf_warning("nx=%d nz=%d ", nx,nz);

    /* source definition */
    int ixs, izs;
    ixs=nxv/2;
    izs=nzv/2;

    /* setup I/O files */
    oRSF Elasticx("out");
    oRSF Elasticz("Elasticz");
    oRSF ElasticPx("ElasticPx");
    oRSF ElasticPz("ElasticPz");
    oRSF ElasticSx("ElasticSx");
    oRSF ElasticSz("ElasticSz");

    Elasticx.put("n1",nkz);
    Elasticx.put("n2",nkx);
    Elasticx.put("d1",dz/1000);
    Elasticx.put("d2",dx/1000);
    Elasticx.put("o1",fz/1000);
    Elasticx.put("o2",fx/1000);

    Elasticz.put("n1",nkz);
    Elasticz.put("n2",nkx);
    Elasticz.put("d1",dz/1000);
    Elasticz.put("d2",dx/1000);
    Elasticz.put("o1",fz/1000);
    Elasticz.put("o2",fx/1000);

    ElasticPx.put("n1",nkz);
    ElasticPx.put("n2",nkx);
    ElasticPx.put("d1",dz/1000);
    ElasticPx.put("d2",dx/1000);
    ElasticPx.put("o1",fz/1000);
    ElasticPx.put("o2",fx/1000);

    ElasticPz.put("n1",nkz);
    ElasticPz.put("n2",nkx);
    ElasticPz.put("d1",dz/1000);
    ElasticPz.put("d2",dx/1000);
    ElasticPz.put("o1",fz/1000);
    ElasticPz.put("o2",fx/1000);

    ElasticSx.put("n1",nkz);
    ElasticSx.put("n2",nkx);
    ElasticSx.put("d1",dz/1000);
    ElasticSx.put("d2",dx/1000);
    ElasticSx.put("o1",fz/1000);
    ElasticSx.put("o2",fx/1000);

    ElasticSz.put("n1",nkz);
    ElasticSz.put("n2",nkx);
    ElasticSz.put("d1",dz/1000);
    ElasticSz.put("d2",dx/1000);
    ElasticSz.put("o1",fz/1000);
    ElasticSz.put("o2",fx/1000);

    /********************* wavefield extrapolation *************************/
    double *px1=(double*)malloc(sizeof(double)*nxz);
    double *px2=(double*)malloc(sizeof(double)*nxz);
    double *px3=(double*)malloc(sizeof(double)*nxz);
    double *pz1=(double*)malloc(sizeof(double)*nxz);
    double *pz2=(double*)malloc(sizeof(double)*nxz);
    double *pz3=(double*)malloc(sizeof(double)*nxz);

    double *svx1=(double*)malloc(sizeof(double)*nxz);
    double *svx2=(double*)malloc(sizeof(double)*nxz);
    double *svx3=(double*)malloc(sizeof(double)*nxz);
    double *svz1=(double*)malloc(sizeof(double)*nxz);
    double *svz2=(double*)malloc(sizeof(double)*nxz);
    double *svz3=(double*)malloc(sizeof(double)*nxz);

    double *ux=(double*)malloc(sizeof(double)*nxz);
    double *uz=(double*)malloc(sizeof(double)*nxz);

    double *pp=(double*)malloc(sizeof(double)*nxz);
    double *ppp=(double*)malloc(sizeof(double)*nxz);

    zero1double(px1, nxz);
    zero1double(px2, nxz);
    zero1double(px3, nxz);
    zero1double(pz1, nxz);
    zero1double(pz2, nxz);
    zero1double(pz3, nxz);

    zero1double(svx1, nxz);
    zero1double(svx2, nxz);
    zero1double(svx3, nxz);
    zero1double(svz1, nxz);
    zero1double(svz2, nxz);
    zero1double(svz3, nxz);

    zero1double(ux, nxz);
    zero1double(uz, nxz);

    int *ijkx = sf_intalloc(nkx);
    int *ijkz = sf_intalloc(nkz);

    ikxikz(ijkx, ijkz, nkx, nkz);

    std::valarray<float> x(nxz);

    /* Setting Stability Conditions, by Chenlong Wang & Zedong Wu */
    float fmax = 3*f0;
    float kxm, kzm, kxzm;
    float amin, bmin;
    amin = 99999999999;
    bmin = 99999999999;
    float c11, c33, c44;
    float c1144, c3344;
    i=0;
    for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	{
	    c33 = vp[i] * vp[i];
	    c44 = vs[i] * vs[i];
	    c11 = (1+2*ep[i])*c33;
	    c1144 = c11 + c44;
	    c3344 = c33 + c44;

	    if (c1144<amin)
		amin = c1144;
	    if (c3344<bmin)
		bmin = c3344;
	    i++;
	}
    float kxmax = kx0 + nkx*dkx;
    float kzmax = kz0 + nkz*dkz;
    kxm = 2*sqrt(2)*SF_PI*fmax/sqrt(amin);
    kzm = 2*sqrt(2)*SF_PI*fmax/sqrt(bmin);
    float abmin = MIN(amin, bmin);
    kxzm = 2*sqrt(2)*SF_PI*fmax/sqrt(abmin);

    cerr<<"max kx="<<kxmax<<endl;
    cerr<<"max kz="<<kzmax<<endl;
    cerr<<"kxm="<<kxm<<endl;
    cerr<<"kzm="<<kzm<<endl;
    cerr<<"kxzm="<<kxzm<<endl;

    for(int it=0;it<ns;it++)
    {
        float t=it*dt;

        if(it%100==0)
	    sf_warning("Elastic: it= %d  t=%f(s)",it,t);
 
	// 2D exploding force source
	for(i=-1;i<=1;i++)
	    for(j=-1;j<=1;j++)
	    {
		if(fabs(i)+fabs(j)==2)
		{
		    ux[(ixs+i)*nz+(izs+j)]+=i*Ricker(t, f0, t0, A);
		    uz[(ixs+i)*nz+(izs+j)]+=j*Ricker(t, f0, t0, A);

		    if(it==2||it%100==0) sf_warning("ux=%f uz=%f ",ux[(ixs+i)*nz+(izs+j)],uz[(ixs+i)*nz+(izs+j)]);
		}
	    }
        /* extrapolation of Upx-componet */
	zero1double(pp, nxz);
	// extrapolation without stability consideration
        //fwpvti2delowrank(ldataopx1,rdataopx1,fmidopx1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2opx1,n2opx1);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataopx1,rdataopx1,fmidopx1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2opx1,n2opx1, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) ppp[i] = pp[i];
	zero1double(pp, nxz);
        //fwpvti2delowrank(ldataopx2,rdataopx2,fmidopx2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2opx2,n2opx2);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataopx2,rdataopx2,fmidopx2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2opx2,n2opx2, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) px3[i] = ppp[i] + pp[i] - px1[i];

        /* extrapolation of Upz-componet */
	zero1double(pp, nxz);
        //fwpvti2delowrank(ldataopz1,rdataopz1,fmidopz1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2opz1,n2opz1);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataopz1,rdataopz1,fmidopz1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2opz1,n2opz1, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) ppp[i] = pp[i];
	zero1double(pp, nxz);
        //fwpvti2delowrank(ldataopz2,rdataopz2,fmidopz2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2opz2,n2opz2);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataopz2,rdataopz2,fmidopz2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2opz2,n2opz2, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) pz3[i] = ppp[i] + pp[i] - pz1[i];

        /* extrapolation of Usvx-componet */
	zero1double(pp, nxz);
        //fwpvti2delowrank(ldataosx1,rdataosx1,fmidosx1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2osx1,n2osx1);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataosx1,rdataosx1,fmidosx1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2osx1,n2osx1, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) ppp[i] = pp[i];
	zero1double(pp, nxz);
        //fwpvti2delowrank(ldataosx2,rdataosx2,fmidosx2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2osx2,n2osx2);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataosx2,rdataosx2,fmidosx2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2osx2,n2osx2, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) svx3[i] = ppp[i] + pp[i] - svx1[i];

        /* extrapolation of Usvz-componet */
	zero1double(pp, nxz);
        //fwpvti2delowrank(ldataosz1,rdataosz1,fmidosz1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2osz1,n2osz1);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataosz1,rdataosz1,fmidosz1,pp,ux,ijkx,ijkz,nx,nz,nxz,nk,m2osz1,n2osz1, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) ppp[i] = pp[i];
	zero1double(pp, nxz);
        //fwpvti2delowrank(ldataosz2,rdataosz2,fmidosz2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2osz2,n2osz2);
	// extrapolation with stability consideration
        fwpvti2delowranksvd_double(ldataosz2,rdataosz2,fmidosz2,pp,uz,ijkx,ijkz,nx,nz,nxz,nk,m2osz2,n2osz2, kxm, kzm, kxzm, akx, akz);
        for(i=0;i<nxz;i++) svz3[i] = ppp[i] + pp[i] - svz1[i];

        /******* coupling the qP & qSV wavefield ********/
        for(i=0;i<nxz;i++){
	    px1[i]=px2[i];
	    px2[i]=px3[i];
	    pz1[i]=pz2[i];
	    pz2[i]=pz3[i];
	    svx1[i]=svx2[i];
	    svx2[i]=svx3[i];
	    svz1[i]=svz2[i];
	    svz2[i]=svz3[i];
	    ux[i] = px3[i]+svx3[i];
	    uz[i] = pz3[i]+svz3[i];
        }
	/*
	  if(it%100==0){
	  float e;
	  e = energy(ux,nxz)+energy(uz,nxz);
	  sf_warning("====================Total energy: %f",e);
	  }*/
        /******* output wavefields: components******/
        if(it==ns-1)
        {
	    for(i=0;i<nxz;i++) x[i]=ux[i];
	    Elasticx<<x;
	    for(i=0;i<nxz;i++) x[i]=uz[i];
	    Elasticz<<x;
	    for(i=0;i<nxz;i++) x[i]=px3[i];
	    ElasticPx<<x;
	    for(i=0;i<nxz;i++) x[i]=pz3[i];
	    ElasticPz<<x;
	    for(i=0;i<nxz;i++) x[i]=svx3[i];
	    ElasticSx<<x;
	    for(i=0;i<nxz;i++) x[i]=svz3[i];
	    ElasticSz<<x;
        }
        if(it%1==0&&it>=12000)
        {
            FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
            char dataxnumber[20],dataznumber[20];
            char datasxnumber[20],datapznumber[20];
            char datapxnumber[20],datasznumber[20];
            memset(dataxnumber,0,10);
            memset(dataznumber,0,10);
            memset(datapxnumber,0,10);
            memset(datasxnumber,0,10);
            memset(datapznumber,0,10);
            memset(datasznumber,0,10);

            sprintf(dataznumber,"dataz%d",it/1);
            sprintf(dataxnumber,"datax%d",it/1);
            sprintf(datapxnumber,"datapx%d",it/1);
            sprintf(datasxnumber,"datasx%d",it/1);
            sprintf(datapznumber,"datapz%d",it/1);
            sprintf(datasznumber,"datasz%d",it/1);

            fp1 = fopen(dataxnumber,"wb");
            fp2 = fopen(dataznumber,"wb");
            fp3 = fopen(datapxnumber,"wb");
            fp4 = fopen(datapznumber,"wb");
            fp5 = fopen(datasxnumber,"wb");
            fp6 = fopen(datasznumber,"wb");
            for(i=0;i<nxz;i++) x[i]=ux[i];
	    fwrite(&x[0],sizeof(float),nxz,fp1);
            for(i=0;i<nxz;i++) x[i]=uz[i];
            fwrite(&x[0],sizeof(float),nxz,fp2);
            for(i=0;i<nxz;i++) x[i]=px3[i];
            fwrite(&x[0],sizeof(float),nxz,fp3);
            for(i=0;i<nxz;i++) x[i]=pz3[i];
            fwrite(&x[0],sizeof(float),nxz,fp4);
            for(i=0;i<nxz;i++) x[i]=svx3[i];
            fwrite(&x[0],sizeof(float),nxz,fp5);
            for(i=0;i<nxz;i++) x[i]=svz3[i];
            fwrite(&x[0],sizeof(float),nxz,fp6);
            fclose(fp1);
	    fclose(fp2);
	    fclose(fp3);
	    fclose(fp4);
	    fclose(fp5);
	    fclose(fp6);
	}

    } //* it loop */

    t3=clock();
    timespent=(float)(t3-t2)/CLOCKS_PER_SEC;
    sf_warning("CPU time for wavefield extrapolation.: %f(second)",timespent);

    free(ldataopx1);
    free(fmidopx1);
    free(rdataopx1);

    free(ldataopx2);
    free(fmidopx2);
    free(rdataopx2);

    free(ldataopz1);
    free(fmidopz1);
    free(rdataopz1);

    free(ldataopz2);
    free(fmidopz2);
    free(rdataopz2);

    free(ldataosx1);
    free(fmidosx1);
    free(rdataosx1);

    free(ldataosx2);
    free(fmidosx2);
    free(rdataosx2);

    free(ldataosz1);
    free(fmidosz1);
    free(rdataosz1);

    free(ldataosz2);
    free(fmidosz2);
    free(rdataosz2);

    free(px1);
    free(px2);
    free(px3);
    free(pz1);
    free(pz2);
    free(pz3);

    free(svx1);
    free(svx2);
    free(svx3);
    free(svz1);
    free(svz2);
    free(svz3);

    free(pp);
    free(ppp);

    free(ux);
    free(uz);

    free(ijkx);
    free(ijkz);

    free(akx);
    free(akz);

    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qP-wave low rank decomposition of operator BpxAx + BpxzAxz applying to ux **********/
int sampleopx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, ax, axz, u1, u2, v1, v2, u1v1, u2v2,u1_2,u2_2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }
            u1_2= u1*u1;
            u2_2= u2*u2;

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
            if(v1*cx - v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;

            /* kxz = sx*cx*k2;
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;
            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            ax=2.0 - dt2*( u1_2*va[0]*sinclam1_2 + u2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //ax = 2.0*cos(dt1*sqrt(c11*sx2+c44*cx2));
            //ax = 2.0-dt2*(c11*sx2+c44*cx2);
            //axz = dt2*c13c44*kxz;

            resx(a,b) = ax*u1*u1-axz*u1*u2;
              
	}// b loop
    }// a loop

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qP-wave low rank decomposition of operator BpxAxz + BpxzAz applying to uz **********/
int sampleopx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, az, axz, u1, u2,v1, v2, u1v1, u2v2,v1_2, v2_2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
            if(v1*cx - v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;
            v1_2 = v1*v1;
            v2_2 = v2*v2;

            /* kxz = sx*cx*k2;
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;

            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            az = 2.0 - dt2*(v1_2*va[0]*sinclam1_2 + v2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //az = 2.0*cos(dt1*sqrt(c44*sx2+c33*cx2));
	    //az = 2.0-dt2*(c44*sx2+c33*cx2);
	    //axz = dt2*c13c44*kxz;

            resx(a,b) =-axz*u1*u1+az*u1*u2;
              
	}// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qP-wave low rank decomposition of operator BpxzAx + BpzAxz applying to ux **********/
int sampleopz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, ax, axz, u1, u2,v1,v2,u1_2,u2_2,u1v1,u2v2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }
            u1_2= u1*u1;
            u2_2= u2*u2;

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
            if(v1*cx - v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;

            /* kxz = sx*cx*k2;
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;

            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            ax = 2.0 - dt2*(u1_2*va[0]*sinclam1_2 + u2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //ax = 2.0*cos(dt1*sqrt(c11*sx2+c44*cx2));
	    //ax = 2.0-dt2*(c11*sx2+c44*cx2);
            //axz = dt2*c13c44*kxz;

            resx(a,b) = ax*u1*u2-axz*u2*u2;
              
	}// b loop
    }// a loop

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qP-wave low rank decomposition of operator BpxzAxz + BpzAz applying to uz **********/
int sampleopz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, az, axz, u1, u2,v1, v2, v1_2, v2_2,u1v1,u2v2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;


    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
            if(v1*cx - v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;
            v1_2 = v1*v1;
            v2_2 = v2*v2;

            /* kxz = sx*cx*k2;
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;

            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            az = 2.0 - dt2*(v1_2*va[0]*sinclam1_2 + v2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //az = 2.0*cos(dt1*sqrt(c44*sx2+c33*cx2));
            //az = 2.0-dt2*(c44*sx2+c33*cx2);
            //axz = dt2*c13c44*kxz;

            resx(a,b) =-axz*u1*u2+az*u2*u2;
              
	}// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qSV-wave low rank decomposition of operator BsxAx + BsxzAxz applying to ux **********/
int sampleosx1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, ax, axz, u1, u2,v1, v2, u1v1, u2v2,u1_2,u2_2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }
            u1_2= u1*u1;
            u2_2= u2*u2;

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
            if(v1*cx - v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;

            /* kxz = sx*cx*k2;
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;

            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            ax = 2.0 -dt2*(u1_2*va[0]*sinclam1_2 + u2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //ax = 2.0*cos(dt1*sqrt(c11*sx2+c44*cx2));
            //ax = 2.0-dt2*(c11*sx2+c44*cx2);
            //axz = dt2*c13c44*kxz;

            resx(a,b) = ax*v1*v1-axz*v1*v2; 
              
	}// b loop
    }// a loop

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qSV-wave low rank decomposition of operator BsxAxz + BsxzAz applying to uz **********/
int sampleosx2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, az, axz, u1, u2,v1, v2, u1v1, u2v2,v1_2, v2_2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
            if(v1*cx - v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;
            v1_2 = v1*v1;
            v2_2 = v2*v2;

            /* kxz = sx*cx*k2; 
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;

            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            az = 2.0 - dt2*(v1_2*va[0]*sinclam1_2 + v2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //az = 2.0*cos(dt1*sqrt(c44*sx2+c33*cx2));
            //az = 2.0-dt2*(c44*sx2+c33*cx2);
            //axz = dt2*c13c44*kxz;

            resx(a,b) =-axz*v1*v1+az*v1*v2;
              
	}// b loop
    }// a loop

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qSV-wave low rank decomposition of operator BsxzAx + BszAxz applying to ux **********/
int sampleosz1(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, ax, axz, u1, u2,v1, v2, u1v1, u2v2,u1_2,u2_2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }
            u1_2= u1*u1;
            u2_2= u2*u2;

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
	    if(v1*cx -  v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;

            /* kxz = sx*cx*k2;
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;

            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            ax=2.0 - dt2*( u1_2*va[0]*sinclam1_2 + u2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //ax = 2.0*cos(dt1*sqrt(c11*sx2+c44*cx2));
            //ax = 2.0-dt2*(c11*sx2+c44*cx2);
            //axz = dt2*c13c44*kxz;

            resx(a,b) = ax*v1*v2-axz*v2*v2;
              
	}// b loop
    }// a loop

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/********* qSV-wave low rank decomposition of operator BsxzAxz + BszAz applying to uz **********/
int sampleosz2(vector<int>& rs, vector<int>& cs, DblNumMat& resx)
{
    int nr = rs.size();
    int nc = cs.size();

    resx.resize(nr,nc);

    setvalue(resx,0.0);

    double  aa[2][2],ve[2][2],va[2];  /*matrix, eigeinvector and eigeinvalues*/
    double  c44, c11, c33, c13c44, a11, a12, a22;
    double  sx, cx, az, axz, u1, u2,v1, v2, u1v1, u2v2,v1_2, v2_2;
    double  lam1, lam2, sinclam1, sinclam2, sinclam1_2, sinclam2_2;

    for(int a=0; a<nr; a++) 
    {
        int i=rs[a];
        double vp2 = vp[i]*vp[i];
        double vs2 = vs[i]*vs[i];
        double ep2 = 1.0+2*ep[i];
        double de2 = 1.0+2*de[i];
	double coss=cos(th[i]);
        double sins=sin(th[i]);

        for(int b=0; b<nc; b++)
        {
            double s = sinx[cs[b]];
            double c = cosx[cs[b]];
            double k2 = rk2[cs[b]];
            if(s==0&&c==0)
            {
		resx(a,b) = 0.0;
		continue;
            }

            c33=vp2;
            c44=vs2;
            c11=ep2*c33;
	    c13c44=sqrt((de2*c33-c44)*(c33-c44));

	    // rotatiing according to tilted symmetry axis
	    sx=s*coss+c*sins;
	    cx=c*coss-s*sins;

            // vector decomposition operators based on polarization
            a11= c11*sx*sx+c44*cx*cx;
            a12= c13c44*sx*cx;
            a22= c44*sx*sx+c33*cx*cx;

            aa[0][0] = a11;
            aa[0][1] = a12;
            aa[1][0] = a12;
            aa[1][1] = a22;

            dsolveSymmetric22(aa, ve, va);

            /* qP-wave's polarization vector */
            u1=ve[0][0];
            u2=ve[0][1];
            if(u1*sx + u2*cx <0) {
		u1 = -u1;
		u2 = -u2;
            }

            /* qS-wave's polarization vector */
            v1=ve[1][0];
            v2=ve[1][1];
            if(v1*cx -  v2*sx <0) {
		v1 = -v1;
		v2 = -v2;
            }
            u1v1= u1*v1;
            u2v2= u2*v2;
            v1_2 = v1*v1;
            v2_2 = v2*v2;

            /* kxz = sx*cx*k2;
	       sx2=sx*sx*k2;
	       cx2=cx*cx*k2; */

	    va[0] = va[0]*k2;
	    va[1] = va[1]*k2;
            // wavefield extrapolator
            double dt2=dt1*dt1;

            lam1=sqrt(va[0])*0.5*dt1;
            lam2=sqrt(va[1])*0.5*dt1;
            sinclam1=sin(lam1)/lam1;
            sinclam2=sin(lam2)/lam2;
            sinclam1_2=sinclam1*sinclam1;
            sinclam2_2=sinclam2*sinclam2;
            az = 2.0 - dt2*(v1_2*va[0]*sinclam1_2 + v2_2*va[1]*sinclam2_2);
            axz = dt2*(u1v1*va[0]*sinclam1_2 + u2v2*va[1]*sinclam2_2);

            //az = 2.0*cos(dt1*sqrt(c44*sx2+c33*cx2));
            //az = 2.0-dt2*(c44*sx2+c33*cx2);
            //axz = dt2*c13c44*kxz;

            resx(a,b) =-axz*v1*v2+az*v2*v2;
              
	}// b loop
    }// a loop

    return 0;
}

static void map2d1d(double *d, DblNumMat mat, int m, int n)
{
    int i, j, k;
    k=0;
    for (i=0; i < m; i++)
	for (j=0; j < n; j++)
	{
	    d[k] = (double)mat(i,j);
	    k++;
	}

}
