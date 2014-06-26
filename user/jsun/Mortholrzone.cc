// Lowrank decomposition for 3-D orthorhombic wave propagation with zone's approximation. 
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

static std::valarray<float>  C11,C12,C13,C22,C23,C33,C44,C55,C66,q1,q2;
static std::valarray<float> kx, ky, kz;
static float dt;
static int mode;
static bool approx,reduce;
static float slope=0.8530492393, intcp=0.145965755;

int sample(vector<int>& rs, vector<int>& cs, FltNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0f);
    //    float con2 = pow(2.0,1/3.0);
    for(int a=0; a<nr; a++) {
        int i=rs[a];
	double c11 = C11[i];
	double c12 = C12[i];
	double c13 = C13[i];
	double c22 = C22[i];
	double c23 = C23[i];
	double c33 = C33[i];
	double c44 = C44[i];
	double c55 = C55[i];
	double c66 = C66[i];
	double ss1 = sin(q1[i]);
        double ss2 = sin(q2[i]);
        double cc1 = cos(q1[i]);
        double cc2 = cos(q2[i]);

	for(int b=0; b<nc; b++) {
           int j=cs[b];

	   double x0 = kx[j];
	   double y0 = ky[j];
	   double z0 = kz[j];
	    // rotation of coordinates
	   double x = x0*cc2+y0*ss2;
           double y =-x0*ss2*cc1+y0*cc2*cc1+z0*ss1;
	   double z = x0*ss2*ss1-y0*cc2*ss1+z0*cc1;
           double x2 = x*x;
           double y2 = y*y;
           double z2 = z*z;
	   double xy=x*y;
	   double xz=x*z;
	   double yz=y*z;
	   double r;
	   if (!approx) {
	       double H11= c11*x2+c66*y2+c55*z2;
	       double H22= c66*x2+c22*y2+c44*z2;
	       double H33= c55*x2+c44*y2+c33*z2;
	       double H12= (c12+c66)*xy;
	       double H13= (c13+c55)*xz;
	       double H23= (c23+c44)*yz;
	       double aa = -(H11+H22+H33);
	       double bb = H11*H22+H11*H33+H22*H33-H12*H12-H13*H13-H23*H23;
	       double cc = H11*H23*H23+H22*H13*H13+H33*H12*H12-H11*H22*H33-2*H12*H23*H13;
	       double dd = bb-aa*aa/3.0;
	       double qq = 2.0*aa*aa*aa/27.0-aa*bb/3.0+cc;
	       double Q  = pow(dd/3,3)+pow(qq/2,2);
	       if (Q>0) sf_warning ("!!Q is positive!! Q=%g \n", Q);
	       double cv,vv;
	       if (abs(dd)<0.0000001) {
		   r=0;
	       } else {
		   cv = -qq/(2*sqrt(abs(-dd*dd*dd/27.0)));
		   vv = acos(cv);
		   if (mode==1)
		       r = 2*sqrt(abs(-dd/3.0))*cos(vv/3.0+2.0*SF_PI/3.0)-aa/3.0;
		   else if (mode==2) 
		       r = 2*sqrt(abs(-dd/3.0))*cos(vv/3.0+2.0*2.0*SF_PI/3.0)-aa/3.0;
		   else
		       r = 2*sqrt(abs(-dd/3.0))*cos(vv/3.0)-aa/3.0;
		   
		   r = sqrt(abs(r));
	       }
	   } else {
	       if (mode!=0) sf_error("Only apporximation for qP mode!");
	       double e  = c11*x2 + c22*y2 + c33*z2;
	       double k  = hypot(hypot(x,y),z);
//	       if (k!=hypot(hypot(x0,y0),z0)) sf_warning("wavenumbers not equal!!!");
	       if (k!=0) {
		   double n1= x/k; double n2= y/k; double n3= z/k;
		   n1*=n1; n2*=n2; n3*=n3;
		   if (n1!=1 && n2!=1 && n3!=1) {
		       if (c22*(c33-c44)==0 || c33*(c22-c44)==0 || c11*(c33-c55)==0 || c33*(c11-c55)==0 || c22*(c11-c66)==0 || c11*(c22-c66)==0 )
			   sf_warning("Dividing zero when calculating qv&qh!");
		       double qv1= (pow(c23+c44,2)+c44*(c33-c44))/(c22*(c33-c44));
		       double qh1= reduce ? slope*qv1+intcp : (pow(c23+c44,2)+c44*(c22-c44))/(c33*(c22-c44));
		       double qv2= (pow(c13+c55,2)+c55*(c33-c55))/(c11*(c33-c55));
		       double qh2= reduce ? slope*qv2+intcp : (pow(c13+c55,2)+c55*(c11-c55))/(c33*(c11-c55));
		       double qv3= (pow(c12+c66,2)+c66*(c11-c66))/(c22*(c11-c66));
		       double qh3= reduce ? qv3 : (pow(c12+c66,2)+c66*(c22-c66))/(c11*(c22-c66));
		       double q1 = (qv1*n3+qh1*n2)/(1.-n1);
		       double q2 = (qv2*n3+qh2*n1)/(1.-n2);
		       double q3 = (qv3*n1+qh3*n2)/(1.-n3);
		       double qm = (q1-1.)*c22*c33*y2*z2 + (q2-1.)*c11*c33*x2*z2 + (q3-1.)*c11*c22*x2*y2;
		       if ( 2*(c33*(1.-qh1)+c22*(qv1-1.))*(c33*(qv1-qh1)+c22*(qh1+qv1-pow(qv1,2)-1.)) ==0 ||
			    2*(c33*(1.-qh1)+c22*(qv1-1.))*(c22*(qh1-qv1)+c33*(qh1+qv1-pow(qh1,2)-1.)) ==0 ||
			    2*(c33*(1.-qh2)+c11*(qv2-1.))*(c33*(qv2-qh2)+c11*(qh2+qv2-pow(qv2,2)-1.)) ==0 ||
			    2*(c33*(1.-qh2)+c11*(qv2-1.))*(c11*(qh2-qv2)+c33*(qh2+qv2-pow(qh2,2)-1.)) ==0 ||
			    2*(c11*(1.-qh3)+c22*(qv3-1.))*(c11*(qv3-qh3)+c22*(qh3+qv3-pow(qv3,2)-1.)) ==0 ||
			    2*(c11*(1.-qh3)+c22*(qv3-1.))*(c22*(qh3-qv3)+c11*(qh3+qv3-pow(qh3,2)-1.)) ==0 )
			   sf_warning("Dividing zero when calculating sv&sh!");
		       double sv1= c22*(c33-c22)*(qh1-1.)*pow(qv1-1.,2)/(2*(c33*(1.-qh1)+c22*(qv1-1.))*(c33*(qv1-qh1)+c22*(qh1+qv1-pow(qv1,2)-1.)));
		       double sh1= c33*(c33-c22)*(qv1-1.)*pow(qh1-1.,2)/(2*(c33*(1.-qh1)+c22*(qv1-1.))*(c22*(qh1-qv1)+c33*(qh1+qv1-pow(qh1,2)-1.)));
		       double sv2= c11*(c33-c11)*(qh2-1.)*pow(qv2-1.,2)/(2*(c33*(1.-qh2)+c11*(qv2-1.))*(c33*(qv2-qh2)+c11*(qh2+qv2-pow(qv2,2)-1.)));
		       double sh2= c33*(c33-c11)*(qv2-1.)*pow(qh2-1.,2)/(2*(c33*(1.-qh2)+c11*(qv2-1.))*(c11*(qh2-qv2)+c33*(qh2+qv2-pow(qh2,2)-1.)));
		       double sv3= c22*(c11-c22)*(qh3-1.)*pow(qv3-1.,2)/(2*(c11*(1.-qh3)+c22*(qv3-1.))*(c11*(qv3-qh3)+c22*(qh3+qv3-pow(qv3,2)-1.)));
		       double sh3= c11*(c11-c22)*(qv3-1.)*pow(qh3-1.,2)/(2*(c11*(1.-qh3)+c22*(qv3-1.))*(c22*(qh3-qv3)+c11*(qh3+qv3-pow(qh3,2)-1.)));
		       double s1 = (sh2*n3+sv3*n2)/(1.-n1);
		       double s2 = (sh1*n3+sh3*n1)/(1.-n2);
		       double s3 = (sv1*n1+sv2*n2)/(1.-n3);
		       double sm = s1*n1 + s2*n2 + s3*n3;
		       if (sm==0) sf_warning("Dividing zero: sm=%f !",sm);
		       r = sqrt(e*(1.-sm) + sm*sqrt(e*e + 2.*qm/sm));
		       if (r!=r) sf_warning("r=,n1=,n2=,n3=",r,n1,n2,n3);
		   } else r = sqrt(e);
	       } else r = 0.;
	   }
	   res(a,b) = 2*(cos(r*dt)-1); 
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
    par.get("eps",eps,1.e-4); // tolerance

    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step

    par.get("mode",mode,0);  // '0' means quasi-P (default), '1' means quasi-S, '2' means quasi-S2
    if (mode==0) sf_warning(">>>>> Using quasi-P mode! <<<<<");
    else if (mode==1) sf_warning(">>>>> Using quasi-S mode! <<<<<");
    else if (mode==2) sf_warning(">>>>> Using quasi-S2 mode! <<<<<");
    else sf_warning(">>>>> Invalid mode parameter, using default (P)! <<<<<");

    bool tilt;
    par.get("tilt",tilt,false);
    par.get("approx",approx,false); // if true, use zone's approximation instead of exact phase velocity
    par.get("reduce",reduce,false); // if true, use the expirical linear relationship between qv and qh

    if (approx) sf_warning("Using zone's approximation!");
    if (reduce) sf_warning("Using slope and intercept approximation!");    

    iRSF c11, c12("c12"), c13("c13"), c22("c22"), c23("c23"), c33("c33"), c44("c44"), c55("c55"), c66("c66");
    
    int nz,nx,ny;
    c11.get("n1",nz);
    c11.get("n2",nx);
    c11.get("n3",ny);
    int m = nx*ny*nz;

    C11.resize(m);
    C12.resize(m);
    C13.resize(m);
    C22.resize(m);
    C23.resize(m);
    C33.resize(m);
    C44.resize(m);
    C55.resize(m);
    C66.resize(m);
    q1.resize(m);
    q2.resize(m);
    c11 >> C11;
    c12 >> C12;
    c13 >> C13;
    c22 >> C22;
    c23 >> C23;
    c33 >> C33;
    c44 >> C44;
    c55 >> C55;
    c66 >> C66;

    if (tilt) {
	iRSF seta1("seta1");
	seta1 >> q1;
	seta1 >> q2;
	/* from degrees to radians */
	for (int im=0; im < m; im++) {
	    q1[im] *= SF_PI/180.;
	    q2[im] *= SF_PI/180.;
	}
    } else {
	sf_warning(">>>>> No tilting! <<<<<");
	for (int im=0; im < m; im++) {
	    q1[im] = 0.;
	    q2[im] = 0.;
	}
    }
    
    iRSF fft("fft");

    int nkz,nkx,nky;
    fft.get("n1",nkz);
    fft.get("n2",nkx);
    fft.get("n3",nky);

    float dkz,dkx,dky;
    fft.get("d1",dkz);
    fft.get("d2",dkx);
    fft.get("d3",dky);
    
    float kz0,kx0,ky0;
    fft.get("o1",kz0);
    fft.get("o2",kx0);
    fft.get("o3",ky0);


    int n = nkx*nky*nkz;
    kx.resize(n);
    ky.resize(n);
    kz.resize(n);
    int i = 0;
    for (int iy=0; iy < nky; iy++) {
        for (int ix=0; ix < nkx; ix++) {
            for (int iz=0; iz < nkz; iz++) {
                kx[i] = 2*SF_PI*(kx0+ix*dkx);
                ky[i] = 2*SF_PI*(ky0+iy*dky);
                kz[i] = 2*SF_PI*(kz0+iz*dkz);
                i++;
            }
        }
    }


    vector<int> lidx, ridx;
    FltNumMat mid;

    iC( lowrank(m,n,sample,(float) eps,npk,lidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    FltNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );

    FltNumMat lmat2(m,n2);
    iC( dgemm(1.0, lmat, mid, 0.0, lmat2) );

    float *ldat = lmat2.data();
    std::valarray<float> ldata(m*n2);
    for (int k=0; k < m*n2; k++) 
	ldata[k] = ldat[k];
    oRSF left("left");
    left.put("n1",m);
    left.put("n2",n2);
    left.put("n3",1);
    left << ldata;

    FltNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );

    float *rdat = rmat.data();
    std::valarray<float> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = rdat[k];
    oRSF right;
    right.put("n1",n2);
    right.put("n2",n);
    right.put("n3",1);
    right << rdata;

    exit(0);
}
