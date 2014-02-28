//   BFIO::setup2    
//   BFIO::setup32    
//   BFIO::setup3   
//
//   BFIO::kernel2   
//   BFIO::apkernel2 
//   BFIO::kernel3   
//   BFIO::dikernel3 
//   BFIO::kernel34  
//   
//   BFIO::check2    
//   BFIO::apcheck2  
//   BFIO::check3    
//   BFIO::check34   
//
//   Copyright (C) 2011 University of Texas at Austin
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

#include "bfio.hh"
#include "serialize.hh"

using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::set;
using std::queue;
using std::cerr;

//---------------------------------------
int serialize(const Entry& e, ostream& os, const vector<int>& mask)
{
  iC( serialize(e._grid, os, mask) );
  iC( serialize(e._mats, os, mask) );
  return 0;
}

int deserialize(Entry& e, istream& is, const vector<int>& mask)
{
  iC( deserialize(e._grid, is, mask) );
  iC( deserialize(e._mats, is, mask) );
  return 0;
}

//---------------------------------------
int BFIO::setup2(iRSF& par, iRSF& inp)
{
  vector<int> all(1,1);
  
  par.get("EPSx1", _EPSx1); // number of chebyshev points
  par.get("EPSx2", _EPSx2);
  par.get("EPSk1", _EPSk1);
  par.get("EPSk2", _EPSk2);
  par.get("fi", _fi);
  par.get("EL", _EL);

  ifstream fin("bfio.bin");

  iC( deserialize(_e2dmap, fin, all) );
  
  int nw, nx;
  inp.get("n1",nw);
  inp.get("n2",nx);

  float w0, dw;
  inp.get("o1",w0);
  inp.get("d1",dw);  

  float x0, dx;
  inp.get("o2",x0);
  inp.get("d2",dx);

  wmin = w0;
  wmax = w0+nw*dw;

  xmin = x0;
  xmax = x0+nx*dx;

  int ntau;
  float tau0, dtau;
  par.get("ntau",ntau);
  par.get("tau0",tau0);
  par.get("dtau",dtau);
  taumin = tau0;
  taumax = tau0+ntau*dtau;

  int np;
  float p0, dp;
  par.get("np",np);
  par.get("p0",p0);
  par.get("dp",dp);
  pmin = p0;
  pmax = p0+np*dp;
  
  cerr<<"nw "<<nw<<" nx "<<nx<<endl;
  cerr<<"wmin "<<wmin<<" wmax "<<wmax<<endl;
  cerr<<"xmin "<<xmin<<" xmax "<<xmax<<endl;
  cerr<<"ntau "<<ntau<<" np "<<np<<endl;
  cerr<<"taumin "<<taumin<<" taumax "<<taumax<<endl;
  cerr<<"pmin "<<pmin<<" pmax "<<pmax<<endl;

  cerr<<"fi "<<_fi<<endl;
  return 0;
}


//---------------------------------------
int BFIO::setup32(iRSF& par, iRSF& inp)
{
  vector<int> all(1,1);
  
  par.get("EPSx1", _EPSx1); // number of chebyshev points
  par.get("EPSx2", _EPSx2);
  par.get("EPSk1", _EPSk1);
  par.get("EPSk2", _EPSk2);
  par.get("fi", _fi);
  par.get("EL", _EL);
  
  ifstream fin("bfio.bin");

  iC( deserialize(_e2dmap, fin, all) );
  
  int nw, nx1, nx2;
  inp.get("n1",nw);
  inp.get("n2",nx1);
  inp.get("n3",nx2);

  float w0, dw;
  inp.get("o1",w0);
  inp.get("d1",dw);  

  float x10, dx1;
  inp.get("o2",x10);
  inp.get("d2",dx1);

  float x20, dx2;
  inp.get("o3",x20);
  inp.get("d3",dx2);

  wmin = w0;
  wmax = w0+nw*dw;

  float x1min = x10;
  float x1max = x10+nx1*dx1;
  float x2min = x20;
  float x2max = x20+nx2*dx2;
  xmin = 0.0;
  xmax = sqrt(x1max*x1max+x2max*x2max);

  int ntau;
  float tau0, dtau;
  par.get("ntau",ntau);
  par.get("tau0",tau0);
  par.get("dtau",dtau);
  taumin = tau0;
  taumax = tau0+ntau*dtau;

  int np;
  float p0, dp;
  par.get("np",np);
  par.get("p0",p0);
  par.get("dp",dp);
  pmin = p0;
  pmax = p0+np*dp;
  
  cerr<<"nw "<<nw<<" nx1*nx2 "<<nx1*nx2<<endl;
  cerr<<"wmin "<<wmin<<" wmax "<<wmax<<endl;
  cerr<<"x1min "<<x1min<<" x1max "<<x1max<<endl;
  cerr<<"x2min "<<x2min<<" x2max "<<x2max<<endl;
  cerr<<"xmin "<<xmin<<" xmax "<<xmax<<endl;
  cerr<<"ntau "<<ntau<<" np "<<np<<endl;
  cerr<<"taumin "<<taumin<<" taumax "<<taumax<<endl;
  cerr<<"pmin "<<pmin<<" pmax "<<pmax<<endl;

  cerr<<"fi "<<_fi<<endl;
  return 0;
}


//---------------------------------------
int BFIO::setup3(iRSF& par, iRSF& inp)
{
  vector<int> all(1,1);
  
  par.get("EPSx1", _EPSx1); // number of chebyshev points
  par.get("EPSx2", _EPSx2);
  par.get("EPSx3", _EPSx3);
  par.get("EPSk1", _EPSk1);
  par.get("EPSk2", _EPSk2);
  par.get("EPSk3", _EPSk3);
  par.get("fi", _fi);
  par.get("EL", _EL);

  ifstream fin("bfio.bin");

  iC( deserialize(_e2dmap, fin, all) );
  
  int nw, nx, ny;
  inp.get("n1",nw);
  inp.get("n2",nx);
  inp.get("n3",ny);

  float w0, dw;
  inp.get("o1",w0);
  inp.get("d1",dw);  

  float x0, dx;
  inp.get("o2",x0);
  inp.get("d2",dx);

  float y0, dy;
  inp.get("o3",y0);
  inp.get("d3",dy);

  wmin = w0;
  wmax = w0+nw*dw;

  xmin = x0;
  xmax = x0+nx*dx;

  ymin = y0;
  ymax = y0+ny*dy;

  int ntau;
  float tau0, dtau;
  par.get("ntau",ntau);
  par.get("tau0",tau0);
  par.get("dtau",dtau);
  taumin = tau0;
  taumax = tau0+ntau*dtau;

  int np, nq;
  float p0, q0;
  float dp, dq;
  par.get("np",np);
  par.get("nq",nq);

  par.get("p0",p0);
  par.get("q0",q0);

  par.get("dp",dp);
  par.get("dq",dq);

  pmin = p0;
  pmax = p0+np*dp;
  qmin = q0;
  qmax = q0+nq*dq; 
   
  cerr<<"nw "<<nw<<" nx "<<nx<<" ny "<<ny<<endl;
  cerr<<"wmin "<<wmin<<" wmax "<<wmax<<endl;
  cerr<<"xmin "<<xmin<<" xmax "<<xmax<<endl;
  cerr<<"ymin "<<ymin<<" ymax "<<ymax<<endl;
  cerr<<"ntau "<<ntau<<" np "<<np<<" nq "<<nq<<endl;
  cerr<<"taumin "<<taumin<<" taumax "<<taumax<<endl;
  cerr<<"pmin "<<pmin<<" pmax "<<pmax<<endl;
  cerr<<"qmin "<<qmin<<" qmax "<<qmax<<endl;

  cerr<<"fi "<<_fi<<endl;
  return 0;
}

//---------------------------------------
int BFIO::kernel2(int N, vector<Point2>& trg, vector<Point2>& src, CpxNumMat& res)
{
  if(_fi==1) {
    // hyperbolic Radon 
    //--------------------------
    int m = trg.size();
    int n = src.size();
    vector<float> taus(m), ps(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    vector<float> ws(n), xs(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
	float pz = ps[i]*xs[j];
	phs(i,j) = COEF * (sqrt(taus[i]*taus[i] + pz*pz)) * (ws[j]);
      }
    FltNumMat ss(m,n), cc(m,n);
    //int TTL = m*n;
    //vdsincos_(&TTL, phs.data(), ss.data(), cc.data());
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
	//sincos(phs(i,j), &(ss(i,j)), &(cc(i,j)));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else if(_fi==2) {
    // adjoint of hyperbolic Radon 
    //--------------------------
    int m = trg.size();
    int n = src.size();
    vector<float> taus(m), ps(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    vector<float> ws(n), xs(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
	float pz = ps[i]*xs[j];
	phs(i,j) = -COEF * (sqrt(ws[j]*ws[j] + pz*pz)) * (taus[i]);
      }
    FltNumMat ss(m,n), cc(m,n);
    //int TTL = m*n;
    //vdsincos_(&TTL, phs.data(), ss.data(), cc.data());
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
	//sincos(phs(i,j), &(ss(i,j)), &(cc(i,j)));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );   
  } else if(_fi==3) {
    // x*k
    //--------------------------
    int m = trg.size();
    int n = src.size();
    vector<float> taus(m), ps(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    vector<float> ws(n), xs(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
	phs(i,j) = COEF * (taus[i]*ws[j]+ps[i]*xs[j]);
      }
    FltNumMat ss(m,n), cc(m,n);
    //int TTL = m*n;
    //vdsincos_(&TTL, phs.data(), ss.data(), cc.data());
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
	//sincos(phs(i,j), &(ss(i,j)), &(cc(i,j)));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );   
  } else if(_fi==4) {
    // -x*k
    //--------------------------
    int m = trg.size();
    int n = src.size();
    vector<float> taus(m), ps(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    vector<float> ws(n), xs(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
	phs(i,j) = -COEF * (taus[i]*ws[j]+ps[i]*xs[j]);
      }
    FltNumMat ss(m,n), cc(m,n);
    //int TTL = m*n;
    //vdsincos_(&TTL, phs.data(), ss.data(), cc.data());
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
	//sincos(phs(i,j), &(ss(i,j)), &(cc(i,j)));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );   
  } else {
    //--------------------------
    iA(0);
  }
  return 0;
}

//---------------------------------------
int BFIO::apkernel2(int N, vector<Point2>& trg, vector<Point2>& src, CpxNumMat& res, const float xx)
{
  if(_fi==1) {
    // apex shifted hyperbolic Radon 
    //--------------------------
    int m = trg.size();
    int n = src.size();
    vector<float> taus(m), ps(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    vector<float> ws(n), xs(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
	float pz = ps[i]*(xs[j]-xx);
	phs(i,j) = COEF * (sqrt(taus[i]*taus[i] + pz*pz)) * (ws[j]);
      }
    FltNumMat ss(m,n), cc(m,n);
    //int TTL = m*n;
    //vdsincos_(&TTL, phs.data(), ss.data(), cc.data());
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
	//sincos(phs(i,j), &(ss(i,j)), &(cc(i,j)));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else {
    //--------------------------
    iA(0);
  }
  return 0;
}

//---------------------------------------
int BFIO::kernel3(int N, vector<Point3>& trg, vector<Point3>& src, CpxNumMat& res)
{
  if(_fi==0) {
    // linear Radon 
    //--------------------------
    int m = trg.size();
    int n = src.size();
    //
    vector<float> taus(m), ps(m), qs(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    for(int i=0; i<m; i++)      qs[i] = trg[i](2)*(qmax-qmin) + qmin;
    //
    vector<float> ws(n), xs(n), ys(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    for(int i=0; i<n; i++)      ys[i] = src[i](2)*(ymax-ymin) + ymin; 
    //
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
	phs(i,j) = COEF * (taus[i] + ps[i]*xs[j] + qs[i]*ys[j]) * ws[j];
      }
    FltNumMat ss(m,n), cc(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else if(_fi==1) {
    // reflection Radon
    // ws --> w; xs --> ax; ys --> ay
    // taus --> tau0; ps --> a0x; qs --> a0y
    int m = trg.size();
    int n = src.size();
    //
    vector<float> taus(m), ps(m), qs(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    for(int i=0; i<m; i++)      qs[i] = trg[i](2)*(qmax-qmin) + qmin;
    //
    vector<float> ws(n), xs(n), ys(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    for(int i=0; i<n; i++)      ys[i] = src[i](2)*(ymax-ymin) + ymin; 
    //
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    vector<float> tana0x(m), tana0y(m);
    vector<float> tanax(n), tanay(n);
    for(int i=0; i<m; i++) {
      tana0x[i] = tan(ps[i]*M_PI/180);
      tana0y[i] = tan(qs[i]*M_PI/180);
    }
    for(int i=0; i<n; i++) {
      tanax[i] = tan(xs[i]*M_PI/180);
      tanay[i] = tan(ys[i]*M_PI/180);
    }
    float a, b, c;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
        a = sqrt(1 + tana0x[i]*tana0x[i] + tana0y[i]*tana0y[i]);
        b = sqrt(1 + tanax[j]*tanax[j] + tanay[j]*tanay[j]);
        c = tana0x[i]*tanax[j] + tana0y[i]*tanay[j];
	phs(i,j) = COEF * ws[j] * taus[i] / (a*b-c);
      }
    FltNumMat ss(m,n), cc(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else if(_fi==2) {
    // diffraction Radon
    // ws --> w; xs --> ax; ys --> ay
    // taus --> tau0; ps --> kix; qs --> kiy
    int m = trg.size();
    int n = src.size();
    //
    vector<float> taus(m), ps(m), qs(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    for(int i=0; i<m; i++)      qs[i] = trg[i](2)*(qmax-qmin) + qmin;
    //
    vector<float> ws(n), xs(n), ys(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    for(int i=0; i<n; i++)      ys[i] = src[i](2)*(ymax-ymin) + ymin; 
    //
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    vector<float> tanax(n), tanay(n);
    for(int i=0; i<n; i++) {
      tanax[i] = tan(xs[i]*M_PI/180);
      tanay[i] = tan(ys[i]*M_PI/180);
    }
    float K1, K;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
        K1 = ps[i]*tanax[j] + qs[i]*tanay[j];
        K = sqrt(K1*K1 + ps[i]*ps[i]+ qs[i]*qs[i] + 1);
	phs(i,j) = COEF * ws[j] * taus[i] * (K1+K);
      }
    FltNumMat ss(m,n), cc(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else if(_fi==3) {
    // adjoint of reflection Radon 
    // only add - to phi
    //--------------------------
    int m = trg.size();
    int n = src.size();
    //
    vector<float> taus(m), ps(m), qs(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    for(int i=0; i<m; i++)      qs[i] = trg[i](2)*(qmax-qmin) + qmin;
    //
    vector<float> ws(n), xs(n), ys(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    for(int i=0; i<n; i++)      ys[i] = src[i](2)*(ymax-ymin) + ymin; 
    //
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    vector<float> tana0x(m), tana0y(m);
    vector<float> tanax(n), tanay(n);
    for(int i=0; i<m; i++) {
      tana0x[i] = tan(ps[i]*M_PI/180);
      tana0y[i] = tan(qs[i]*M_PI/180);
    }
    for(int i=0; i<n; i++) {
      tanax[i] = tan(xs[i]*M_PI/180);
      tanay[i] = tan(ys[i]*M_PI/180);
    }
    float a, b, c;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
        a = sqrt(1 + tana0x[i]*tana0x[i] + tana0y[i]*tana0y[i]);
        b = sqrt(1 + tanax[j]*tanax[j] + tanay[j]*tanay[j]);
        c = tana0x[i]*tanax[j] + tana0y[i]*tanay[j];
	phs(i,j) = -COEF * ws[j] * taus[i] / (a*b-c);
      }
    FltNumMat ss(m,n), cc(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else if(_fi==4) {
    // adjoint of diffraction Radon 
    //--------------------------
    int m = trg.size();
    int n = src.size();
    //
    vector<float> taus(m), ps(m), qs(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    for(int i=0; i<m; i++)      qs[i] = trg[i](2)*(qmax-qmin) + qmin;
    //
    vector<float> ws(n), xs(n), ys(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    for(int i=0; i<n; i++)      ys[i] = src[i](2)*(ymax-ymin) + ymin; 
    //
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    vector<float> tanax(m), tanay(m);
    for(int i=0; i<m; i++) {
      tanax[i] = tan(ps[i]*M_PI/180);
      tanay[i] = tan(qs[i]*M_PI/180);
    }
    float K1, K;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
        K1 = xs[j]*tanax[i] + ys[j]*tanay[i];
        K = sqrt(K1*K1 + xs[j]*xs[j]+ ys[j]*ys[j] + 1);
	phs(i,j) = -COEF * ws[j] * taus[i] * (K1+K);
      }
    FltNumMat ss(m,n), cc(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else {
    //--------------------------
    iA(0);
  }
  return 0;
}

//---------------------------------------
int BFIO::dikernel3(const int fi, const float tau, const float p, const float q, const float x, const float y, float& t)
{
  if(fi==1) {
    // reflection Radon
    // tau --> tau0; p --> a0x; q --> a0y
    // t --> tau; x --> ax; y --> ay
    float tana0x = tan(p*M_PI/180); 
    float tana0y = tan(q*M_PI/180);
    float tanax = tan(x*M_PI/180);
    float tanay = tan(y*M_PI/180);
    float a = sqrt(1 + tana0x*tana0x + tana0y*tana0y);
    float b = sqrt(1 + tanax*tanax + tanay*tanay);
    float c = tana0x*tanax + tana0y*tanay;
    t = tau / (a*b-c);
  } else if(fi==2) {
    // diffraction Radon
    // tau --> tau0; p --> kix; q --> kiy
    // t --> tau; x --> ax; y --> ay
    float tanax = tan(x*M_PI/180);
    float tanay = tan(y*M_PI/180);
    float K1 = p*tanax + q*tanay;
    float K = sqrt(K1*K1 + p*p + q*q + 1);
    t = tau * (K1+K);
  } else {
    //--------------------------
    iA(0);
  }
  return 0;
}

//---------------------------------------
int BFIO::kernel34(int N, vector<Point3>& trg, vector<Point3>& src, CpxNumMat& res, const float xx)
{
  if(_fi==1) {
    // ps=sqrt(W11);  qs=sqrt(W22);  xx=W12;
    //--------------------------
    int m = trg.size();
    int n = src.size();
    //
    vector<float> taus(m), ps(m), qs(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    for(int i=0; i<m; i++)      qs[i] = trg[i](2)*(qmax-qmin) + qmin;
    //
    vector<float> ws(n), xs(n), ys(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    for(int i=0; i<n; i++)      ys[i] = src[i](2)*(ymax-ymin) + ymin; 
    //
    FltNumMat phs(m,n);
    FltNumMat ss(m,n), cc(m,n);
    res.resize(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
        phs(i,j) = taus[i]*taus[i] + ps[i]*ps[i]*xs[j]*xs[j] + qs[i]*qs[i]*ys[j]*ys[j] + 2*xx*xs[j]*ys[j];
        if (phs(i,j)>=0) {
	  phs(i,j) = COEF * sqrt(phs(i,j)) * ws[j];
          ss(i,j) = sin(phs(i,j));
	  cc(i,j) = cos(phs(i,j));
          res(i,j) = cpx( cc(i,j), ss(i,j) );
	} else {
          res(i,j) = cpx(0,0);
        }
      }
 } else if(_fi==2) {
    // ps=Wcos;  qs=Wsin;  xx=Wavg;
    //--------------------------
    int m = trg.size();
    int n = src.size();
    //
    vector<float> taus(m), ps(m), qs(m);
    for(int i=0; i<m; i++)      taus[i] = trg[i](0)*(taumax-taumin) + taumin;
    for(int i=0; i<m; i++)      ps[i] = trg[i](1)*(pmax-pmin) + pmin;
    for(int i=0; i<m; i++)      qs[i] = trg[i](2)*(qmax-qmin) + qmin;
    //
    vector<float> ws(n), xs(n), ys(n);
    for(int i=0; i<n; i++)      ws[i] = src[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      xs[i] = src[i](1)*(xmax-xmin) + xmin; 
    for(int i=0; i<n; i++)      ys[i] = src[i](2)*(ymax-ymin) + ymin; 
    //
    FltNumMat phs(m,n);
    FltNumMat ss(m,n), cc(m,n);
    res.resize(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) 
      for(int i=0; i<m; i++) {
	phs(i,j) = taus[i]*taus[i] + xx*(xs[j]*xs[j]+ys[j]*ys[j]) + ps[i]*(xs[j]*xs[j]-ys[j]*ys[j]) + 2*qs[i]*(xs[j]*ys[j]);
        if (phs(i,j)>=0) {
	  phs(i,j) = COEF * sqrt(phs(i,j)) * ws[j];
          ss(i,j) = sin(phs(i,j));
	  cc(i,j) = cos(phs(i,j));
          res(i,j) = cpx( cc(i,j), ss(i,j) );
	} else {
          res(i,j) = cpx(0,0);
        }
      }
  } else {
    //--------------------------
    iA(0);
  }
  return 0;
}


//---------------------------------------
int BFIO::check2(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& x, const CpxNumMat& u, const FltNumVec& tau, const FltNumVec& p, int NC, float& relerr)
{
  int N1 = f.m();
  int N2 = f.n();
  int M1 = u.m();
  int M2 = u.n();
  vector<Point2> src;
  for(int j=0; j<N2; j++)
    for(int i=0; i<N1; i++)
      src.push_back( Point2((w(i)-wmin)/(wmax-wmin), (x(j)-xmin)/(xmax-xmin)) );
  vector<cpx> app(NC);
  vector<cpx> dir(NC);
  //
  for(int g=0; g<NC; g++) {
    int x1 = int( floor(drand48()*M1) );
    int x2 = int( floor(drand48()*M2) );
    //
    app[g] = u(x1,x2);
    //
    vector<Point2> trg;  trg.push_back( Point2((tau(x1)-taumin)/(taumax-taumin), (p(x2)-pmin)/(pmax-pmin)) );
    CpxNumMat res(1,N1*N2);  iC( kernel2(N, trg, src, res) );
    CpxNumMat resaux(N1,N2,false,res.data());
    cpx ttl(0,0);
    for(int j=0; j<N2; j++)
      for(int i=0; i<N1; i++)
	ttl = ttl + resaux(i,j) * f(i,j);
    dir[g] = ttl;
  }
  vector<cpx> err(NC);
  for(int g=0; g<NC; g++)
    err[g] = app[g] - dir[g];
  float dn = 0;
  float en = 0;
  for(int g=0; g<NC; g++) {
    dn += abs(dir[g])*abs(dir[g]);
    en += abs(err[g])*abs(err[g]);
  }
  dn = sqrt(dn);
  en = sqrt(en);
  relerr = en/dn;
  //
  return 0;
}

//---------------------------------------
int BFIO::apcheck2(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& x, const CpxNumMat& u, const FltNumVec& tau, const FltNumVec& p, const float xx, int NC, float& relerr)
{
  int N1 = f.m();
  int N2 = f.n();
  int M1 = u.m();
  int M2 = u.n();
  vector<Point2> src;
  for(int j=0; j<N2; j++)
    for(int i=0; i<N1; i++)
      src.push_back( Point2((w(i)-wmin)/(wmax-wmin), (x(j)-xmin)/(xmax-xmin)) );
  vector<cpx> app(NC);
  vector<cpx> dir(NC);
  //
  for(int g=0; g<NC; g++) {
    int x1 = int( floor(drand48()*M1) );
    int x2 = int( floor(drand48()*M2) );
    //
    app[g] = u(x1,x2);
    //
    vector<Point2> trg;  trg.push_back( Point2((tau(x1)-taumin)/(taumax-taumin), (p(x2)-pmin)/(pmax-pmin)) );
    CpxNumMat res(1,N1*N2);  iC( apkernel2(N, trg, src, res, xx) );
    CpxNumMat resaux(N1,N2,false,res.data());
    cpx ttl(0,0);
    for(int j=0; j<N2; j++)
      for(int i=0; i<N1; i++)
	ttl = ttl + resaux(i,j) * f(i,j);
    dir[g] = ttl;
  }
  vector<cpx> err(NC);
  for(int g=0; g<NC; g++)
    err[g] = app[g] - dir[g];
  float dn = 0;
  float en = 0;
  for(int g=0; g<NC; g++) {
    dn += abs(dir[g])*abs(dir[g]);
    en += abs(err[g])*abs(err[g]);
  }
  dn = sqrt(dn);
  en = sqrt(en);
  relerr = en/dn;
  //
  return 0;
}

//---------------------------------------
int BFIO::check3(int N, const CpxNumTns& f, const FltNumVec& w, const FltNumVec& x, const FltNumVec& y, const CpxNumTns& u, const FltNumVec& tau, const FltNumVec& p, const FltNumVec& q, int NC, float& relerr)
{
  int N1 = f.m();
  int N2 = f.n();
  int N3 = f.p();
  int M1 = u.m();
  int M2 = u.n();
  int M3 = u.p();
  //
  vector<Point3> src;
  for(int k=0; k<N3; k++)
    for(int j=0; j<N2; j++)
      for(int i=0; i<N1; i++)
        src.push_back( Point3((w(i)-wmin)/(wmax-wmin), (x(j)-xmin)/(xmax-xmin), (y(k)-ymin)/(ymax-ymin)) );
  //
  vector<cpx> app(NC);
  vector<cpx> dir(NC);
  //
  for(int g=0; g<NC; g++) {
    int x1 = int( floor(drand48()*M1) );
    int x2 = int( floor(drand48()*M2) );
    int x3 = int( floor(drand48()*M3) );
    //
    app[g] = u(x1,x2,x3);
    //
    vector<Point3> trg;  trg.push_back( Point3((tau(x1)-taumin)/(taumax-taumin), (p(x2)-pmin)/(pmax-pmin), (q(x3)-qmin)/(qmax-qmin)) );
    CpxNumMat res(1,N1*N2*N3);  iC( kernel3(N, trg, src, res) );
    CpxNumTns resaux(N1,N2,N3,false,res.data());
    cpx ttl(0,0);
    for(int k=0; k<N3; k++)
      for(int j=0; j<N2; j++)
        for(int i=0; i<N1; i++)
	  ttl = ttl + resaux(i,j,k) * f(i,j,k);
    dir[g] = ttl;
  }
  vector<cpx> err(NC);
  for(int g=0; g<NC; g++)
    err[g] = app[g] - dir[g];
  float dn = 0;
  float en = 0;
  for(int g=0; g<NC; g++) {
    dn += abs(dir[g])*abs(dir[g]);
    en += abs(err[g])*abs(err[g]);
  }
  dn = sqrt(dn);
  en = sqrt(en);
  relerr = en/dn;
  //
  return 0;
}

//---------------------------------------
int BFIO::check34(int N, const CpxNumTns& f, const FltNumVec& w, const FltNumVec& x, const FltNumVec& y, const CpxNumTns& u, const FltNumVec& tau, const FltNumVec& p, const FltNumVec& q, const float xx, int NC, float& relerr)
{
  int N1 = f.m();
  int N2 = f.n();
  int N3 = f.p();
  int M1 = u.m();
  int M2 = u.n();
  int M3 = u.p();
  //
  vector<Point3> src;
  for(int k=0; k<N3; k++)
    for(int j=0; j<N2; j++)
      for(int i=0; i<N1; i++)
        src.push_back( Point3((w(i)-wmin)/(wmax-wmin), (x(j)-xmin)/(xmax-xmin), (y(k)-ymin)/(ymax-ymin)) );
  //
  vector<cpx> app(NC);
  vector<cpx> dir(NC);
  //
  for(int g=0; g<NC; g++) {
    int x1 = int( floor(drand48()*M1) );
    int x2 = int( floor(drand48()*M2) );
    int x3 = int( floor(drand48()*M3) );
    //
    app[g] = u(x1,x2,x3);
    //
    vector<Point3> trg;  trg.push_back( Point3((tau(x1)-taumin)/(taumax-taumin), (p(x2)-pmin)/(pmax-pmin), (q(x3)-qmin)/(qmax-qmin)) );
    CpxNumMat res(1,N1*N2*N3);  iC( kernel34(N, trg, src, res, xx) );
    CpxNumTns resaux(N1,N2,N3,false,res.data());
    cpx ttl(0,0);
    for(int k=0; k<N3; k++)
      for(int j=0; j<N2; j++)
        for(int i=0; i<N1; i++)
	  ttl = ttl + resaux(i,j,k) * f(i,j,k);
    dir[g] = ttl;
  }
  vector<cpx> err(NC);
  for(int g=0; g<NC; g++)
    err[g] = app[g] - dir[g];
  float dn = 0;
  float en = 0;
  for(int g=0; g<NC; g++) {
    dn += abs(dir[g])*abs(dir[g]);
    en += abs(err[g])*abs(err[g]);
  }
  dn = sqrt(dn);
  en = sqrt(en);
  relerr = en/dn;
  //
  return 0;
}




