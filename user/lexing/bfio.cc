// Butterfly algorithm
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

extern "C"
{
  void vdsqrt_(int* n, double*, double*);
  void vdsincos_(int* n, double*, double*, double*);
}

//---------------------------------------
int serialize(const Entry& e, ostream& os, const vector<int>& mask)
{
  iC( serialize(e._grid, os, mask) );
  iC( serialize(e._mats, os, mask) );
  iC( serialize(e._dirc, os, mask) );
  return 0;
}

int deserialize(Entry& e, istream& is, const vector<int>& mask)
{
  iC( deserialize(e._grid, is, mask) );
  iC( deserialize(e._mats, is, mask) );
  iC( deserialize(e._dirc, is, mask) );
  return 0;
}

#define CHEB_POINTS 13 /* number of chebyshev points */

//---------------------------------------
int BFIO::setup(iRSF &par, iRSF &inp) // map<string,string>& opts)
{
    vector<int> all(1,1);

    par.get("EPSx1",_EPSx1,CHEB_POINTS); 
    par.get("EPSx2",_EPSx2,CHEB_POINTS);
    par.get("EPSk1",_EPSk1,CHEB_POINTS);
    par.get("EPSk2",_EPSk2,CHEB_POINTS);
    par.get("fi",_fi);

    ifstream fin("bfio.bin");

    iC( deserialize(_e2dmap, fin, all) );

    int nw;
    float w0,dw;
    inp.get("n1",nw);
    inp.get("o1",w0);
    inp.get("d1",dw);
    wmin = w0;
    wmax = w0+(nw-1)*dw;

    int nx;
    float x0,dx;
    inp.get("n1",nx);
    inp.get("o1",x0);
    inp.get("d1",dx);
    zmin = x0;
    zmax = x0+(nx-1)*dx;

    int nt;
    float t0,dt;
    par.get("nt",nt);
    par.get("t0",t0);
    par.get("dt",dt);
    tmin = t0;
    tmax = t0+(nt-1)*dt;

    int np;
    float p0,dp;
    par.get("np",np);
    par.get("p0",p0);
    par.get("dp",dp);
    pmin = p0;
    pmax = p0+(np-1)*dp;

    return 0;
}


//---------------------------------------
int BFIO::kernel(int N, vector<Point2>& xs, vector<Point2>& ks, CpxNumMat& res)
{
  if(       _fi==0) {
    //--------------------------
    int m = xs.size();
    int n = ks.size();
    for(int i=0; i<m; i++)      xs[i] = xs[i]/double(N); //LEXING: IMPORTANT
    //
    vector<double> cs(m);
    for(int i=0; i<m; i++) {
      Point2 x = xs[i];
      cs[i] = (2+sin(2*M_PI*x(0))*sin(2*M_PI*x(1)))/3.0;
    }
    vector<double> rs(n);
    for(int j=0; j<n; j++) {
      Point2 k = ks[j];
      rs[j] = sqrt(k(0)*k(0) + k(1)*k(1));
    }
    DblNumMat phs(m,n);
    double COEF = 2*M_PI;
    for(int j=0; j<n; j++) {
      Point2 k = ks[j];
      for(int i=0; i<m; i++) {
	Point2 x = xs[i];
	phs(i,j) = COEF * ( (x(0)*k(0)+x(1)*k(1)) + cs[i]*rs[j] );
      }
    }
    DblNumMat ss(m,n), cc(m,n);
//    int TTL = m*n;
    //vdsincos_(&TTL, phs.data(), ss.data(), cc.data());
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
	//sincos(phs(i,j), &(ss(i,j)), &(cc(i,j)));
	ss(i,j) = sin(phs(i,j));
	cc(i,j) = cos(phs(i,j));
      }
    res.resize(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)
	res(i,j) = cpx( cc(i,j), ss(i,j) );
  } else if(_fi==1)  {
    //--------------------------
    int m = xs.size();
    int n = ks.size();
    for(int i=0; i<m; i++)      xs[i] = xs[i]/double(N); //LEXING: IMPORTANT
    //
    vector<double> tmp1(m), tmp2(m);
    for(int i=0; i<m; i++) {
      Point2 x = xs[i];
      tmp1[i] = (2+sin(2*M_PI*x(0))*sin(2*M_PI*x(1)))/3.0;
      tmp2[i] = (2+cos(2*M_PI*x(0))*cos(2*M_PI*x(1)))/3.0;
    }
    DblNumMat tmp(m,n);
    for(int j=0; j<n; j++) {
      Point2 k = ks[j];
      for(int i=0; i<m; i++) {
	tmp(i,j) = tmp1[i]*(k(0)*k(0)) + tmp2[i]*(k(1)*k(1));
      }
    }
//    int TTL = m*n;
    DblNumMat phs(m,n);
    //vdsqrt_(&TTL, tmp.data(), phs.data());
    for(int i=0; i<m; i++)      for(int j=0; j<n; j++)	phs(i,j) = sqrt(tmp(i,j));
    double COEF = 2*M_PI;
    for(int j=0; j<n; j++) {
      Point2 k = ks[j];
      for(int i=0; i<m; i++) {
	Point2 x = xs[i];
	phs(i,j) = COEF * ( (x(0)*k(0)+x(1)*k(1)) + phs(i,j) );
      }
    }
    DblNumMat ss(m,n), cc(m,n);
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
    //--------------------------
    int m = xs.size();
    int n = ks.size();
    for(int i=0; i<m; i++)      xs[i] = xs[i]/double(N); //LEXING: IMPORTANT
    DblNumMat phs(m,n);
    double COEF = 2*M_PI;
    for(int j=0; j<n; j++) {
      Point2 k = ks[j];
      for(int i=0; i<m; i++) {
	Point2 x = xs[i];
	phs(i,j) = COEF * ( x(0)*k(0)+x(1)*k(1) );
      }
    }
    DblNumMat ss(m,n), cc(m,n);
//    int TTL = m*n;
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
    //--------------------------
    int m = xs.size();
    int n = ks.size();

    vector<double> ts(m), ps(m);
    for(int i=0; i<m; i++)      ts[i] = (xs[i](0)/double(N))*(tmax-tmin) + tmin;
    for(int i=0; i<m; i++)      ps[i] = (xs[i](1)/double(N))*(pmax-pmin) + pmin;
    vector<double> ws(n), zs(n);
    for(int i=0; i<n; i++)      ws[i] = (ks[i](0)/double(N))*(wmax-wmin) + (wmax+wmin)/double(2);
    for(int i=0; i<n; i++)      zs[i] = (ks[i](1)/double(N))*(zmax-zmin) + (zmax+zmin)/double(2); //z is x
    DblNumMat phs(m,n);
    double COEF = 2*M_PI;
    for(int j=0; j<n; j++) {
      for(int i=0; i<m; i++) {
	double pz = ps[i]*zs[j];
	phs(i,j) = COEF * sqrt(ts[i]*ts[i] + pz*pz) * ws[j];
      }
    }
    DblNumMat ss(m,n), cc(m,n);
//    int TTL = m*n;
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
int BFIO::check(int N, const CpxNumMat& f, const CpxNumMat& u, int NC, double& relerr)
{
  int N1 = f.m();
  int N2 = f.n();
  int M1 = u.m();
  int M2 = u.n();
  vector<Point2> src;
  for(int j=0; j<N2; j++)
    for(int i=0; i<N1; i++)
      src.push_back( Point2(double(i*N)/double(N1)-double(N/2), double(j*N)/double(N2)-double(N/2)) );
  vector<cpx> app(NC);
  vector<cpx> dir(NC);
  //
  for(int g=0; g<NC; g++) {
    int x1 = int( floor(drand48()*M1) );
    int x2 = int( floor(drand48()*M2) );
    //
    app[g] = u(x1,x2);
    //
    vector<Point2> trg;    trg.push_back( Point2(double(x1*N)/double(M1),double(x2*N)/double(M2)) );
    CpxNumMat res(1,N1*N2); iC( kernel(N, trg, src, res) );
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
  double dn = 0;
  double en = 0;
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

