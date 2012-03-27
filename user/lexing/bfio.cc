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
int BFIO::setup(iRSF& par, iRSF& inp)
{
  vector<int> all(1,1);
  
  par.get("EPSx1", _EPSx1); // number of chebyshev points
  par.get("EPSx2", _EPSx2);
  par.get("EPSp1", _EPSp1);
  par.get("EPSp2", _EPSp2);
  par.get("fi", _fi);
  par.get("EL", _EL);
  
  //cerr<<"EPSx1 "<<_EPSx1<<" EPSx2 "<<_EPSx2<<endl;
  //cerr<<"EPSp1 "<<_EPSp1<<" EPSp2 "<<_EPSp2<<endl;
  //cerr<<"fi "<<_fi<<endl;
  //cerr<<"EL "<<_EL<<endl;

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

  zmin = x0;
  zmax = x0+nx*dx;

  int nt;
  float t0, dt;
  par.get("nt",nt);
  par.get("t0",t0);
  par.get("dt",dt);
  tmin = t0;
  tmax = t0+nt*dt;

  int np;
  float p0, dp;
  par.get("np",np);
  par.get("p0",p0);
  par.get("dp",dp);
  pmin = p0;
  pmax = p0+np*dp;
  
  cerr<<"nw "<<nw<<" nx "<<nx<<endl;
  cerr<<"nt "<<nt<<" np "<<np<<endl;
  cerr<<"wmin "<<wmin<<" wmax "<<wmax<<endl;
  cerr<<"zmin "<<zmin<<" zmax "<<zmax<<endl;
  cerr<<"taumin "<<tmin<<" taumax "<<tmax<<endl;
  cerr<<"pmin "<<pmin<<" pmax "<<pmax<<endl;
  return 0;
}


//---------------------------------------
int BFIO::setup32(iRSF& par, iRSF& inp)
{
  vector<int> all(1,1);
  
  par.get("EPSx1", _EPSx1); // number of chebyshev points
  par.get("EPSx2", _EPSx2);
  par.get("EPSp1", _EPSp1);
  par.get("EPSp2", _EPSp2);
  par.get("fi", _fi);
  par.get("EL", _EL);
  
  //cerr<<"EPSx1 "<<_EPSx1<<" EPSx2 "<<_EPSx2<<endl;
  //cerr<<"EPSp1 "<<_EPSp1<<" EPSp2 "<<_EPSp2<<endl;
  //cerr<<"fi "<<_fi<<endl;
  //cerr<<"EL "<<_EL<<endl;

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

  float xmin = x0;
  float xmax = x0+nx*dx;
  float ymin = y0;
  float ymax = y0+ny*dy;
  zmin = 0.0;
  zmax = sqrt(xmax*xmax+ymax*ymax);

  int nt;
  float t0, dt;
  par.get("nt",nt);
  par.get("t0",t0);
  par.get("dt",dt);
  tmin = t0;
  tmax = t0+nt*dt;

  int np;
  float p0, dp;
  par.get("np",np);
  par.get("p0",p0);
  par.get("dp",dp);
  pmin = p0;
  pmax = p0+np*dp;
  
  cerr<<"nw "<<nw<<" nx*ny "<<nx*ny<<endl;
  cerr<<"nt "<<nt<<" np "<<np<<endl;
  cerr<<"wmin "<<wmin<<" wmax "<<wmax<<endl;
  cerr<<"xmin "<<xmin<<" xmax "<<xmax<<endl;
  cerr<<"ymin "<<ymin<<" ymax "<<ymax<<endl;
  cerr<<"zmin "<<zmin<<" zmax "<<zmax<<endl;
  cerr<<"taumin "<<tmin<<" taumax "<<tmax<<endl;
  cerr<<"pmin "<<pmin<<" pmax "<<pmax<<endl;
  return 0;
}

//---------------------------------------
int BFIO::setup23(iRSF& par, iRSF& inp)
{
  vector<int> all(1,1);
  
  par.get("EPSx1", _EPSx1); // number of chebyshev points
  par.get("EPSx2", _EPSx2);
  par.get("EPSp1", _EPSp1);
  par.get("EPSp2", _EPSp2);
  par.get("fi", _fi);
  par.get("EL", _EL);
  
  //cerr<<"EPSx1 "<<_EPSx1<<" EPSx2 "<<_EPSx2<<endl;
  //cerr<<"EPSp1 "<<_EPSp1<<" EPSp2 "<<_EPSp2<<endl;
  //cerr<<"fi "<<_fi<<endl;
  //cerr<<"EL "<<_EL<<endl;

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

  zmin = x0;
  zmax = x0+nx*dx;

  int nt;
  float t0, dt;
  par.get("nt",nt);
  par.get("t0",t0);
  par.get("dt",dt);
  tmin = t0;
  tmax = t0+nt*dt;

  int np1, np2;
  float p10, p20;
  float dp1, dp2;
  par.get("np1",np1);
  par.get("np2",np2);

  par.get("p10",p10);
  par.get("p20",p20);

  par.get("dp1",dp1);
  par.get("dp2",dp2);

  float p1min = p10;
  float p1max = p10+np1*dp1;

  float p2min = p20;
  float p2max = p20+np2*dp2;
  
  pmin = 0.0;
  pmax = sqrt(p1max*p1max+p2max*p2max);
  
  cerr<<"nw "<<nw<<" nx "<<nx<<endl;
  cerr<<"nt "<<nt<<" np1*np2 "<<np1*np2<<endl;
  cerr<<"wmin "<<wmin<<" wmax "<<wmax<<endl;
  cerr<<"zmin "<<zmin<<" zmax "<<zmax<<endl;
  cerr<<"taumin "<<tmin<<" taumax "<<tmax<<endl;
  cerr<<"p1min "<<p1min<<" p1max "<<p1max<<endl;
  cerr<<"p2min "<<p2min<<" p2max "<<p2max<<endl;
  cerr<<"pmin "<<pmin<<" pmax "<<pmax<<endl;
  return 0;
}



//---------------------------------------
int BFIO::kernel(int N, vector<Point2>& xs, vector<Point2>& ks, CpxNumMat& res)
{
  if(_fi==0) {
    //--------------------------
    int m = xs.size();
    int n = ks.size();
    vector<float> ts(m), ps(m);
    for(int i=0; i<m; i++)      ts[i] = xs[i](0)*(tmax-tmin) + tmin;
    for(int i=0; i<m; i++)      ps[i] = xs[i](1)*(pmax-pmin) + pmin;
    vector<float> ws(n), zs(n);
    for(int i=0; i<n; i++)      ws[i] = ks[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      zs[i] = ks[i](1)*(zmax-zmin) + zmin; 
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    for(int j=0; j<n; j++) {
      for(int i=0; i<m; i++) {
	float pz = ps[i]*zs[j];
	phs(i,j) = COEF * (sqrt(ts[i]*ts[i] + pz*pz)) * (ws[j]);
      }
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
  } else if(_fi==1) {
    //--------------------------
    int m = xs.size();
    int n = ks.size();
    vector<float> ts(m), ps(m);
    for(int i=0; i<m; i++)      ts[i] = xs[i](0)*(tmax-tmin) + tmin;
    for(int i=0; i<m; i++)      ps[i] = xs[i](1)*(pmax-pmin) + pmin;
    vector<float> ws(n), zs(n);
    for(int i=0; i<n; i++)      ws[i] = ks[i](0)*(wmax-wmin) + wmin;
    for(int i=0; i<n; i++)      zs[i] = ks[i](1)*(zmax-zmin) + zmin; 
    FltNumMat phs(m,n);
    float COEF = 2*M_PI;
    FltNumMat ss(m,n), cc(m,n);
    res.resize(m,n);
    for(int j=0; j<n; j++) {
      for(int i=0; i<m; i++) {
	float pz = ps[i]*zs[j];
        float in = ts[i]*ts[i] - pz*pz;
        if( in>0 ) {
	  phs(i,j) = COEF * (sqrt(in)) * (ws[j]);
          ss(i,j) = sin(phs(i,j));
	  cc(i,j) = cos(phs(i,j));
          res(i,j) = cpx( cc(i,j), ss(i,j) );
        }
        else {
	  //res(i,j) = cpx(0,0);
          phs(i,j) = 0;
          ss(i,j) = sin(phs(i,j));
	  cc(i,j) = cos(phs(i,j));
          res(i,j) = cpx( cc(i,j), ss(i,j) );
        }
      }
    }
  } else {
    //--------------------------
    iA(0);
  }
  return 0;
}

//---------------------------------------
int BFIO::check(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& z, const CpxNumMat& u, const FltNumVec& t, const FltNumVec& p, int NC, float& relerr)
{
  int N1 = f.m();
  int N2 = f.n();
  int M1 = u.m();
  int M2 = u.n();
  vector<Point2> src;
  for(int j=0; j<N2; j++)
    for(int i=0; i<N1; i++)
      src.push_back( Point2((w(i)-wmin)/(wmax-wmin), (z(j)-zmin)/(zmax-zmin)) );
  vector<cpx> app(NC);
  vector<cpx> dir(NC);
  //
  for(int g=0; g<NC; g++) {
    int x1 = int( floor(drand48()*M1) );
    int x2 = int( floor(drand48()*M2) );
    //
    app[g] = u(x1,x2);
    //
    vector<Point2> trg;  trg.push_back( Point2((t(x1)-tmin)/(tmax-tmin), (p(x2)-pmin)/(pmax-pmin)) );
    CpxNumMat res(1,N1*N2);  iC( kernel(N, trg, src, res) );
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

