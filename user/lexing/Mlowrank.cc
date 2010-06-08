// Lowrank approximation for wave propagation. 

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

#include <rsf.hh>

#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

int optionsCreate(int argc, char** argv, map<string,string>& options)
{
  options.clear();
  for(int k=1; k<argc; k=k+2) {
    options[ string(argv[k]) ] = string(argv[k+1]);
  }
  return 0;
}

//
int nz=1200;
double dz=0.00762;
int nx=2133;
double dx=0.0143;
double dt = 0.001;
int nkz=1200;
double dkz=0.109361;
double kz0=-65.6168;
int nkx=2160;
double dkx=0.032375;
double kx0=-34.965;

DblNumVec vs;
DblNumVec ks;

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
  int nr = rs.size();
  int nc = cs.size();
  res.resize(nr,nc);  setvalue(res,0.0);
  for(int a=0; a<nr; a++)
    for(int b=0; b<nc; b++) {
      res(a,b) = cos(vs(rs[a])*ks(cs[b])*dt)-1.0;
    }
  return 0;
}

int main(int argc, char** argv)
{
  clock_t ck0, ck1;
  time_t t0, t1;
  srand48(time(NULL));
  map<string,string> opts;
  optionsCreate(argc, argv, opts);
  
  //1. read stuff
  map<string,string>::iterator mi;
  vector<int> all(1,1);
  //
  mi = opts.find("-vsfile");  assert(mi!=opts.end());
  char vsfile[100];  {istringstream ss((*mi).second);  ss>>vsfile;}
  //DblNumVec vs;
  {ifstream fin(vsfile);  iC( deserialize(vs, fin, all) );  }
  cerr<<vs.m()<<endl;
  iA(vs.m()==nx*nz);
  //
  ks.resize(nkz*nkx);  setvalue(ks,0.0);
  int cnt = 0;
  for(int ix=0; ix<nkx; ix++)
    for(int iz=0; iz<nkz; iz++) {
      double tmpz = kz0 + iz*dkz;
      double tmpx = kx0 + ix*dkx;
      double tmp = sqrt(tmpz*tmpz+tmpx*tmpx);
      ks(cnt) = tmp;      cnt++;
    }
  //
  //cerr<<sqrt(energy(vs))<<" "<<sqrt(energy(ks))<<endl;
  int m = vs.m();
  int n = ks.m();
  double eps = 1e-4;
  int npk = 20;
  vector<int> cidx;
  vector<int> ridx;
  DblNumMat mid;
  t0 = time(0);
  iC( lowrank(m,n,&sample,eps,npk,cidx,ridx,mid) );
  t1 = time(0);  cout<<"lowrank used "<<difftime(t1,t0)<<"secs "<<endl;
  for(int k=0; k<cidx.size(); k++)
    cerr<<cidx[k]<<" ";
  cerr<<endl;
  for(int k=0; k<ridx.size(); k++)
    cerr<<ridx[k]<<" ";
  cerr<<endl;
  
  return 0;
}





