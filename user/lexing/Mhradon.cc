// 2-D hyperbolic Radon transform
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

using namespace std;
using std::cerr;

int main(int argc, char** argv)
{
  srand48(time(NULL));
  clock_t ck0, ck1;
  time_t tt0, tt1;
    
  sf_init(argc,argv); // Initialize RSF
    
  // Get input
  iRSF input;

  int nw, nx;
  input.get("n1",nw);
  input.get("n2",nx);

  float w0, dw;
  input.get("o1",w0);
  input.get("d1",dw);  

  float x0, dx;
  input.get("o2",x0);
  input.get("d2",dx);


  std::valarray<sf_complex> fdata(nw*nx);
  //fdata.resize(nw*nx);

  input >> fdata;

  CpxNumMat f(nw,nx);  setvalue(f,cpx(0,0));
  for (int i=0; i<nw; i++)
    for (int j=0; j<nx; j++)
      f(i,j) = cpx(crealf(fdata[nw*j+i]),cimagf(fdata[nw*j+i]));
  FltNumVec w(nw);
  for (int i=0; i<nw; i++)  w(i) = w0+i*dw;
  FltNumVec z(nx);   
  for (int j=0; j<nx; j++)  z(j) = x0+j*dx;
 

  // Set output
  iRSF par(0);
  int nt, np;
  par.get("nt",nt); 
  par.get("np",np); 

  float t0, dt;
  par.get("t0",t0);
  par.get("dt",dt);  

  float p0, dp;
  par.get("p0",p0);
  par.get("dp",dp);
 
  CpxNumMat u(nt,np);  setvalue(u,cpx(0,0));
  FltNumVec t(nt);   
  for (int i=0; i<nt; i++)  t(i) = t0+i*dt;
  FltNumVec p(np);
  for (int j=0; j<np; j++)  p(j) = p0+j*dp;

  oRSF output;
  output.put("n1",nt);
  output.put("n2",np);

  output.put("o1",t0);
  output.put("d1",dt);

  output.put("o2",p0);
  output.put("d2",dp);

  output.type(SF_FLOAT);

  // BFIO setup
  BFIO bfio("bfio_");
  iC( bfio.setup(par,input) );

  int N;
  par.get("N",N); // number of partitions
  cerr<<"N "<<N<<endl;


  float time_eval;

  if(N<=256) {
    ck0 = clock();
    iC( bfio.eval(N,f,w,z,u,t,p) );
    ck1 = clock();    
    time_eval = float(ck1-ck0)/CLOCKS_PER_SEC;
  } else {
    tt0 = time(0);
    iC( bfio.eval(N,f,w,z,u,t,p) );
    tt1 = time(0);    
    time_eval = difftime(tt1,tt0);
  }
  //
  float relerr = 0;
  int NC = 64;
  ck0 = clock();
  iC( bfio.check(N,f,w,z,u,t,p,NC,relerr) );
  ck1 = clock();
  float time_chck = float(ck1-ck0)/CLOCKS_PER_SEC*float(nt)*float(np)/float(NC);
  //
  cerr<<"Ta "<<time_eval<<endl;
  cerr<<"Td "<<time_chck<<endl;
  cerr<<"Rt "<<time_chck/time_eval<<endl;
  cerr<<"Ea "<<relerr<<endl;
  //
  
  std::valarray<float> udata(nt*np);
  //udata.resize(nt*np);
    
  for (int i=0; i<nt; i++)
    for (int j=0; j<np; j++)
      udata[nt*j+i]=real(u(i,j));

  output << udata;

  exit(0);
}
