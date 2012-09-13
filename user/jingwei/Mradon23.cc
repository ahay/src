//   special 2-D to 3-D Radon transform 
//   Input f(w,x) complex  --- f(tau,p)
//   Output u(tau,p1,p2) complex --- u(w,x1,x2)
//   Call bfio.setup23 bfio.kernel2 bfio.check2 bfio.eval2
//   In bfio.kernel2: fi=2 adjoint of hyper Radon
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
#include "serialize1.hh"

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
  FltNumVec x(nx);   
  for (int j=0; j<nx; j++)  x(j) = x0+j*dx;

  
  // Set output
  iRSF par(0);
  int ntau, np1, np2;
  par.get("ntau",ntau); 
  par.get("np1",np1); 
  par.get("np2",np2); 

  float tau0, dtau;
  par.get("tau0",tau0);
  par.get("dtau",dtau);  

  float p10, dp1;
  par.get("p10",p10);
  par.get("dp1",dp1);

  float p20, dp2;
  par.get("p20",p20);
  par.get("dp2",dp2);
 
  CpxNumMat u(ntau,np1*np2);  setvalue(u,cpx(0,0));
  FltNumVec tau(ntau);   
  for (int i=0; i<ntau; i++)  tau(i) = tau0+i*dtau;
  FltNumVec p(np1*np2);
  for (int j=0; j<np1; j++) 
    for (int k=0; k<np2; k++) {
      float pp1 = p10+j*dp1;
      float pp2 = p20+k*dp2;  
      p(np1*k+j) = sqrt(pp1*pp1+pp2*pp2);
    }

  oRSF output;
  output.put("n1",ntau);
  output.put("n2",np1);
  output.put("n3",np2);

  output.put("o1",tau0);
  output.put("d1",dtau);

  output.put("o2",p10);
  output.put("d2",dp1);

  output.put("o3",p20);
  output.put("d3",dp2);

  //output.type(SF_FLOAT);
  // this has be there if the input and output types are different


  // BFIO setup
  BFIO bfio("bfio_");
  iC( bfio.setup23(par,input) );

  int N;
  par.get("N",N); // number of partitions
  cerr<<"N "<<N<<endl;


  float time_eval;

  if(N<=256) {
    ck0 = clock();
    iC( bfio.eval2(N,f,w,x,u,tau,p) );
    ck1 = clock();    
    time_eval = float(ck1-ck0)/CLOCKS_PER_SEC;
  } else {
    tt0 = time(0);
    iC( bfio.eval2(N,f,w,x,u,tau,p) );
    tt1 = time(0);    
    time_eval = difftime(tt1,tt0);
  }
  //
  float relerr = 0;
  int NC = 64;
  ck0 = clock();
  iC( bfio.check2(N,f,w,x,u,tau,p,NC,relerr) );
  ck1 = clock();
  float time_chck = float(ck1-ck0)/CLOCKS_PER_SEC*float(ntau)*float(np1*np2)/float(NC);
  //
  cerr<<"Ta "<<time_eval<<endl;
  cerr<<"Td "<<time_chck<<endl;
  cerr<<"Rt "<<time_chck/time_eval<<endl;
  cerr<<"Ea "<<relerr<<endl;
  //
  
  std::valarray<sf_complex> udata(ntau*np1*np2);
  //std::valarray<float> udata(ntau*np1*np2);
  //udata.resize(ntau*np1*np2);
    
  for (int i=0; i<ntau; i++)
    for (int j=0; j<np1; j++)
      for (int k=0; k<np2; k++)
        udata[ntau*np1*k+ntau*j+i]=sf_cmplx(real(u(i,np1*k+j)),imag(u(i,np1*k+j)));
        //udata[ntau*np1*k+ntau*j+i]=real(u(i,np1*k+j));

  output << udata;

  exit(0);
}
