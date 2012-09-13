//   Slow 3-D to 3-D linear Radon transform
//   Input f(w,x,y) complex
//   Output u(w,p,q) complex
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

  sf_init(argc,argv); // Initialize RSF
    
  // Get input
  iRSF input;

  int nw, nx, ny;
  input.get("n1",nw);
  input.get("n2",nx);
  input.get("n3",ny);

  float w0, dw;
  input.get("o1",w0);
  input.get("d1",dw);  

  float x0, dx;
  input.get("o2",x0);
  input.get("d2",dx);

  float y0, dy;
  input.get("o3",y0);
  input.get("d3",dy);


  std::valarray<sf_complex> fdata(nw*nx*ny);

  input >> fdata;

  CpxNumTns f(nw,nx,ny);    setvalue(f,cpx(0,0));
  for (int i=0; i<nw; i++)
    for (int j=0; j<nx; j++)
      for (int k=0; k<ny; k++)
        f(i,j,k) = cpx(crealf(fdata[nw*nx*k+nw*j+i]),cimagf(fdata[nw*nx*k+nw*j+i]));
  FltNumVec w(nw);
  for (int i=0; i<nw; i++)  w(i) = w0+i*dw;
  FltNumVec x(nx);   
  for (int j=0; j<nx; j++)  x(j) = x0+j*dx;
  FltNumVec y(ny);   
  for (int k=0; k<ny; k++)  y(k) = y0+k*dy;

  
  // Set output
  iRSF par(0);
  int np, nq;
  par.get("np",np); 
  par.get("nq",nq); 

  float p0, dp;
  par.get("p0",p0);
  par.get("dp",dp);

  float q0, dq;
  par.get("q0",q0);
  par.get("dq",dq);
 
  CpxNumTns u(nw,np,nq);  setvalue(u,cpx(0,0));
  FltNumVec p(np);
  for (int j=0; j<np; j++)  p(j) = p0+j*dp;
  FltNumVec q(nq);
  for (int k=0; k<nq; k++)  q(k) = q0+k*dq;

  oRSF output;
  output.put("n1",nw);
  output.put("n2",np);
  output.put("n3",nq);

  output.put("o1",w0);
  output.put("d1",dw);

  output.put("o2",p0);
  output.put("d2",dp);

  output.put("o3",q0);
  output.put("d3",dq);

  clock_t ck0, ck1; 
  time_t tt0, tt1; 
  float clockt_eval, timet_eval;

  ck0 = clock();
  tt0 = time(0);
  for(int i=0; i<nw; i++)
    for(int j=0; j<np; j++)
      for(int k=0; k<nq; k++)
        for(int m=0; m<nx; m++)
	  for(int n=0; n<ny; n++) {
	    float phs = 2*M_PI*w(i)*(p(j)*x(m)+q(k)*y(n));
            cpx res(0,0);
            float cc=cos(phs);
            float ss=sin(phs);
            res = cpx(cc,ss);
            u(i,j,k) += res*f(i,m,n);
	  }
  ck1 = clock();
  tt1 = time(0);
  clockt_eval = float(ck1-ck0)/CLOCKS_PER_SEC;   
  timet_eval = difftime(tt1,tt0);
  cerr<<"slow lradon3 clockt_eval"<<clockt_eval<<endl;
  cerr<<"slow lradon3 timet_eval"<<timet_eval<<endl;


  std::valarray<sf_complex> udata(nw*np*nq);
  
  for (int i=0; i<nw; i++)
    for (int j=0; j<np; j++)
      for (int k=0; k<nq; k++) 
        udata[nw*np*k+nw*j+i]=sf_cmplx(real(u(i,j,k)),imag(u(i,j,k)));

  output << udata;

  exit(0);
}
