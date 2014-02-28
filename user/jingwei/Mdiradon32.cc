//   direct azimuthally isotropic 3to2 hyper Radon transform
//   real f(t,x1,x2) --> real u(tau,p)
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

  int nt, nx1, nx2;
  input.get("n1",nt);
  input.get("n2",nx1);
  input.get("n3",nx2);

  float t0, dt;
  input.get("o1",t0);
  input.get("d1",dt);  

  float x10, dx1;
  input.get("o2",x10);
  input.get("d2",dx1);

  float x20, dx2;
  input.get("o3",x20);
  input.get("d3",dx2);

  std::valarray<float> fdata(nt*nx1*nx2);

  input >> fdata;

  // Set output
  iRSF par(0);

  int ntau, np;
  par.get("ntau",ntau); 
  par.get("np",np); 

  float tau0, dtau;
  par.get("tau0",tau0);
  par.get("dtau",dtau);  

  float p0, dp;
  par.get("p0",p0);
  par.get("dp",dp);

  oRSF output;
  output.put("n1",ntau);
  output.put("n2",np);
  output.put("n3",1);

  output.put("o1",tau0);
  output.put("d1",dtau);

  output.put("o2",p0);
  output.put("d2",dp);

  std::valarray<float> udata(ntau*np);

  cerr<<"tmin "<<t0<<" tmax "<<t0+nt*dt<<endl;
  cerr<<"x1min "<<x10<<" x1max "<<x10+nx1*dx1<<endl;
  cerr<<"x2min "<<x20<<" x2max "<<x20+nx2*dx2<<endl;
  cerr<<"taumin "<<tau0<<" taumax "<<tau0+ntau*dtau<<endl;
  cerr<<"pmin "<<p0<<" pmax "<<p0+np*dp<<endl;


  float time_eval1, time_eval2;
  float t, x1, x2, tau, p; 
  int l;

  ck0 = clock(); 
  tt0 = time(0);
  for (int i=0; i<ntau; i++)
    for (int j=0; j<np; j++) {
      udata[ntau*j+i] = 0.0;
      for (int m=0; m<nx1; m++) 
        for (int n=0; n<nx2; n++) {
          tau = tau0 + i*dtau;
          p = p0 + j*dp;
          x1 = x10 + m*dx1;
          x2 = x20 + n*dx2;
          t = sqrt(tau*tau+p*p*(x1*x1+x2*x2));
          l = int(round((t-t0)/dt));
          if (l>=0 && l<nt)
            udata[ntau*j+i] += fdata[nt*nx1*n+nt*m+l];   
        }
    }
  ck1 = clock();    
  tt1 = time(0);  
  time_eval1 = float(ck1-ck0)/CLOCKS_PER_SEC;
  time_eval2 = difftime(tt1,tt0);
  cerr<<"Teval1 "<<time_eval1<<endl;
  cerr<<"Teval2 "<<time_eval2<<endl;
  
  output << udata;

  exit(0);
}
