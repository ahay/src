//   direct adjoint 2to2 hyper Radon transform (single integral, nearest point interpolation)
//   real f(tau,p) --> real u(t,x)
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

  int ntau, np;
  input.get("n1",ntau);
  input.get("n2",np);

  float tau0, dtau;
  input.get("o1",tau0);
  input.get("d1",dtau);  

  float p0, dp;
  input.get("o2",p0);
  input.get("d2",dp);

  std::valarray<float> fdata(ntau*np);

  input >> fdata;

  // Set output
  iRSF par(0);

  int nt, nx;
  par.get("nt",nt); 
  par.get("nx",nx); 

  float t0, dt;
  par.get("t0",t0);
  par.get("dt",dt);  

  float x0, dx;
  par.get("x0",x0);
  par.get("dx",dx);

  oRSF output;
  output.put("n1",nt);
  output.put("n2",nx);

  output.put("o1",t0);
  output.put("d1",dt);

  output.put("o2",x0);
  output.put("d2",dx);

  std::valarray<float> udata(nt*nx);

  cerr<<"taumin "<<tau0<<" taumax "<<tau0+ntau*dtau<<endl;
  cerr<<"pmin "<<p0<<" pmax "<<p0+np*dp<<endl;
  cerr<<"tmin "<<t0<<" tmax "<<t0+nt*dt<<endl;
  cerr<<"xmin "<<x0<<" xmax "<<x0+nx*dx<<endl;

  
  float time_eval1, time_eval2;
  float t, x, tau, p; 
  int m;

  ck0 = clock(); 
  tt0 = time(0);
  for (int i=0; i<nt; i++)
    for (int j=0; j<nx; j++) {
      udata[nt*j+i] = 0.0;
      for (int n=0; n<np; n++) {
        t = t0 + i*dt;
        x = x0 + j*dx;
        p = p0 + n*dp;
        tau = t*t-p*p*x*x;
        if (tau>=0) {
	  tau = sqrt(tau);
          m = int(round((tau-tau0)/dtau));  
          if (m>=0 && m<ntau) 
            udata[nt*j+i] += fdata[ntau*n+m];    
        }
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
