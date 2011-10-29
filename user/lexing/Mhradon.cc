// 2-D hyperbolic Radon transform.
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

int main(int argc, char** argv)
{
  //0. init and get options
  srand48(time(NULL));
  clock_t ck0, ck1;
  time_t t0, t1;
  // map<string,string> opts;
  // optionsCreate(argc, argv, opts);

  sf_init(argc,argv); // Initialize RSF

  //1. data
  vector<int> all(1,1);
  map<string,string>::iterator mi;

  iRSF par(0);

  iRSF input;

  int nw, nx;
  input.get("n1",nw);
  input.get("n2",nx);
  if (nw != nx) sf_error("Need n1=n2"); // Need arbitrary nt and nx!

  std::valarray<sf_complex> fdata;
  fdata.resize(nx*nw);

  input >> fdata;

  CpxNumMat f(nw,nx);
  cpx *fd = f.data();
  for (int k=0; k < nw*nx; k++) {
      cpx ck(crealf(fdata[k]),cimagf(fdata[k]));
      fd[k] = ck;
  }
  
  CpxNumMat u(nw,nx);  setvalue(u,cpx(0,0));
  //2. bfio
  BFIO bfio("bfio_");
  iC( bfio.setup(par) );
  //
  double time_eval;
  if(nw<=256) {
    ck0 = clock();
    iC( bfio.eval(f,u) );
    ck1 = clock();    time_eval = double(ck1-ck0)/CLOCKS_PER_SEC;
  } else {
    t0 = time(0);
    iC( bfio.eval(f,u) );
    t1 = time(0);    time_eval = difftime(t1,t0);
  }
  //
  double relerr = 0;
  int NC = 64;
  ck0 = clock();
  iC( bfio.check(f,u,NC,relerr) );
  ck1 = clock();
  //  double time_chck = double(ck1-ck0)/CLOCKS_PER_SEC * nw * nw / double(NC);
  //
  // printf("RESULT\n");
  // printf("N  %d\n", N);
  // printf("Ta %.2e\n", time_eval);
  // printf("Td %.2e\n", time_chck);
  // printf("Rt %.2e\n", time_chck/time_eval);
  // printf("Ea %.2e\n", relerr);
  //

  cpx *ud = u.data();
  for (int k=0; k < nw*nx; k++) {
      fdata[k] = sf_cmplx(real(ud[k]),imag(ud[k]));
  }

  oRSF output;
  output << fdata;

  return 0;
}
