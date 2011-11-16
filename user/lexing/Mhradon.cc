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

    // Get input

    iRSF input;

    int nw, nx;
    input.get("n1",nw);
    input.get("n2",nx);

    std::valarray<sf_complex> fdata;
    fdata.resize(nx*nw);

    input >> fdata;

    CpxNumMat f(nw,nx);
    cpx *fd = f.data();
    for (int k=0; k < nw*nx; k++) {
	cpx ck(crealf(fdata[k]),cimagf(fdata[k]));
	fd[k] = ck;
    }
  
    // Set output
    int nt, np;
    par.get("nt",nt); // time samples
    par.get("np",np); // velocity samples

    CpxNumMat u(nt,np);  setvalue(u,cpx(0,0));


    oRSF output;
    output.put("n1",nt);
    output.put("n2",np);

    float o, d;

    par.get("t0",o);
    par.get("dt",d);
    output.put("o1",o);
    output.put("d1",d);

    par.get("p0",o);
    par.get("dp",d);
    output.put("o2",o);
    output.put("d2",d);

    output.type(SF_FLOAT);

    // BFIO setup
    BFIO bfio("bfio_");

    iC( bfio.setup(par,input) );
    //
    double time_eval;

    int n;
    par.get("n",n); // number of partitions

    if(n<=256) {
	ck0 = clock();
	iC( bfio.eval(n,f,u) );
	ck1 = clock();    time_eval = double(ck1-ck0)/CLOCKS_PER_SEC;
    } else {
	t0 = time(0);
	iC( bfio.eval(n,f,u) );
	t1 = time(0);    time_eval = difftime(t1,t0);
    }
    //
    double relerr = 0;
    int NC = 64;
    ck0 = clock();
    iC( bfio.check(n,f,u,NC,relerr) );
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

    std::valarray<float> udata;
    udata.resize(nt*np);

    for (int k=0; k < nt*np; k++) {
	udata[k] = real(ud[k]);
    }

    output << udata;

    exit(0);
}
