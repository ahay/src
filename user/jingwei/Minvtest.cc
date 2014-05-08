//   Test inverse rank-1 approximation for lowrank wave propagation: prop1, prop2, prop3, prop4
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

extern "C" {
#include "prop1.h"
#include "prop2.h"
#include "prop3.h"
#include "prop4.h"
#include "cfft2nsps.h"
}

using namespace std;
using std::cerr;

static int nz, nx, nzx, nkz, nkx, nkzx, m2, flag;
static float dz, dx, z0, x0, dkz, dkx, kz0, kx0, reg;

//------------------------------------------------------------

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0); // Get parameters
    par.get("flag",flag);
    par.get("reg",reg);

    par.get("nz",nz); 
    par.get("dz",dz); 
    par.get("z0",z0); 
    par.get("nx",nx);    
    par.get("dx",dx);
    par.get("x0",x0);  

    iRSF input, fft("fft"), left("left"), right("right"), alpha("alpha"), beta("beta"); // Get input
    fft.get("n1",nkz);
    fft.get("d1",dkz);
    fft.get("o1",kz0);
    fft.get("n2",nkx);
    fft.get("d2",dkx);
    fft.get("o2",kx0);
    left.get("n2",m2);

    nzx = nz*nx; 
    nkzx = nkz*nkx;

    // Read input
    std::valarray<sf_complex> ipt(nzx), lt(nzx*m2), rt(m2*nkzx), al(nzx), be(nkzx);
    input >> ipt;
    left >> lt;
    right >> rt;
    alpha >> al;
    beta >> be;

    sf_complex *cipt, *clt, *crt, *cal, *cbe;
    cipt = sf_complexalloc(nzx);    
    clt = sf_complexalloc(nzx*m2);    
    crt = sf_complexalloc(m2*nkzx);
    cal = sf_complexalloc(nzx);
    cbe = sf_complexalloc(nkzx);
    for (int i=0; i<nzx; i++) cipt[i] = ipt[i];
    for (int i=0; i<nzx*m2; i++) clt[i] = lt[i];
    for (int i=0; i<m2*nkzx; i++) crt[i] = rt[i];
    for (int i=0; i<nzx; i++) cal[i] = al[i];
    for (int i=0; i<nkzx; i++) cbe[i] = be[i];


    // first apply prop to input to get output1
    sf_complex *copt1;
    copt1 = sf_complexalloc(nzx);
    if (flag==1) {
	prop1( cipt, copt1, clt, crt, nz, nx, nkzx, m2, reg);  
    } else if (flag==2) {
	prop2( cipt, copt1, clt, crt, nz, nx, nkzx, m2, reg);  
    } else if (flag==3) {
	prop3( cipt, copt1, clt, crt, nz, nx, nkzx, m2, reg);  
    } else if (flag==4) {
	prop4( cipt, copt1, clt, crt, nz, nx, nkzx, m2, reg);  
    } else {
	cerr<<"Need to provide flag#"<<endl;
    }

    // then apply approximate rank-1 inverse to output1 to get output2
    sf_complex *cdat1, *cdat2, *cfdat2, *cfdat3, *cdat3, *copt2;
    cdat1 = sf_complexalloc(nzx);
    cdat2 = sf_complexalloc(nkzx);
    cfdat2 = sf_complexalloc(nkzx);
    cfdat3 = sf_complexalloc(nkzx);
    cdat3 = sf_complexalloc(nkzx);
    copt2 = sf_complexalloc(nzx);    

    // divide by alpha
    for (int i=0; i<nzx; i++) cdat1[i] = sf_cdiv(copt1[i],cal[i]);    
   
    // forward fft
    int nz2, nx2, nk; 
    nk = cfft2_init(1,nz,nx,&nz2,&nx2);
    if(nk!=nkzx) cerr<<"nk discrepancy "<<endl;   
    for (int ix = 0; ix < nx2; ix++) {
	for (int iz = 0; iz < nz2; iz++) {
	    int i = iz+ix*nz;
	    int j = iz+ix*nz2;
	    if (ix<nx && iz<nz)
		cdat2[j] = cdat1[i];
	    else 
		cdat2[j] = sf_cmplx(0.,0.);
	}
    }    
    cfft2(cdat2,cfdat2);
    cfft2_finalize();  
    
    // divide by beta
    for (int i=0; i<nkzx; i++) cfdat3[i] = sf_cdiv(cfdat2[i],cbe[i]);
    
    // inverse fft
    nk = cfft2_init(1,nz,nx,&nz2,&nx2);
    if(nk!=nkzx) cerr<<"nk discrepancy "<<endl;    
    icfft2_allocate(cfdat3);    
    icfft2(cdat3,cfdat3);
    cfft2_finalize();
    for (int ix = 0; ix < nx; ix++) {
	for (int iz = 0; iz < nz; iz++) {
	    int i = iz+ix*nz;
	    int j = iz+ix*nz2;
	    copt2[i] = cdat3[j];
	}
    }

    // Write output1 and ouput2
    std::valarray<sf_complex> opt1(nzx), opt2(nzx);
    for (int i=0; i<nzx; i++) {
        opt1[i] = copt1[i];
        opt2[i] = copt2[i];
    }

    oRSF output1, output2("output2");
    output1.type(SF_COMPLEX);
    output1.put("n1",nz);
    output1.put("n2",nx);
    output1 << opt1;  

    output2.type(SF_COMPLEX);
    output2.put("n1",nz);
    output2.put("n2",nx);
    output2 << opt2;

    exit(0);
}
