//   Test rank-1 approximation for lowrank decomposition matrices
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
#include "rankoneapprox.hh"

using namespace std;
using std::cerr;

static int nz, nx, nzx, nkz, nkx, nkzx, m2;
static float dz, dx, z0, x0, dkz, dkx, kz0, kx0;

//------------------------------------------------------------

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0); // Get parameters
    par.get("nz",nz); 
    par.get("dz",dz); 
    par.get("z0",z0); 
    par.get("nx",nx);    
    par.get("dx",dx);
    par.get("x0",x0);  

    int npk;
    par.get("npk",npk); // maximum sample rows/columns
    
    iRSF fft, left("left"), right("right"); // Get input
    fft.get("n1",nkz);
    fft.get("d1",dkz);
    fft.get("o1",kz0);
    fft.get("n2",nkx);
    fft.get("d2",dkx);
    fft.get("o2",kx0);
    left.get("n2",m2);

    nzx = nz*nx; 
    nkzx = nkz*nkx;

    // Read lowrank matrices
    std::valarray<sf_complex> lt(nzx*m2), rt(m2*nkzx);
    left >> lt;
    right >> rt;

    CpxNumMat cleft(nzx,m2), cright(m2,nkzx);
    for (int i=0; i<nzx; i++) 
	for (int j=0; j<m2; j++)
	    cleft(i,j) = cpx(crealf(lt[j*nzx+i]),cimagf(lt[j*nzx+i]));
    for (int i=0; i<m2; i++) 
	for (int j=0; j<nkzx; j++)
	    cright(i,j) = cpx(crealf(rt[j*m2+i]),cimagf(rt[j*m2+i]));
    CpxNumVec aa, bb;
    iC( rankoneapprox(cleft, cright, aa, bb, npk) );
    
    // Write alpha and beta
    std::valarray<sf_complex> adata(nzx), bdata(nkzx);
    for (int i=0; i<nzx; i++)
        adata[i] = sf_cmplx(real(aa(i)),imag(aa(i)));
    for (int i=0; i<nkzx; i++) 
        bdata[i] = sf_cmplx(real(bb(i)),imag(bb(i)));

    oRSF alpha, beta("beta");
    alpha.type(SF_COMPLEX);
    alpha.put("n1",nzx);
    alpha.put("n2",1);
    alpha << adata;  

    beta.type(SF_COMPLEX);
    beta.put("n1",nkzx);
    beta.put("n2",1);
    beta << bdata;

    exit(0);
}
