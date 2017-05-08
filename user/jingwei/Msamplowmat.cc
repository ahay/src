//   Sample the lowrank propagation matrix 
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

using namespace std;
using std::cerr;

static int nz, nx, nzx, nkz, nkx, nkzx;
static float dz, dx, z0, x0, dkz, dkx, kz0, kx0, dt;

static std::valarray<float> zs, xs, kzs, kxs;   
static std::valarray<float> vels;

//------------------------------------------------------------

int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    
    for(int a=0; a<nr; a++) 
	for(int b=0; b<nc; b++) {
	    int iz = rs[a] % nz;
	    int ix = (int) rs[a]/nz;
            int ikz = cs[b] % nkz;
            int ikx = (int) cs[b]/nkz;
	    float hypk = hypot(kzs[ikz],kxs[ikx]);
            float phase = 2*SF_PI*vels[rs[a]]*hypk*dt;
	    res(a,b) = cpx(cos(phase),sin(phase));
	}
    return 0;
}
//------------------------------------------------------------

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0); // Get parameters
    par.get("dt",dt);

    par.get("nz",nz); 
    par.get("dz",dz); 
    par.get("z0",z0); 
    par.get("nx",nx);    
    par.get("dx",dx);
    par.get("x0",x0);  

    iRSF velocity, fft("fft"); // Get input
    fft.get("n1",nkz);
    fft.get("d1",dkz);
    fft.get("o1",kz0);
    fft.get("n2",nkx);
    fft.get("d2",dkx);
    fft.get("o2",kx0);

    nzx = nz*nx; 
    nkzx = nkz*nkx;

    // Read velocity information
    std::valarray<float> vel(nzx);
    velocity >> vel;
    vels.resize(nzx);
    vels = vel;

    // Provide axis information
    std::valarray<float> zz(nz), xx(nx), kz(nkz), kx(nkx);
    for(int i=0; i<nz; i++) zz[i] = z0+i*dz;
    for(int i=0; i<nx; i++) xx[i] = x0+i*dx;
    for(int i=0; i<nkz; i++) kz[i] = kz0+i*dkz;
    for(int i=0; i<nkx; i++) kx[i] = kx0+i*dkx;
    
    zs.resize(nz);
    xs.resize(nx);
    kzs.resize(nkz);
    kxs.resize(nkx);
    zs = zz;
    xs = xx;
    kzs = kz;
    kxs = kx;

    // Get the whole matrix
    vector<int> rs(nzx), cs(nkzx);
    for(int i=0; i<nzx; i++) rs[i] = i;
    for(int i=0; i<nkzx; i++) cs[i] = i;
    CpxNumMat res(nzx,nkzx);

    iC( sample(rs, cs, res) ); 
   
    std::valarray<sf_complex> cres(nzx*nkzx);
    for (int i=0; i<nzx; i++)
	for (int j=0; j<nkzx; j++) {
	    cres[i+j*nzx] = sf_cmplx(real(res(i,j)),imag(res(i,j)));
	}

    oRSF mat;
    mat.type(SF_COMPLEX);
    mat.put("n1",nzx);
    mat.put("n2",nkzx);
    mat << cres;  

    exit(0);
}
