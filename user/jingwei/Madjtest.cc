//   Ajoint test of prop1, prop2, prop3, prop4
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

    iRSF fft, model("model"), data("data"), left("left"), right("right"); // Get input
    fft.get("n1",nkz);
    fft.get("d1",dkz);
    fft.get("o1",kz0);
    fft.get("n2",nkx);
    fft.get("d2",dkx);
    fft.get("o2",kx0);
    left.get("n2",m2);

    model.get("n1",nz);
    model.get("d1",dz);
    model.get("o1",z0);
    model.get("n2",nx);
    model.get("d2",dx);
    model.get("o2",x0);

    nzx = nz*nx; 
    nkzx = nkz*nkx;

    std::valarray<sf_complex> mod(nzx), dat(nzx), lt(nzx*m2), rt(m2*nkzx);
    model >> mod;
    data >> dat;
    left >> lt;
    right >> rt;

    sf_complex *cleft, *cright, *cmod, *cdat, *cpmod, *cpdat;
    cleft = sf_complexalloc(nzx*m2);    
    cright = sf_complexalloc(m2*nkzx);
    cmod = sf_complexalloc(nzx);
    cdat = sf_complexalloc(nzx);   
    cpmod = sf_complexalloc(nzx);   
    cpdat = sf_complexalloc(nzx);  
    
    for (int i=0; i<nzx*m2; i++) cleft[i] = lt[i];
    for (int i=0; i<m2*nkzx; i++) cright[i] = rt[i];
    for (int i=0; i<nzx; i++) cmod[i] = mod[i];
    for (int i=0; i<nzx; i++) cdat[i] = dat[i];


    // Check inner product <cpmod,cdat>?=<cmod,cpdat>
    if (flag==1) {
	prop1( cmod, cpmod, cleft, cright, nz, nx, nkzx, m2, reg);  
	prop1( cdat, cpdat, cleft, cright, nz, nx, nkzx, m2, reg);  
    } else if (flag==2) {
	prop2( cmod, cpmod, cleft, cright, nz, nx, nkzx, m2, reg);  
	prop2( cdat, cpdat, cleft, cright, nz, nx, nkzx, m2, reg);  
    } else if (flag==3) {
	prop3( cmod, cpmod, cleft, cright, nz, nx, nkzx, m2, reg);  
	prop3( cdat, cpdat, cleft, cright, nz, nx, nkzx, m2, reg);  
    } else if (flag==4) {
	prop4( cmod, cpmod, cleft, cright, nz, nx, nkzx, m2, reg);  
	prop4( cdat, cpdat, cleft, cright, nz, nx, nkzx, m2, reg);  
    } else {
	cerr<<"Need to provide flag#"<<endl;
    }

    double aar=0., aai=0., bbr=0., bbi=0.;
    for (int i=0; i<nzx; i++) {
        aar += crealf(sf_cmul(cpmod[i],conjf(cdat[i])));
        aai += cimagf(sf_cmul(cpmod[i],conjf(cdat[i])));
        bbr += crealf(sf_cmul(cmod[i],conjf(cpdat[i])));
	bbi += cimagf(sf_cmul(cmod[i],conjf(cpdat[i])));
    }

    cerr<<"aa= "<<aar<<" "<<aai<<endl;
    cerr<<"bb= "<<bbr<<" "<<bbi<<endl;

    exit(0);
}
