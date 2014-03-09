//   Test sample function
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
#include "cfft2nsps.h"
}

using namespace std;
using std::cerr;

static int nz, nx, nzx, nkz, nkx, nkzx, m2;
static float dz, dx, z0, x0, dkz, dkx, kz0, kx0;

static sf_complex *cleft, *cright;
static std::valarray<float> zs, xs, kzs, kxs;   

//------------------------------------------------------------

// Sample rows and columns of the symbol
int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res, int hu)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    
    if (hu==1) {
	// sample columns
	for(int b=0; b<nc; b++) {
	    // for each fixed k0
	    int ikz = cs[b] % nkz;
	    int ikx = (int) cs[b]/nkz;

	    sf_complex *cin, *cout;
	    cin = sf_complexalloc(nzx);
	    cout = sf_complexalloc(nzx);

	    // construct cin = exp(2*pi*i*x*k0) for all x
	    float phs;
	    for(int i=0; i<nz; i++)
		for(int j=0; j<nx; j++) {
		    phs = 2*SF_PI*(zs[i]*kzs[ikz]+xs[j]*kxs[ikx]);
		    cin[j*nz+i] = sf_cmplx(cos(phs),sin(phs));
		}
	    
	    // apply cout = PP*cin
	    //prop1( cin, cout, cleft, cright, nz, nx, nkzx, m2);
            // apply cout = NN*cin
            prop2( cin, cout, cleft, cright, nz, nx, nkzx, m2);

	    for(int a=0; a<nr; a++) {
                int izx = rs[a];
                cout[izx] = sf_cmul(cout[izx], conjf(cin[izx]));
		res(a,b) = cpx(crealf(cout[izx]),cimagf(cout[izx]));
	    }
        }
    } else {
	// sample rows
	for(int a=0; a<nr; a++) {
	    // for each fixed x0
	    int iz = rs[a] % nz;
	    int ix = (int) rs[a]/nz;

	    sf_complex *cin, *cout;
	    cin = sf_complexalloc(nzx);
	    cout = sf_complexalloc(nzx);
            int nz2, nx2, nk;
            sf_complex *cfin, *cin1, *cout1, *cfout;
	    cfin = sf_complexalloc(nkzx);
	    cin1 = sf_complexalloc(nkzx);
            cout1 = sf_complexalloc(nkzx);
            cfout = sf_complexalloc(nkzx);
             
	    /*
            //--------construct cin = delta(x-x0) for all x------------
	    for(int i=0; i<nz; i++)
	    	for(int j=0; j<nx; j++) {
	    	    if(i==iz && j==ix)
	    		cin[j*nz+i] = sf_cmplx(1.0,0.0);
	    	    else 
	    		cin[j*nz+i] = sf_cmplx(0.0,0.0);
	    	}
            //---------------------------------------------------------
	    */
	    
            // construct cfin = exp(-2*pi*i*x0*k) for all k
            float phs;
	    for(int i=0; i<nkz; i++)
		for(int j=0; j<nkx; j++) {
		    phs = 2*SF_PI*(zs[iz]*kzs[i]+xs[ix]*kxs[j]);
		    cfin[j*nkz+i] = sf_cmplx(cos(phs),-sin(phs));
		}

            
            //--------construct cin in a different way----------------
            // inverse fft on cfin to get cin1
	    nk = cfft2_init(1,nz,nx,&nz2,&nx2);
	    if(nk!=nkzx) cerr<<"nk discrepancy "<<endl;
            icfft2_allocate(cfin);
	    icfft2(cin1,cfin);
	    cfft2_finalize();
            // truncate cin1 to get cin
	    for (int ix = 0; ix < nx; ix++) {
		for (int iz = 0; iz < nz; iz++) {
		    int i = iz+ix*nz;
		    int j = iz+ix*nz2;
		    cin[i] = cin1[j];
		}
	    }
	    //---------------------------------------------------------
	    

	    // apply cout = PP*cin
	    //prop1( cin, cout, cleft, cright, nz, nx, nkzx, m2);
            // apply cout = NN*cin
            prop2( cin, cout, cleft, cright, nz, nx, nkzx, m2);
       
            // pad cout to get cout1
            nk = cfft2_init(1,nz,nx,&nz2,&nx2);
	    if(nk!=nkzx) cerr<<"nk discrepancy "<<endl;
            for (int ix = 0; ix < nx2; ix++) {
		for (int iz = 0; iz < nz2; iz++) {
		    int i = iz+ix*nz;
		    int j = iz+ix*nz2;
		    if (ix<nx && iz<nz)
			cout1[j] = cout[i];
		    else 
			cout1[j] = sf_cmplx(0.,0.);
		}
	    }

	    // forward fft on cout1 to get cfout  
	    cfft2(cout1,cfout);
	    cfft2_finalize();  
           
	    for(int b=0; b<nc; b++) {
                int ikzx = cs[b];
		cfout[ikzx] = sf_cmul(conjf(cfout[ikzx]), cfin[ikzx]);
		res(a,b) = cpx(crealf(cfout[ikzx]),cimagf(cfout[ikzx]));
	    }
	}
    }
    return 0;
}
//------------------------------------------------------------

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0); // Get parameters

    int zx0, kzx0;
    par.get("zx0",zx0);
    par.get("kzx0",kzx0);

    par.get("nz",nz); 
    par.get("dz",dz); 
    par.get("z0",z0); 
    par.get("nx",nx);    
    par.get("dx",dx);
    par.get("x0",x0);  

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

    cleft = sf_complexalloc(nzx*m2);    
    cright = sf_complexalloc(m2*nkzx);
    for (int i=0; i<nzx*m2; i++) cleft[i] = lt[i];
    for (int i=0; i<m2*nkzx; i++) cright[i] = rt[i];

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

    // Get value of the symbol at (zx0, kzx0)
    vector<int> rs(1), cs(1);
    rs[0] = zx0;
    cs[0] = kzx0;

    CpxNumMat res(1,1);   
    iC( sample(rs, cs, res, 1) ); 
    cerr<<"from sample column: "<<res<<endl;

    iC( sample(rs, cs, res, 0) ); 
    cerr<<"from sample row: "<<res<<endl;

    exit(0);
}
