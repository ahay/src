#include <rsf.h>

#include "mask6.h"

void mask32 (int nw, int nj1, int nj2, int nx, int ny, int nz, 
	     float ***yy, bool ***m1, bool ***m2)
{
    int ix, iy, iz, iw, is;
    bool ***xx;

    xx = sf_boolalloc3(nx,ny,nz);

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		xx[iz][iy][ix] = (yy[iz][iy][ix] == 0.);
		m1[iz][iy][ix] = false;
		m2[iz][iy][ix] = false;
	    }
	}
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny-1; iy++) {
	    for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj1;		  
		    m1[iz][iy][ix] = m1[iz][iy][ix] || 
			xx[iz][iy+1][ix+is] || xx[iz][iy][ix-is];
		}
	    }
	}
    }

    for (iz=0; iz < nz-1; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj2;		  
		    m2[iz][iy][ix] = m2[iz][iy][ix] ||
			xx[iz][iy+1][ix+is] || xx[iz][iy][ix-is];
		}
	    }
	}
    }
	
    free(xx[0][0]);
    free(xx[0]);
    free(xx); 
}

void mask3 (int nw, int nj, int nx, int ny, float **yy, bool **mm) 
{
    int ix, iy, iw, is;
    bool **xx;

    xx = sf_boolalloc2(nx,ny);
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    xx[iy][ix] = (yy[iy][ix] == 0.);
	    mm[iy][ix] = false;
	}
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj; ix < nx-nw*nj; ix++) {
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj;
		mm[iy][ix] = mm[iy][ix] || xx[iy+1][ix+is] || xx[iy][ix-is];
	    }
	}
    }
    
    free(xx[0]);
    free(xx);
}

void mask6 (int nw, int nj1, int nj2, int nx, int ny, float **yy, bool **mm) 
{
    int ix, iy, iw, is;
    bool **xx;

    xx = sf_boolalloc2(nx,ny);
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    mm[iy][ix] = (yy[iy][ix] == 0.);
	    xx[iy][ix] = false;
	}
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj1;
		xx[iy][ix] = xx[iy][ix] || mm[iy+1][ix+is] || mm[iy][ix-is];
	    }
	}
    }
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    mm[iy][ix] = false;
	}
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj2;
		mm[iy][ix] = mm[iy][ix] || xx[iy+1][ix+is] || xx[iy][ix-is];
	    }
	}
    }

    free(xx[0]);
    free(xx);
}
