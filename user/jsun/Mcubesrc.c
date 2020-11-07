/* Simple 2-D wave propagation with multi-threaded fftw3 */
/*
  Copyright (C) 2009 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    int it,iz,ix;
    int nt,nz,nx,nzx;

    float *ww,*rr,*ss;  /* I/O arrays*/
    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    sf_init(argc,argv);

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at);
    az = sf_iaxa(Fr,1); nz = sf_n(az); 
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); 
    nzx = nz*nx;

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3); 
    
    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt);  sf_floatread(ww,nt ,Fw);
    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);
    ss=sf_floatalloc(nzx*nt);

    /* MAIN LOOP */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
    for (it=0; it<nt; it++)
	for (ix = 0; ix < nx; ix++)
	    for (iz=0; iz < nz; iz++)
                ss[(it*nx+ix)*nz+iz] = ww[it]*rr[ix*nz+iz]; 

    /* write final wavefield to output */
    sf_floatwrite(ss,nzx*nt,Fo); 

    exit (0);
}
