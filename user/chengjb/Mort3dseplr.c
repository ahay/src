/* Correct projection deviation for 3-D pseudo-pure P-wave field in ORT media.
   using low-rank symbol approximation

   Refernces:
             Cheng and Kang (SEG Abstract, 2012);
             Wang et al.(SEG Abstract, 2012)      

   Copyright (C) 2012 Tongji University, Shanghai, China 

   Authors: Jiubing Cheng
     
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
#include <assert.h>

/* low rank decomposition  */
#include "vecmatop.hh"
#include "serialize.hh"

/* prepared head files by myself */
#include "_cjb.h"

/* head files aumatically produced from *.c */
#include "puthead.h"
#include "zero.h"

/* spectra multiply for wave-mode separation*/
#include "kykxkztaper.h"
#include "spectramultipy.h"

int main(int  argc,char **argv)
{
    sf_init(argc,argv);

    sf_file Fi1, Fi2, Fi3;
    sf_file Fo1, Fo2;

    clock_t t1, t2;
    float   timespent;

    t1=clock();

   int isep;

   if (!sf_getint("isep",&isep)) isep=1;

    /* setup I/O files */
    Fi1 = sf_input("in");           /* pseudo-pure P-wave x-component */
    Fi2 = sf_input("PseudoPurePy"); /* pseudo-pure P-wave y-component */
    Fi3 = sf_input("PseudoPurePz"); /* pseudo-pure P-wave z-component */
    Fo1 = sf_output("out");               /* pseudo-pure scalar P-wave */
    Fo2 = sf_output("PseudoPureSepP");    /* separated scalar P-wave */

    int     nx, ny, nz, i, j, k;
    float   fx, fy, fz, dx, dy, dz;

    /* Read/Write axes */
    sf_axis az, ax, ay;

    az = sf_iaxa(Fi1,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
    ax = sf_iaxa(Fi1,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
    ay = sf_iaxa(Fi1,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
    fy=sf_o(ay)*1000.0;
    fx=sf_o(ax)*1000.0;
    fz=sf_o(az)*1000.0;

    int nxyz= nx*ny*nz, nxz=nx*nz;

    sf_warning("nx=%d ny=%d nz=%d dx=%f dy=%f dz=%f",nx,ny,nz,dx,dy,dz);

    puthead3x(Fo1, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
    puthead3x(Fo2, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);

    float* px=sf_floatalloc(nxyz);
    float* py=sf_floatalloc(nxyz);
    float* pz=sf_floatalloc(nxyz);

    float* p=sf_floatalloc(nxyz);

    int iy, ix, iz;

    i = 0;
    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
    {
      sf_floatread(&px[i],nz,Fi1);
      sf_floatread(&py[i],nz,Fi2);
      sf_floatread(&pz[i],nz,Fi3);
      i += nz;
    }
    
    for(i=0;i<nxyz;i++)
        p[i] = px[i] + py[i] + pz[i];

    i = 0;
    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
    {
      sf_floatwrite(&p[i],nz,Fo1);
      i += nz;
    }
    
    // separate qP-wave
    seplowrank3d(ldataxp,rdataxp,fmidxp,px,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2xp,n2xp,iflag);
    seplowrank3d(ldatayp,rdatayp,fmidyp,py,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2yp,n2yp,iflag);
    seplowrank3d(ldatazp,rdatazp,fmidzp,pz,ijkx,ijky,ijkz,nx,ny,nz,nxyz,nk,m2zp,n2zp,iflag);
   
    for(i=0;i<nxyz;i++)
        p[i] = px[i] + py[i] + pz[i];

    i = 0;
    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
    {
      sf_floatwrite(&p[i],nz,Fo2);
      i += nz;
    }

    free(**px);
    free(**py);
    free(**pz);
    free(**p);
		
    t2=clock();
    timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
    sf_warning("Computation time (Filtering): %f (second)",timespent);

    return 0;
}
