/* Complex-valued compensation (between two wavefields) */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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
#include "compensate.h"

int main(int argc, char* argv[])
{
    bool verb,cmplx;
    int i; /* index variables */
    int nz,nx,nzx;
    float perc,vmax,eps;

    /* I/O arrays*/
    sf_complex *num,*num2,*den,*wht;
    /* I/O files */
    sf_file Fnum,Fden,Fres;
    /* cube axes */
    sf_axis at,az,ax;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(!sf_getbool("cmplx",&cmplx)) cmplx=true; /* use complex i/o */
    if(!sf_getfloat("perc",&perc)) perc=0.1; /* precentage (of max) for protection when dividing */
    perc /= 100.;

    /* setup I/O files */
    Fnum = sf_input("in" );
    Fden = sf_input("den");
    Fres = sf_output("out");

    if (SF_COMPLEX != sf_gettype(Fnum)) sf_error("Need complex input");
    if (SF_COMPLEX != sf_gettype(Fden)) sf_error("Need complex den");

    sf_settype(Fres,SF_COMPLEX);

    /* Read/Write axes */
    az = sf_iaxa(Fnum,1); nz = sf_n(az); 
    ax = sf_iaxa(Fnum,2); nx = sf_n(ax); 
    at = sf_iaxa(Fnum,3);

    sf_oaxa(Fres,az,1); 
    sf_oaxa(Fres,ax,2); 
    sf_oaxa(Fres,at,3);
    
#ifdef _OPENMP
    int nth;
#pragma omp parallel
    {
        nth = omp_get_num_threads();
    }
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    nzx = nz*nx;

    /* Input wavefields */
    num  = sf_complexalloc(nzx);
    num2 = sf_complexalloc(nzx);
    den  = sf_complexalloc(nzx);
    wht  = sf_complexalloc(nzx);

    sf_complexread(num,nzx,Fnum);
    sf_complexread(den,nzx,Fden);

    sf_fileclose(Fnum);
    sf_fileclose(Fden);

    /***************************************/
    /* start the work */

    /* stable division in the space domain: num/den */
    vmax = find_cmax(nzx, den);
    eps = (vmax==0) ? SF_EPS : vmax*vmax*perc;
    if (verb) sf_warning("Space domain: vmax=%f, eps=%f",vmax,eps);
    stable_cdiv(nzx, eps, num, den, wht);

    /* apply scaling */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < nzx; i++) {
#ifdef SF_HAS_COMPLEX_H
        num2[i] = num[i]*wht[i];
#else
        num2[i] = sf_cmul(num[i],wht[i]);
#endif
    }

    /* output result */
    sf_complexwrite(num2,nzx,Fres);

    exit (0);
}
