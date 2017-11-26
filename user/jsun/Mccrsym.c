/* determine symmetry using correlation */

/*
  Copyright (C) 2015 University of Texas at Austin
  
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
#include <math.h>
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

float ccr(int nt, int size, int offset, float pad, float *dat);

int main(int argc, char * argv[])
{
    
    int nt, i,nth,size;
    sf_file in, out;
    float *dat0,*dat,pad;
    sf_axis at;

    /* init RSF */
    sf_init (argc, argv);
    in = sf_input("in");
    out= sf_output("out");

    if (!sf_getint("size",&size)) size=0; /* sliding window radius */
    if (!sf_getfloat("pad",&pad)) pad=SF_EPS; /* pad for stable devision */

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
    sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    at = sf_iaxa(in, 1); nt = sf_n(at);
    sf_oaxa(out, at, 1);

    dat0 = sf_floatalloc(nt);
    dat  = sf_floatalloc(nt);
    sf_floatread(dat0, nt, in);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<nt; i++) {
        dat[i] = ccr(nt, size, i, pad, dat0);
    }

    sf_floatwrite(dat, nt, out);

    exit(0);
}

/* maximum absolute value and variance (optional) */
float ccr(int nt, int size, int offset, float pad, float *dat)
{
    float num=0., den1=0., den2=0.,d1,d2,res;
    int i;

    if (size>offset) size = offset;
    if (offset+size>=nt) size = nt-1-offset;

    for (i=1; i<=size; i++) { 
        d1 = dat[offset-i];
        d2 = dat[offset+i];
        num += d1*d2;
        den1 += d1*d1;
        den2 += d2*d2;
    }

    if (den1==0. || den2 ==0.) res=0.;
    else {
        res = num/(sqrtf(den1)*sqrtf(den2)+pad);
    }

    return res;
}
