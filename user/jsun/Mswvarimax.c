/* Sliding window varimax */

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

float varimax(int n1, int n2, float thres, float *dat);

int main(int argc, char * argv[])
{
    
    int fz,fx,nz0,nz,nx0,nx,nt0,nt,nzx0,nzx,nzxt0,nzxt;
    int it,iz,ix,nth,size;
    float thres,term;
    bool sw;
    sf_file in, out;
    float *dat0,*dat,*vari;
    sf_axis az, ax, at;

    /* init RSF */
    sf_init (argc, argv);
    in = sf_input("in");
    out= sf_output("out");

    if (!sf_getbool("sw",&sw)) sw=false; /* sliding window */
    if (!sf_getint("size",&size)) size=0; /* sliding window radius */
    if (!sf_getfloat("thres",&thres)) thres=0.; /* variance threshold (normalized) */
    if (!sf_getfloat("term",&term)) term=100.; /* variance threshold (normalized) */
    term /= 100.;

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
    sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    az = sf_iaxa(in, 1); nz0 = sf_n(az);
    ax = sf_iaxa(in, 2); nx0 = sf_n(ax);
    at = sf_iaxa(in, 3); nt0 = sf_n(at);
    nt = (int) (nt0*term);

    if (!sf_getint("f1",&fz)) fz=0; /* sliding window radius */
    if (!sf_getint("f2",&fx)) fx=0; /* sliding window radius */
    if (!sf_getint("n1",&nz)) nz=nz0-fz; /* sliding window radius */
    if (!sf_getint("n2",&nx)) nx=nx0-fx; /* sliding window radius */

    nzx0 = nz0*nx0;
    nzx  = nz *nx ;
    nzxt0= nzx0*nt0;
    nzxt = nzx *nt ;
    
    sf_oaxa(out, at, 1);
    sf_putint(out,"n2",1);
    sf_putint(out,"n3",1);

    dat0 = sf_floatalloc(nzxt0);
    dat  = sf_floatalloc(nzxt);
    vari = sf_floatalloc(nt0);
    sf_floatread(dat0, nzxt0, in);

    for         (it=0; it<nt; it++) {
        for     (ix=0; ix<nx; ix++) {
            for (iz=0; iz<nz; iz++) {
                dat[(it*nx+ix)*nz+iz] = dat0[(it*nx0+(fx+ix))*nz0+fz+iz];
            }
        }
    }

    if(!sw) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it)
#endif
        for (it=0; it<nt; it++) {
            sf_warning("it = %d/%d;",it,nt);
            vari[it] = varimax(nzx,1,thres,dat+it*nzx);
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel default(shared) 
#endif
        {
#ifdef _OPENMP
#pragma omp for private(it)
#endif
            for (it=0; it<size; it++) {
                sf_warning("it = %d/%d;",it,nt);
                vari[it] = varimax(nzx,it+1+size,thres,dat);
            }

#ifdef _OPENMP
#pragma omp for private(it)
#endif
            for (it=size; it<nt-size; it++) {
                sf_warning("it = %d/%d;",it,nt);
                vari[it] = varimax(nzx,2*size+1,thres,dat+(it-size)*nzx);
            }

#ifdef _OPENMP
#pragma omp for private(it)
#endif
            for (it=nt-size; it<nt; it++) {
                sf_warning("it = %d/%d;",it,nt);
                vari[it] = varimax(nzx,size+1+(nt-1-it),thres,dat+(it-size)*nzx);
            }
        }

    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it)
#endif
    for (it=nt; it<nt0; it++) {
        vari[it] = 0.;
    }

    sf_floatwrite(vari, nt0, out);

    exit(0);
}

/* maximum absolute value and variance (optional) */
float varimax(int n1, int n2, float thres, float *dat)
{
    float dd, vari, sum = 0, sum2 = 0;
    int i,j;

    for (i=0; i<n2; i++)
        for (j=0; j<n1; j++)
        { 
            dd = powf(dat[i*n1+j],2);
            sum += dd;
            sum2 += dd*dd;
        }

    if (sum>thres) vari = n1*n2*sum2/sum/sum;
    else       vari = 0;

    return vari;
}
