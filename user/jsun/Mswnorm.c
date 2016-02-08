/* Sliding window normalization */

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

float maxval(int n1, int n2, float *dat);
int scaling(float scale, int n1, int n2, float *dat);
int normalize(float den, int n1, int n2, float *dat);

int main(int argc, char * argv[])
{
    
    int dim,nz,nx,nt,nzx,nzxt,ndims[SF_MAX_DIM];
    int size,i,nth;
    bool sw,logsc;
    sf_file in, out;
    float *dat0,*dat,den,scale,rescale;
    float max_all,perc,thres;
    sf_axis az, ax, at;

    /* init RSF */
    sf_init (argc, argv);
    in = sf_input("in");
    out= sf_output("out");

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
    sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    if (!sf_getint("size",&size)) size=0; /* sliding window radius */
    if (!sf_getbool("sw",&sw)) sw=true; /* sliding window */
    if (!sf_getbool("log",&logsc)) logsc=false; /* log scaling */
    if (!sf_getfloat("perc",&perc)) perc=5; /* threshold percentage of the maximum value */
    perc /= 100.;

    dim = sf_filedims(in,ndims);

    az = sf_iaxa(in, 1); nz = sf_n(az);
    ax = sf_iaxa(in, 2); nx = sf_n(ax);
    at = sf_iaxa(in, 3); nt = sf_n(at);

    nzx = nz*nx;
    nzxt = nzx*nt;
    
    sf_oaxa(out, az, 1);
    sf_oaxa(out, ax, 2);
    sf_oaxa(out, at, 3);

    dat0 = sf_floatalloc(nzxt);
    dat  = sf_floatalloc(nzxt);
    sf_floatread(dat0, nzxt, in);

    max_all = maxval(nzx,nt,dat0);
    thres = max_all*perc;
    
    /* for log scaling (set maximum to 1 for better mapping) */
    if (logsc) {
        scaling(1./thres,nzx,nt,dat0);
        max_all *= 1./thres;
        thres = 1;
        rescale = log10(max_all+1);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<nzxt; i++) dat[i] = dat0[i];

    if (!sw) {

        scaling(1./max_all,nzx,nt,dat);

    } else {

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,scale)
#endif
    for (i=0; i<size; i++) {
        sf_warning("i = %d/%d;",i,nt);
        den = maxval(nzx,i+1+size,dat0);
        if (logsc) scale = log10(den+1.)/rescale;
        else scale = 1.;
        if (den <= thres) den = thres;
        scale /= den;
        scaling(scale,nzx,1,dat+i*nzx);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,scale)
#endif
    for (i=size; i<nt-size; i++) {
        sf_warning("i = %d/%d;",i,nt);
        den = maxval(nzx,2*size+1,dat0+(i-size)*nzx);
        if (logsc) scale = log10(den+1.)/rescale;
        else scale = 1.;
        if (den <= thres) den = thres;
        scale /= den;
        scaling(scale,nzx,1,dat+i*nzx);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,scale)
#endif
    for (i=nt-size; i<nt; i++) {
        sf_warning("i = %d/%d;",i,nt);
        den = maxval(nzx,size+1+(nt-1-i),dat0+(i-size)*nzx);
        if (logsc) scale = log10(den+1.)/rescale;
        else scale = 1.;
        if (den <= thres) den = thres;
        scale /= den;
        scaling(scale,nzx,1,dat+i*nzx);
    }
    sf_warning(".");

    }

    sf_floatwrite(dat, nzxt, out);

    exit(0);
    
}

/* maximum absolute value */
float maxval(int n1, int n2, float *dat)
{
    float max = 0;
    int i,j;

    for (i=0; i<n2; i++)
        for (j=0; j<n1; j++)
        { 
            if (max<fabs(dat[i*n1+j])) max = fabs(dat[i*n1+j]);
            //sf_warning("dat[%d]=%g,max=%g",i*n1+j,dat[i*n1+j],max); 
        }

    return max;
}

int scaling(float scale, int n1, int n2, float *dat)
{
    int i,j;

    for (i=0; i<n2; i++)
        for (j=0; j<n1; j++)
            dat[i*n1+j] *= scale;

    return 0;
}

int normalize(float den, int n1, int n2, float *dat)
{
    int i,j;

    for (i=0; i<n2; i++)
        for (j=0; j<n1; j++)
            dat[i*n1+j] /= den;

    return 0;
}
