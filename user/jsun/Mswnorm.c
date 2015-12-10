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

float maxval(int n1, int n2, float *dat);
int scaling(float scale, int n1, int n2, float *dat);
int normalize(float den, int n1, int n2, float *dat);

int main(int argc, char * argv[])
{
    
    int dim,nz,nx,nt,nzx,nzxt,ndims[SF_MAX_DIM];
    int size,i;
    bool logsc;
    sf_file in, out;
    float *dat,den,scale;
    sf_axis az, ax, at;

    /* init RSF */
    sf_init (argc, argv);
    in = sf_input("in");
    out= sf_output("out");

    if (!sf_getint("size",&size)) size=1; /* sliding window size */
    if (!sf_getbool("log",&logsc)) logsc=true; /* log scaling */

    dim = sf_filedims(in,ndims);

    az = sf_iaxa(in, 1); nz = sf_n(az);
    ax = sf_iaxa(in, 2); nx = sf_n(ax);
    at = sf_iaxa(in, 3); nt = sf_n(at);

    nzx = nz*nx;
    nzxt = nzx*nt;
    
    sf_oaxa(out, az, 1);
    sf_oaxa(out, ax, 2);
    sf_oaxa(out, at, 3);

    dat = sf_floatalloc(nzxt);
    sf_floatread(dat, nzxt, in);

    for (i=0; i<size/2; i++) {
        den = maxval(nzx,size,dat);
        if (den <= SF_EPS) den = SF_EPS;
        if (logsc) scale = log10(den+1)/den;
        else scale = 1/den;
        scaling(scale,nzx,1,dat+i*nzx);
        sf_warning("i = %d;",i);
    }

    for (i=size/2; i<nt-size/2; i++) {
        den = maxval(nzx,size,dat+(i-size/2)*nzx);
        if (den <= SF_EPS) den = SF_EPS;
        if (logsc) scale = log10(den+1)/den;
        else scale = 1/den;
        scaling(scale,nzx,1,dat+i*nzx);
        sf_warning("i = %d;",i);
    }

    for (i=nt-size/2; i<nt; i++) {
        den = maxval(nzx,size,dat+(nt-size-1)*nzx);
        if (den <= SF_EPS) den = SF_EPS;
        if (logsc) scale = log10(den+1)/den;
        else scale = 1/den;
        scaling(scale,nzx,1,dat+i*nzx);
        sf_warning("i = %d;",i);
    }
    sf_warning(".");

    sf_floatwrite(dat, nzxt, out);

    exit(0);
    
}

float maxval(int n1, int n2, float *dat)
{
    float max = 0;
    int i,j;

    for (i=0; i<n2; i++)
        for (j=0; j<n1; j++)
            if (max<dat[i*n1+j]) max = dat[i*n1+j];

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
