/* */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include <math.h>
#include <limits.h>

static float ***Gx, ***Gz;
static int *sxx, *sxz, *szx, *szz;
static int lenx, lenz;
static int marg;

void initsglfdcoe(int nxb, int nzb, int vlenx, int vlenz)
/*< allocation for SG LFD coefficients >*/
{
    lenx = vlenx;
    lenz = vlenz;

    Gx = sf_floatalloc3(nzb, nxb, lenx);
    Gz = sf_floatalloc3(nzb, nxb, lenz);
    
    sxx = sf_intalloc(lenx);
    sxz = sf_intalloc(lenx);
    szx = sf_intalloc(lenz);
    szz = sf_intalloc(lenz);
}

void loadcoe(int nzb, int nxb, sf_file FGx, sf_file FGz)
/*< read SG LFD coefficients >*/
{
    sf_floatread(Gx[0][0], nzb*nxb*lenx, FGx);
    sf_floatread(Gz[0][0], nzb*nxb*lenz, FGz);
}

void loadschm(sf_file Fsxx, sf_file Fsxz, sf_file Fszx, sf_file Fszz)
/*< read SG LFD scheme >*/
{
    float *sxxtmp, *sxztmp, *szxtmp, *szztmp;
    int ix, iz;
    int mx=0, mz=0;

    sxxtmp = sf_floatalloc(lenx);
    sxztmp = sf_floatalloc(lenx);
    szxtmp = sf_floatalloc(lenz);
    szztmp = sf_floatalloc(lenz);

    sf_floatread(sxxtmp, lenx, Fsxx);
    sf_floatread(sxztmp, lenx, Fsxz);
    sf_floatread(szxtmp, lenz, Fszx);
    sf_floatread(szztmp, lenz, Fszz);

    for (ix=0; ix<lenx; ix++) {
	sxx[ix] = (int)sxxtmp[ix];
	sxz[ix] = (int)sxztmp[ix];
	mx = abs(sxx[ix])>mx? abs(sxx[ix]):mx;
    }
    
    for (iz=0; iz<lenz; iz++) {
	szx[iz] = (int)szxtmp[iz];
	szz[iz] = (int)szztmp[iz];
	mz = abs(szz[iz])>mz? abs(szz[iz]):mz;
    }
    marg = mx>mz?mx:mz;
    free(sxxtmp); free(sxztmp); free(szxtmp); free(szztmp);
}

int getmarg(void)
/*< return marg (order of LFD) >*/
{
    return marg;
}

void freesglfdcoe(void)
/*< free allocation >*/
{
    free(**Gx);
    free(*Gx);
    free(Gx);
    free(**Gz);
    free(*Gz);
    free(Gz);
    free(sxx);
    free(sxz);
    free(szx);
    free(szz);
}

float ldx(float **data, int ix, int iz)
/*<Low rank finite difference : d/dx>*/
{
    float res = 0.0;
    int il;
    for (il = 0; il < lenx; il++) {
	res += 0.5*(-1*data[ix-sxx[il]+1][iz-sxz[il]] + data[ix+sxx[il]][iz+sxz[il]])*Gx[il][ix][iz];
    }
    return res;
}

float ldz(float **data, int ix, int iz)
/*<Low rank finite difference : d/dz>*/
{
    float res = 0.0;
    int il;
    for (il = 0; il < lenz; il++) {
	res += 0.5*(-1*data[ix-szx[il]][iz-szz[il]+1] + data[ix+szx[il]][iz+szz[il]])*Gz[il][ix][iz];
    }
    return res;
}
