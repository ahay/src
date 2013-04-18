/* Array operations */

/*
  Copyright (C) 2006 The Board of Trustees of Stanford University
  
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

int i2(int i, int j, int n1)
/*< index-2 >*/
{
    int ij;
    ij=i*n1+j;
    return ij;
}

int i3(int i, int j, int k, int *n)
/*< index-3 >*/
{
    int ijk;
    ijk=i*n[0]+j*n[1]+k;
    return ijk;
}

void d3( int n2,int n3,int *n)
/*< 3-D >*/
{
    n[0]=n2*n3; 
    n[1]=n3;
}


int i4(int i, int j, int k, int l, int *n)
/*< index-4 >*/
{
    int ijkl;
    ijkl=i*n[0]+j*n[1]+k*n[2]+l;
    return ijkl;
    /*
      n[1]=n_j*n_k*n_l; n[2]=n_k*n_l; n[3]=n_l
    */
}

void d4(int n2,int n3,int n4,int*n)
/*< 4-D >*/
{
    n[0]=n2*n3*n4; 
    n[1]=n3*n4; 
    n[2]=n4;
}

int i5(int i, int j, int k, int l, int m, int *n){
    int ijklm;
    ijklm=i*n[0]+j*n[1]+k*n[2]+l*n[3]+m;
    return ijklm;
}

void d5(int n2,int n3,int n4,int n5,int *n){
    n[0]=n2*n3*n4*n5; n[1]=n3*n4*n5; n[2]=n4*n5; n[3]=n5;
}

int i6(int i, int j, int k, int l, int m, int p, int *n){
    int ijklmp;
    ijklmp=i*n[0]+j*n[1]+k*n[2]+l*n[3]+m*n[4]+p;
    return ijklmp;
}


void d6(int n2,int n3,int n4,int n5,int n6, int *n){
    n[0]=n2*n3*n4*n5*n6; n[1]=n3*n4*n5*n6; n[2]=n4*n5*n6; 
    n[3]=n5*n6; n[4]=n6;
}

void rowf( float *a,float *rowofa,int iy, int nx, int forward)
/*< float row >*/
{
    int ix;
    if (forward){
	for(ix=0;ix<nx; ix++)
	    rowofa[ix]=a[i2(iy,ix,nx)];
    }
    else{
	for(ix=0;ix<nx; ix++)
	    a[i2(iy,ix,nx)]=rowofa[ix];
    }
}

void rowc(sf_complex *a,sf_complex *rowofa,int iy, int nx,int forward)
/*< complex row >*/
{
    int ix;
    if (forward){
	for(ix=0;ix<nx; ix++)
	    rowofa[ix]=a[i2(iy,ix,nx)];
    }
    else{
	for(ix=0;ix<nx; ix++)
	    a[i2(iy,ix,nx)]=rowofa[ix];
    }
}

void colc(sf_complex *a,sf_complex *rowofa,int ix, int nx,int ny,int forward)
/*< complex column >*/
{
    int iy;
    if (forward){
	for(iy=0;iy<ny; iy++)
	    rowofa[iy]=a[i2(iy,ix,nx)];
    }
    else{
	for(iy=0;iy<ny; iy++)
	    a[i2(iy,ix,nx)]=rowofa[iy];
    }
}

void rowi( int *a,int *rowofa,int iy, int nx,int forward)
{
    int ix;
    if (forward){
	for(ix=0;ix<nx; ix++)
	    rowofa[ix]=a[i2(iy,ix,nx)];
    }
    else{
	for(ix=0;ix<nx; ix++)
	    a[i2(iy,ix,nx)]=rowofa[ix];
    }

}

float maxval( float *a, int nx)
/*< maximum value >*/
{
    float mxval;
    int ix;
    for(ix=0,mxval=a[0];ix<nx;ix++){
	if ( a[ix] > mxval) mxval=a[ix];
    }
    return mxval;
}

float maxvalc(sf_complex *a,int nx)
/*< maximum value >*/
{
    float mxval;
    int ix;
    for(ix=0,mxval=fabs(crealf(a[0]));ix<nx;ix++){
	if ( fabs(crealf(a[ix])) > mxval) mxval=fabs(crealf(a[ix]));
    }
    return mxval;
}

float minval( float *a,int nx)
/*< minimum value >*/
{
    float mnval;
    int ix;
    for(ix=0,mnval=a[0]; ix<nx; ix++){
	if (a[ix] <mnval) mnval=a[ix];
    }
    return mnval;
}

void matrix_equ_c(sf_complex *a,sf_complex *b,int ny,int nx)
/*< matrix >*/
{
    int iy,ix;
    for(iy=0;iy<ny;iy++)
	for(ix=0;ix<nx;ix++)
	    a[i2(iy,ix,nx)]=b[i2(iy,ix,nx)];
}

void vector_value_f(float *a, float value,int nx)
/*< fill a float array >*/
{
    int i;
    for (i=0;i<nx;i++)
	a[i]=value;
}

void vector_value_c(sf_complex *a,sf_complex value,int nx)
/*< fill a complex array >*/
{
    int i;
    for (i=0;i<nx;i++)
	a[i]=value;
}

void vector_cp_c(sf_complex *a, sf_complex *b, int nx)
/*< copy complex vector >*/
{
    int i;
    for(i=0;i<nx;i++) a[i]=b[i];
}

void vector_cp_f(float *a, float *b, int nx){
    int i;
    for(i=0;i<nx;i++) a[i]=b[i];
}
