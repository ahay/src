/* Interface for Helmholtz discretization */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
#include <umfpack.h>

#include "fdprep5.h"
#include "fdprep9.h"
#include "fdprep9o.h"
#include "fdprep25.h"

#include "fdprep.h"

static char *order;

void fdprep_order(char *order0)
/*< set order >*/
{
    order = order0;
}

SuiteSparse_long fdprep_n(const int pad1, const int pad2)
/*< domain size >*/
{
    SuiteSparse_long n=0;

    switch (order[0]) {
	case '5':
	    n = (pad1-2)*(pad2-2);
	    break;

	case '9':
	    n = (pad1-4)*(pad2-4);
	    break;

	case 'j':
	    n = (pad1-2)*(pad2-2);
	    break;

	case 'c':
	    n = (pad1-4)*(pad2-4);
	    break;
    }

    return n;
}

SuiteSparse_long fdprep_nz(const int pad1, const int pad2)
/*< number of triplets >*/
{
    SuiteSparse_long nz=0;

    switch (order[0]) {
	case '5':
	    nz = 5*(pad1-2)*(pad2-2)
		-2*(pad1-4)-2*(pad2-4)-8;
	    break;

	case '9':
	    nz = 9*(pad1-4)*(pad2-4)
		-4*(pad1-6)-4*(pad2-6)-16
		-2*(pad1-8)-2*(pad2-8)-8;
	    break;

	case 'j':
	    nz = 9*(pad1-2)*(pad2-2)
		-6*(pad1-4)-6*(pad2-4)-20;
	    break;

	case 'c':
	    nz = 25*(pad1-4)*(pad2-4)
		-20*(pad1-8)-20*(pad2-8)-52-64
		-10*(pad1-8)-10*(pad2-8)-36;
	    break;
    }

    return nz;
}

void fdprep(const double omega,
	    const float a0, const float f0,
	    const int n1, const int n2,
	    const float d1, const float d2,
	    float **v,
	    const int npml,
	    const int pad1, const int pad2,
	    SuiteSparse_long n, SuiteSparse_long nz,
	    SuiteSparse_long *Ti, SuiteSparse_long *Tj,
	    double* Tx, double *Tz)
/*< discretization >*/
{
    switch (order[0]) {
	case '5':
	    fdprep5 (omega, a0,f0, n1,n2, d1,d2, v, npml,pad1,pad2, n,nz,Ti,Tj,Tx,Tz);
	    break;

	case '9':
	    fdprep9 (omega, a0,f0, n1,n2, d1,d2, v, npml,pad1,pad2, n,nz,Ti,Tj,Tx,Tz);
	    break;

	case 'j':
	    fdprep9o(omega, a0,f0, n1,n2, d1,d2, v, npml,pad1,pad2, n,nz,Ti,Tj,Tx,Tz);
	    break;

	case 'c':
	    fdprep25(omega, a0,f0, n1,n2, d1,d2, v, npml,pad1,pad2, n,nz,Ti,Tj,Tx,Tz);
	    break;
    }
}

void fdpad(const int npml,
	   const int pad1, const int pad2,
	   sf_complex **dat,
	   double *Bx, double *Bz)
/*< pad >*/
{
    switch (order[0]) {
	case '5':
	    fdpad5 (npml,pad1,pad2, dat,Bx,Bz);
	    break;

	case '9':
	    fdpad9 (npml,pad1,pad2, dat,Bx,Bz);
	    break;

	case 'j':
	    fdpad9o(npml,pad1,pad2, dat,Bx,Bz);
	    break;

	case 'c':
	    fdpad25(npml,pad1,pad2, dat,Bx,Bz);
	    break;
    }
}

void fdcut(const int npml,
	   const int pad1, const int pad2,
	   sf_complex **dat,
	   double *Xx, double *Xz)
/*< cut >*/
{
    switch (order[0]) {
	case '5':
	    fdcut5 (npml,pad1,pad2, dat,Xx,Xz);
	    break;

	case '9':
	    fdcut9 (npml,pad1,pad2, dat,Xx,Xz);
	    break;

	case 'j':
	    fdcut9o(npml,pad1,pad2, dat,Xx,Xz);
	    break;

	case 'c':
	    fdcut25(npml,pad1,pad2, dat,Xx,Xz);
	    break;
    }
}
