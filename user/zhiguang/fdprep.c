/* Interface for Helmholtz discretization */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
/*^*/

#ifndef SuiteSparse_long
#ifdef UF_long
#define SuiteSparse_long UF_long
#endif
#endif
/*^*/

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
	    n = (pad1-2)*(pad2-2);
	    break;

	case 'j':
	    n = (pad1-2)*(pad2-2);
	    break;

	case 'c':
	    n = (pad1-2)*(pad2-2);
	    break;

	default:
	    sf_error("Fail to load discretization scheme.");
    }

    return n;
}

SuiteSparse_long fdprep_nz(const int pad1, const int pad2)
/*< number of triplets >*/
{
    SuiteSparse_long nz=0;
    int i, j;

    switch (order[0]) {
	case '5':
	    for (j=1; j < pad2-1; j++) {
		for (i=1; i < pad1-1; i++) {
		    if (i != 1) nz++;
		    if (i != pad1-2) nz++;
		    if (j != 1) nz++;
		    if (j != pad2-2) nz++;
		    nz++;
		}
	    }
	    break;

	case '9':
	    for (j=1; j < pad2-1; j++) {
		for (i=1; i < pad1-1; i++) {
		    if (i > 2) nz++;
		    if (i > 1) nz++;
		    if (i < pad1-2) nz++;
		    if (i < pad1-3) nz++;
		    if (j > 2) nz++;
		    if (j > 1) nz++;
		    if (j < pad2-2) nz++;
		    if (j < pad2-3) nz++;
		    nz++;
		}
	    }
	    break;

	case 'j':
	    for (j=1; j < pad2-1; j++) {
		for (i=1; i < pad1-1; i++) {
		    if (i != 1) nz++;
		    if (i != pad1-2) nz++;
		    if (j != 1) nz++;
		    if (j != pad2-2) nz++;
		    if (i != 1 && j != 1) nz++;
		    if (i != pad1-2 && j != pad2-2) nz++;
		    if (i != 1 && j != pad2-2) nz++;
		    if (i != pad1-2 && j != 1) nz++;
		    nz++;
		}
	    }
	    break;

	case 'c':
	    for (j=1; j < pad2-1; j++) {
		for (i=1; i < pad1-1; i++) {		    
		    if (i > 1 && j > 1) nz++;
		    if (i > 1) nz++;
		    if (i > 1 && j < pad2-2) nz++;
		    if (j < pad2-2) nz++;
		    if (i < pad1-2 && j < pad2-2) nz++;
		    if (i < pad1-2) nz++;
		    if (i < pad1-2 && j > 1) nz++;
		    if (j > 1) nz++;
		    if (i > 2 && j > 2) nz++;
		    if (i > 2 && j > 1) nz++;
		    if (i > 2) nz++;
		    if (i > 2 && j < pad2-2) nz++;
		    if (i > 2 && j < pad2-3) nz++;
		    if (i > 1 && j < pad2-3) nz++;
		    if (j < pad2-3) nz++;
		    if (i < pad1-2 && j < pad2-3) nz++;
		    if (i < pad1-3 && j < pad2-3) nz++;
		    if (i < pad1-3 && j < pad2-2) nz++;
		    if (i < pad1-3) nz++;
		    if (i < pad1-3 && j > 1) nz++;
		    if (i < pad1-3 && j > 2) nz++;
		    if (i < pad1-2 && j > 2) nz++;
		    if (j > 2) nz++;
		    if (i > 1 && j > 2) nz++;
		    nz++;
		}
	    }
	    break;

	default:
	    sf_error("Fail to load discretization scheme.");
    }

    return nz;
}

void fdprep(const double omega,
	    const int n1, const int n2,
	    const float d1, const float d2,
	    float **v,
	    const int npml,
	    const int pad1, const int pad2,
	    SuiteSparse_long *Ti, SuiteSparse_long *Tj,
	    double* Tx, double *Tz)
/*< discretization >*/
{
    switch (order[0]) {
	case '5':
	    fdprep5 (omega, n1,n2, d1,d2, v, npml,pad1,pad2, Ti,Tj,Tx,Tz);
	    break;

	case '9':
	    fdprep9 (omega, n1,n2, d1,d2, v, npml,pad1,pad2, Ti,Tj,Tx,Tz);
	    break;

	case 'j':
	    fdprep9o(omega, n1,n2, d1,d2, v, npml,pad1,pad2, Ti,Tj,Tx,Tz);
	    break;

	case 'c':
	    fdprep25(omega, n1,n2, d1,d2, v, npml,pad1,pad2, Ti,Tj,Tx,Tz);
	    break;

	default:
	    sf_error("Fail to load discretization scheme.");
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

	default:
	    sf_error("Fail to load discretization scheme.");
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

	default:
	    sf_error("Fail to load discretization scheme.");
    }
}
