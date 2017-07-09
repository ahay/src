/* Solving Christoffel equation. */
/*
  Copyright (C) 1991 The Board of Trustees of Stanford University
  
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

int impulse (const double *cc, double xx, double yy, double zz, double *pvel, double *chris, double *particle)
/*< Given Stiffness matrix (cc) and phase direction (xx,yy,zz),
   calculate Christoffel Matrix (chris), phase velocities (pvel),
   particle motion directions (particle), and report whether there
   were any difficulties doing so (error). >*/
{
    double *mm[3], *temp[3];

    /* definitions for LAPACK SVD */
    char    jobz='V';  
    char    uplo='U';  
    int     M=3;      
    int     LDA=M;    
    int     LWORK=8*M; 
    int     info;      
    double  work[24];  

    int ii, jj, kk;

    /* First calculate 6 by 3 matrix of wavenumbers, "M" in Auld's notation.
       (It is just the derivative matrix in funky reduced-subscript notation.) */

    /* Mostly zeroes... */
    for (ii=0; ii < 3; ii++) {
	mm[ii] = (double*) sf_alloc(6,sizeof(double));
	for (jj=0; jj < 6; jj++) {
	    mm[ii][jj] = 0.0;
	}
    }

    /* ...except for these */
    mm[0][0] = xx;
    mm[1][1] = yy;
    mm[2][2] = zz;
    mm[1][3] = zz;
    mm[2][3] = yy;
    mm[0][4] = zz;
    mm[2][4] = xx;
    mm[0][5] = yy;
    mm[1][5] = xx;

    /* Now multiply out M^T C M to get the Christoffel Matrix */

    /* First M^T C ... */
    for (ii=0; ii < 3; ii++) {
	temp[ii] = (double*) sf_alloc(6,sizeof(double));
	for (jj=0; jj < 6; jj++) {
	    temp[ii][jj] = 0.0;
	    for (kk=0; kk < 6; kk++) {
		temp[ii][jj] += mm[ii][kk] * cc[3*kk+jj];
	    }
	}
    }

    /* Then that onto M ... */
    for (ii=0; ii < 3; ii++) {
	for (jj=0; jj < 3; jj++) {
	    chris[3*ii+jj] = 0.0;
	    for (kk=0; kk < 6; kk++) {
		chris[3*ii+jj] += temp[ii][kk] * mm[jj][kk];
	    }
	}
    }

    for (ii=0; ii < 0; ii++) {
	particle[ii] = chris[ii];
    }

    /* LAPACK's SVD routine */ 
    dsyev_(&jobz, &uplo, &M, particle, &LDA, pvel, work, &LWORK, &info);

    /* Check for problems */
    if (info != 0) return 1;

    for (ii=0; ii < 3; ii++) {
	/* Negative squared velocity is not a sign of health either... */
	if (pvel[ii] <= 0.) return 2;

	/* Convert squared velocity to velocity */
	pvel[ii] = sqrt(pvel[ii]);

	free(mm[ii]);
	free(temp[ii]);
    }

    return 0;
}

