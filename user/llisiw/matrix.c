/* Matrix inversion */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "matrix.h"

void MatrixInversion(float **A, int order, float **Y)
/*< matrix inversion main interface >*/
{
    /* calculate determinant */
    double det = 1.0/CalcDeterminant(A,order);

    /* allocate memory for co-factor matrix */
    float **minor = sf_floatalloc2(order-1,order-1);

    for(int j=0; j < order; j++) {
        for(int i=0; i < order; i++) {

            /* get the co-factor matrix of A(j,i) */
            GetMinor(A,minor,j,i,order);
            Y[i][j] = det*CalcDeterminant(minor,order-1);
            if((i+j)%2 == 1)
                Y[i][j] = -Y[i][j];
        }
    }
}

void GetMinor(float **src, float **dest, int row, int col, int order)
/*< get co-factor matrix >*/
{
    /* indicate which col and row is being copied to dest */
    int colCount=0, rowCount=0;

    for(int i=0; i < order; i++) {
        if(i != row) {

            colCount = 0;
            for(int j=0; j < order; j++) {
                if(j != col) {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }

            rowCount++;
        }
    }
}

double CalcDeterminant(float **mat, int order)
/*< recursive calculation of determinant >*/
{
    /* stop recursion when matrix is a single element */
    if(order == 1)
        return mat[0][0];

    /* the determinant value */
    float det = 0.;

    /* allocate memory for cofactor matrix */
    float **minor = sf_floatalloc2(order-1,order-1);

    for(int i=0; i < order; i++)
    {
        /* get minor of element (0,i) */
        GetMinor(mat,minor,0,i,order);
	
	/* recursion */
        det += (i%2==1?-1.0:1.0) * mat[0][i] * CalcDeterminant(minor,order-1);
    }

    return det;
}
