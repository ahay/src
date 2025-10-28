/* Matrix-vector multiplication using MPI*/
/*
  Copyright (C) 2025 University of Texas at Austin

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


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <rsf.h>
#include <string.h>

int main(int argc, char* argv[])
{
    int rank, size;
    int nrows=0, ncols=0;        // global matrix dimensions
    int myrows=0, offset=0;      // local rows (or local columns when adj)
    bool adj = false;
    sf_file in=NULL, out=NULL, mat=NULL;
    float *x=NULL, *y=NULL;
    float *alocal=NULL, *ylocal=NULL;

    MPI_Init(&argc, &argv);
    sf_init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Root reads matrix and vector info */
    if(rank==0){
        in  = sf_input("--input");   // vector x
        out = sf_output("--output"); // result y
        mat = sf_input("mat");       // matrix A

        if(!sf_histint(in,"n1",&ncols)) sf_error("No n1= in input vector");
        if(!sf_getbool("adj",&adj)) adj=false;

        if(!sf_histint(mat,"n1",&ncols)) sf_error("No n1= in matrix");
        if(!sf_histint(mat,"n2",&nrows)) sf_error("No n2= in matrix");

        sf_putint(out,"n1", adj ? ncols : nrows);
    }

    /* Broadcast global sizes and adj flag */
    MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&adj,   1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    /* Decide local rows/columns for this rank */
    int rows_per_rank, remainder;
    if(adj){
        rows_per_rank = ncols / size;  // distributing columns when adj
        remainder     = ncols % size;
    } else {
        rows_per_rank = nrows / size;  // distributing rows when forward
        remainder     = nrows % size;
    }

    myrows = rows_per_rank + (rank < remainder ? 1 : 0);
    offset = rank*rows_per_rank + (rank < remainder ? rank : remainder);

    /* Allocate buffers */
    if(adj){
        /* alocal stores myrows columns, each column has length nrows:
           contiguous layout: alocal[i*nrows + j] = A[j, offset+i] */
        alocal = sf_floatalloc((size_t)myrows * nrows);
        ylocal = sf_floatalloc(ncols);          // full y (will be reduced)
    } else {
        alocal = sf_floatalloc((size_t)myrows * ncols); // local rows
        ylocal = sf_floatalloc(myrows);                 // local result
    }

    /* IMPORTANT: allocate x with correct length depending on adj */
    int xlen = adj ? nrows : ncols;
    x = sf_floatalloc(xlen);
    if(!alocal || !ylocal || !x) MPI_Abort(MPI_COMM_WORLD,1);

    /* Root reads vector with correct length, then broadcast that many floats */
    if(rank==0){
        if(adj)
            sf_floatread(x, nrows, in);   /* for adjoint x length must be nrows */
        else
            sf_floatread(x, ncols, in);   /* for forward x length must be ncols */
    }
    MPI_Bcast(x, xlen, MPI_FLOAT, 0, MPI_COMM_WORLD);

    /* Root reads full matrix and scatters (keeps your column-send approach for adj) */
    if(rank==0){
        float *Aall = sf_floatalloc((size_t)nrows * ncols);
        sf_floatread(Aall, nrows*ncols, mat);

        for(int r=0;r<size;r++){
            int rrows = rows_per_rank + (r < remainder ? 1 : 0);
            int roffs = r*rows_per_rank + (r < remainder ? r : remainder);

            if(adj){
                /* send rrows columns to rank r; each column has nrows floats */
                for(int i=0;i<rrows;i++){
                    if(r==0){
                        for(int j=0;j<nrows;j++)
                            alocal[i*(size_t)nrows + j] = Aall[j*(size_t)ncols + roffs + i];
                    } else {
                        float *tmp = malloc((size_t)nrows * sizeof(float));
                        for(int j=0;j<nrows;j++)
                            tmp[j] = Aall[j*(size_t)ncols + roffs + i];
                        MPI_Send(tmp, nrows, MPI_FLOAT, r, 0, MPI_COMM_WORLD);
                        free(tmp);
                    }
                }
            } else {
                if(r==0){
                    memcpy(alocal, Aall + roffs*(size_t)ncols, (size_t)rrows*(size_t)ncols*sizeof(float));
                } else {
                    MPI_Send(Aall + roffs*(size_t)ncols, rrows*ncols, MPI_FLOAT, r, 0, MPI_COMM_WORLD);
                }
            }
        }
        free(Aall);
    } else {
        if(adj){
            for(int i=0;i<myrows;i++){
                MPI_Recv(alocal + i*(size_t)nrows, nrows, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            MPI_Recv(alocal, myrows*ncols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    /* Local mat-vec multiply */
    if(adj){
        /* y = A^T x
           alocal is stored per-local-column: alocal[i*nrows + j] == A[j, offset+i]
           We'll compute local contribution to full y (length ncols),
           so initialize ylocal to zeros and add into ylocal[offset+i]. */
        for(int j=0;j<ncols;j++) ylocal[j] = 0.0f;   /* zero full-length buffer */

        for(int i=0;i<myrows;i++){                   /* local columns (0..myrows-1) */
            int col = offset + i;                    /* global column index */
            float accum = 0.0f;
            for(int j=0;j<nrows;j++){
                /* alocal[i*nrows + j] == A[j, col] */
                /* x[j] corresponds to row j of A (x length is nrows for adjoint) */
                accum += alocal[i*(size_t)nrows + j] * x[j];
            }
            ylocal[col] = accum; /* write contribution into the correct position */
        }

        /* reduce contributions across ranks: each rank has nonzero ylocal only
           on indices [offset .. offset+myrows-1], other entries zero */
        y = sf_floatalloc(ncols);
        MPI_Reduce(ylocal, y, ncols, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        /* y = A x (forward). alocal stores rows: alocal[i*ncols + j] = A[offset+i, j] */
        for(int i=0;i<myrows;i++){
            float sum = 0.0f;
            for(int j=0;j<ncols;j++){
                sum += alocal[i*(size_t)ncols + j] * x[j];
            }
            ylocal[i] = sum;
        }
    }

    /* Gather / write results */
    if(!adj){
        if(rank==0){
            y = sf_floatalloc(nrows);
            int pos = 0;
            for(int r=0;r<size;r++){
                int rrows = rows_per_rank + (r < remainder ? 1 : 0);
                if(r==0){
                    memcpy(y, ylocal, rrows * sizeof(float));
                } else {
                    MPI_Recv(y + pos, rrows, MPI_FLOAT, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                pos += rrows;
            }
            sf_floatwrite(y, nrows, out);
            free(y);
        } else {
            MPI_Send(ylocal, myrows, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    } else {
        if(rank==0){
            sf_floatwrite(y, ncols, out);
            free(y);
        }
    }

    free(alocal);
    free(ylocal);
    free(x);

    MPI_Finalize();
    return 0;
}
