/* MPI example: compute PI using MPI, study the basic use of collections
   NB: PI=\int_{x} 4/(1+x*x) dx
   Compiling: mpicc -o mpipi mpipi.c; mpirun -np 4 ./mpipi; enter 400000;
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <math.h>
#include <mpi.h>

float f(float x)
{
  return 4.0/(1.0+x*x);
}

int main(int argc, char* argv[])
{
  int n, i;
  float t0, t1, h, sum, mypi, pi, x;
  int rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0){
    printf("computing pi using %d processes!\n",size);
    printf("enter the number of intervals:");
    scanf("%d",&n);
  }
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  t0=MPI_Wtime();
  h=1.0/(float)n;
  sum=0.0;
  for(i=rank; i<n; i+=size){
    x=h*((float)i+0.5);
    sum+=f(x);
  }
  mypi=h*sum;

  MPI_Reduce(&mypi, &pi, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  t1=MPI_Wtime();
  if(rank==0){
    printf("elapsed time is %.4f seconds\n", t1-t0);
    printf("PI is approximately %.7f\n", pi);
  }

  MPI_Finalize();
  return 0;
}
