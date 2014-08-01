/* Large rectangular matrix in-place transpose with MPI */
/* matrix size does not have to be dividable by processor numbers */
/*
   Copyright (C) 2013 King Abduallah University of Science and Technology

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
  
  /*
   * :TODO: 
   *  1.Extend to transposing arbitrary two axes for higher dimensional tensors 
   * ========================================================================== 
   * I/O File   :                         Description
   * --------------------------------------------------------------------------
   *  indat     :    [n2][n1]             input matrix 
   *  transp    :    [n1][n2]             transposed matrix
   * ========================================================================== 
   */


#include <rsf.h>
#include <mpi.h>

static void transpose(float *m, int w, int h);
static void transpmpi(float *localmatin, float *localmatout, 
    int n1, int n2,
    int n1local, int n2local, 
    int n1local_start, int n2local_start,
    int nprocs, int rank);

int 
main(int argc, char **argv){

  int nprocs,rank;
  sf_file Fin,Fout;
  int n1,n2,fin_name_str_len,fout_name_str_len;
  float d1,d2,o1,o2;
  char *in_dat_name,*out_dat_name;
  int n1local,n2local;
  int n1local_start,n2local_start;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank==0) {
    sf_init(argc,argv);
    Fin = sf_input("indat");
    Fout = sf_output("transp");
    if (!sf_histint(Fin,"n1",&n1))	sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2))	sf_error("No n2= in input");
    sf_histfloat(Fin,"o1",&o1);
    sf_histfloat(Fin,"d1",&d1);
    sf_histfloat(Fin,"o2",&o2);
    sf_histfloat(Fin,"d2",&d2);
    in_dat_name = sf_histstring(Fin,"in");
    out_dat_name = sf_histstring(Fout,"in");
    fin_name_str_len = strlen(in_dat_name)+1;
    fout_name_str_len = strlen(out_dat_name)+1;
    sf_putint(Fout,"n1",n2);
    sf_putint(Fout,"n2",n1);
    sf_putfloat(Fout,"d1",d2);
    sf_putfloat(Fout,"d2",d1);
    sf_putfloat(Fout,"o1",o2);
    sf_putfloat(Fout,"o2",o1);
    sf_fileflush(Fout,NULL);
  }
  MPI_Bcast(&fin_name_str_len,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&fout_name_str_len,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&n1,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&n2,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&o1,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&o2,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&d1,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&d2,1,MPI_FLOAT,0,MPI_COMM_WORLD);

  if (rank!=0) {
    in_dat_name = (char *) malloc (fin_name_str_len);
    out_dat_name = (char *) malloc (fout_name_str_len);
  }
  MPI_Bcast(in_dat_name,fin_name_str_len,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(out_dat_name,fout_name_str_len,MPI_CHAR,0,MPI_COMM_WORLD);

  n1local = n1/nprocs;
  n2local = n2/nprocs;
  if (rank < (n1%nprocs))
    n1local++;
  if (rank < (n2%nprocs))
    n2local++;

  MPI_Scan(&n2local,&n2local_start,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  n2local_start -= n2local;
  MPI_Scan(&n1local,&n1local_start,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  n1local_start -= n1local;

  float *localmat = (float *) malloc (sizeof(float)*n1*n2local);
  float *localmat_transp = (float *) malloc (sizeof(float)*n1local*n2);

  MPI_File fh0;
  MPI_Status status;
  MPI_Offset file_disp = n1*n2local_start*sizeof(MPI_FLOAT);
  MPI_File_open(MPI_COMM_WORLD,in_dat_name,MPI_MODE_RDONLY, MPI_INFO_NULL, &fh0);
  MPI_File_set_view(fh0,file_disp,MPI_FLOAT,MPI_FLOAT,"native",MPI_INFO_NULL);
  MPI_File_read_all(fh0,localmat,n1*n2local,MPI_FLOAT,&status);
  MPI_File_close(&fh0);

  transpmpi(localmat,localmat_transp,n1,n2,n1local,n2local,n1local_start,n2local_start,nprocs,rank);

  MPI_Barrier(MPI_COMM_WORLD);
  file_disp = n2*n1local_start*sizeof(MPI_FLOAT);
  MPI_File_open(MPI_COMM_WORLD,out_dat_name,MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh0);
  MPI_File_set_view(fh0,file_disp,MPI_FLOAT,MPI_FLOAT,"native",MPI_INFO_NULL);
  MPI_File_write_all(fh0,localmat_transp,n2*n1local,MPI_FLOAT,&status);
  MPI_File_close(&fh0);

  free(localmat);
  MPI_Finalize();

  return 0;
}


static void 
transpmpi( float *localmatin, float *localmatout, 
    int n1, int n2,
    int n1local, int n2local, 
    int n1local_start, int n2local_start,
    int nprocs, int rank) {
  /* Need three steps:
  1. Transpose & pack sendout buf
  2. Sendout buf by Alltoallv (not alltoall because buf-size & displacements are different)
  3. Transpose sub-bufs & transpose whole array */

  int i;
  int scounts[nprocs];
  int sdispls[nprocs];
  int rcounts[nprocs];
  int rdispls[nprocs];
  int rows[nprocs];

  /* Pre-communication phase: Transpose to make sendout buf contiguous */
  transpose(localmatin,n1,n2local);

  MPI_Allgather(&n1local,1,MPI_INT,scounts,1,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgather(&n2local,1,MPI_INT,rcounts,1,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgather(&n1local_start,1,MPI_INT,sdispls,1,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgather(&n2local_start,1,MPI_INT,rdispls,1,MPI_INT,MPI_COMM_WORLD);

  for (i=0;i<nprocs;i++) {
    rows[i] = rcounts[i];
    scounts[i] *= n2local;
    sdispls[i] *= n2local;
    rcounts[i] *= n1local; 
    rdispls[i] *= n1local;
  }

  /* Inter-processor communication to sendout bufs onto each processor */
  MPI_Alltoallv(localmatin,scounts,sdispls,MPI_FLOAT,localmatout,rcounts,rdispls,MPI_FLOAT,MPI_COMM_WORLD);

  /* Post-communication phase: 
       1). transpose sub-local buf arrived from each other processor 
       2). transpose the whole local array which is the output array */
  for (i=0;i<nprocs;i++)
    transpose( &localmatout[rdispls[i]],rows[i],n1local );
  transpose(localmatout,n1local,n2);

  return;
}

static void 
transpose(float *m, int w, int h) {
  /* Following-the-cycle algorithm */
  /* w is fast dimension  */
  int start, next, i;
  float tmp;
  for (start = 0; start <= w * h - 1; start++) {
    next = start;
    i = 0;
    do { 
      i++;
      next = (next % h) * w + next / h;
    } while (next > start);
    if (next < start || i == 1) continue;

    tmp = m[next = start];
    do {
      i = (next % h) * w + next / h;
      m[next] = (i == start) ? tmp : m[i];
      next = i;
    } while (next > start);
  }   
  return;
}







