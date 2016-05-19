/* MPI example, summation of vectors c = a + b */
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
#include <mpi.h>

int main (int argc, char *argv[]) {
    int n1, n2, nc, esize, i, j, k = 0;
    float *a, *b, **c;

    sf_file ain, bin, cout = NULL;

    int cpuid; /* CPU id */
    int ncpu; /* Number of CPUs */
    MPI_Status mpi_stat;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size (MPI_COMM_WORLD, &ncpu);

    sf_init (argc, argv);
    ain = sf_input ("--input"); /* Input vector a */
    bin = sf_input ("b"); /* Input vector b */
    if (SF_FLOAT != sf_gettype (ain) ||
        SF_FLOAT != sf_gettype (bin))
        sf_error ("Need float");
    /* Size of an element */
    if (!sf_histint (ain, "esize", &esize))
        esize = sizeof(float);
    /* Vector size */
    if (!sf_histint (ain, "n1", &n1)) sf_error ("No n1=");
    /* Total number of vectors */
    n2 = sf_leftsize (ain, 1);
    /* Only the first CPU will do output */
    if (0 == cpuid) {
        cout = sf_output ("--output");
        sf_putint (cout, "n1", n1);
        sf_putint (cout, "n2", n2);
        sf_warning ("Running on %d CPUs", ncpu);
    }
    /* Input vectors in memory */
    a = sf_floatalloc (n1);
    b = sf_floatalloc (n1);
    /* How many vectors per CPU */
    nc = (int)(n2/(float)ncpu + 0.5f);
    c = sf_floatalloc2 (n1, nc);
    /* Starting position in input files */
    sf_seek (ain, n1*cpuid*esize, SEEK_CUR);
    sf_seek (bin, n1*cpuid*esize, SEEK_CUR);
    for (i = cpuid; i < n2; i += ncpu, k++) {
        /* Read local portion of input data */
        sf_floatread (a, n1, ain);
        sf_floatread (b, n1, bin);
        /* Parallel summation here */
        for (j = 0; j < n1; j++)
            c[k][j] = a[j] + b[j];
        /* Move on to the next portion */
        sf_seek (ain, n1*(ncpu - 1)*esize, SEEK_CUR);
        sf_seek (bin, n1*(ncpu - 1)*esize, SEEK_CUR);
    }

    if (0 == cpuid) { /* Collect results from all nodes */
        for (i = 0; i < n2; i++) {
            k = i / ncpu; /* Iteration number */
            j = i % ncpu; /* CPU number to receive from */
            if (j) /* Receive from non-zero CPU */
                MPI_Recv (&c[k][0], n1, MPI_FLOAT, j, j,
                          MPI_COMM_WORLD, &mpi_stat);
            sf_floatwrite (c[k], n1, cout);
        }
        sf_fileclose (cout);
    } else { /* Send results to CPU #0 */
        for (i = 0; i < k; i++) /* Vector by vector */
            MPI_Send (&c[i][0], n1, MPI_FLOAT, 0, cpuid,
                      MPI_COMM_WORLD);
    }
    sf_fileclose (ain); sf_fileclose (bin);
    MPI_Finalize ();
    return 0;
}

