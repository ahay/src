/* OpenMP example, summation of vectors c = a + b */
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

int main (int argc, char *argv[]) {
    int n1, n2, i, j;
    float *a, *b, *c;

    sf_file ain, bin, cout = NULL;

    sf_init (argc, argv);
    ain = sf_input ("in"); /* Input vector a */
    bin = sf_input ("b"); /* Input vector b */
    if (SF_FLOAT != sf_gettype (ain) ||
        SF_FLOAT != sf_gettype (bin))
        sf_error ("Need float");
    /* Vector size */
    if (!sf_histint (ain, "n1", &n1)) sf_error ("No n1=");
    /* Number of vectors */
    n2 = sf_leftsize (ain, 1);
    /* Output vector */
    cout = sf_output ("out");
    /* Vectors in memory */
    a = sf_floatalloc (n1); b = sf_floatalloc (n1);
    c = sf_floatalloc (n1);
    /* Outer loop over vectors */
    for (i = 0; i < n2; i++) {
        sf_floatread (a, n1, ain);
        sf_floatread (b, n1, bin);
        /* Parallel summation */
#pragma omp parallel for private(j) shared(a,b,c)
        for (j = 0; j < n1; j++)
            c[j] = a[j] + b[j];
        sf_floatwrite (c, n1, cout);
    }
    sf_fileclose (ain); sf_fileclose (bin);
    sf_fileclose (cout);
    return 0;
}

