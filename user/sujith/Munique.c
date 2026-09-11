/* Unique values in trace
   1-d only, for headers
   USAGE : < tfile.rsf sheadermath output="key" | sfwindow | sfdd type=float | sfunique > out.rsf
*/

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

#include <rsf.h>

int compare_floats(const void *a, const void *b) {
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

int main(int argc, char* argv[])
{
    sf_file in, out;
    int n1, n2, i, count;
    float *data, *unique;

    sf_init(argc, argv);

    in  = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) {
        sf_error("Input must be of type float. Use sfdd to convert if necessary.");
    }

    if (!sf_histint(in, "n1", &n1)) sf_error("No n1 in input");
    
    if (sf_histint(in, "n2", &n2) && n2 > 1) {
        sf_error("Input must be 1D. n2=%d detected.", n2);
    }

    data = sf_floatalloc(n1);
    sf_floatread(data, n1, in);

    qsort(data, n1, sizeof(float), compare_floats);

    unique = sf_floatalloc(n1);
    unique[0] = data[0];
    count = 1;

    for (i = 1; i < n1; i++) {
        if (data[i] != data[i-1]) {
            unique[count] = data[i];
            count++;
        }
    }

    sf_putint(out, "n1", count);
    /* o1 and d1 lose physical meaning after unique, set to defaults */
    sf_putfloat(out, "o1", 0.0);
    sf_putfloat(out, "d1", 1.0);
    sf_putstring(out, "label1", "Unique Values");

    sf_floatwrite(unique, count, out);

    return 0;
}
