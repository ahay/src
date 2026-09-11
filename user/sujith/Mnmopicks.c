/* Manual NMO picks pairs (ti, vi) are interpolated to 1D 
USAGE: 
    sfnmopicks o1=0 d1=.002 n1=501 t=0,0.5,1.0 v=1500,2000,2500 > picks.rsf
*/
/*
  Copyright (C) 2026 University of Texas at Austin
  
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

int parse_floats(const char* key, float** arr) {
    char *str, *str_copy, *token;
    int count = 0, i = 0;
    
    /* Get the comma separated string from the command line */
    str = sf_getstring(key);
    if (str == NULL) return 0;
    
    str_copy = strdup(str);
    if (str_copy == NULL) sf_error("Memory allocation failed");
    
    /* count the number of elements */
    token = strtok(str_copy, ",");
    while (token != NULL) {
        count++;
        token = strtok(NULL, ",");
    }
    free(str_copy);
    
    if (count == 0) return 0;
    
    /* Allocate the exact amount of memory needed */
    *arr = sf_floatalloc(count);
    
    /* extract the floats */
    str_copy = strdup(str);
    token = strtok(str_copy, ",");
    while (token != NULL) {
        (*arr)[i++] = atof(token);
        token = strtok(NULL, ",");
    }
    free(str_copy);
    
    return count;
}

int main(int argc, char* argv[])
{
    int n1, nt, nv, i, k;
    float d1, o1, ti, vi;
    float *t, *v, *vel;
    sf_file out;

    sf_init(argc, argv);

    out = sf_output("out");

    if (!sf_getint("n1", &n1)) sf_error("Need n1=");
    if (!sf_getfloat("d1", &d1)) sf_error("Need d1=");
    if (!sf_getfloat("o1", &o1)) o1 = 0.0; /* default to 0 if not provided */

    nt = parse_floats("t", &t);
    nv = parse_floats("v", &v);

    if (nt == 0 || nv == 0) sf_error("Need t= and v= arrays (e.g., t=0,1,2 v=1500,1600,1700)");
    if (nt != nv) sf_error("t array length (%d) and v array length (%d) must match", nt, nv);

    /* Ensure time picks are strictly increasing */
    for (k = 1; k < nt; k++) {
        if (t[k] <= t[k-1]) sf_error("Time picks must be strictly increasing");
    }

    vel = sf_floatalloc(n1);

    /* Linear Interpolation */
    for (i = 0; i < n1; i++) {
        ti = o1 + i * d1; /* Current time on the grid */

        if (ti <= t[0]) {
            vi = v[0]; /* Constant extrapolation before first pick */
        } else if (ti >= t[nt-1]) {
            vi = v[nt-1]; /* Constant extrapolation after last pick */
        } else {
            /* Find the bounding interval */
            for (k = 0; k < nt - 1; k++) {
                if (ti >= t[k] && ti <= t[k+1]) {
                    /* Standard linear interpolation */
                    vi = v[k] + (v[k+1] - v[k]) * (ti - t[k]) / (t[k+1] - t[k]);
                    break;
                }
            }
        }
        vel[i] = vi;
    }

    sf_putint(out, "n1", n1);
    sf_putfloat(out, "d1", d1);
    sf_putfloat(out, "o1", o1);
    sf_putstring(out, "label1", "Time");
    sf_putstring(out, "unit1", "s");

    sf_floatwrite(vel, n1, out);

    free(t);
    free(v);
    free(vel);

    return 0;
}