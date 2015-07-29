/* Generate fractal fern. */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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
#include <time.h>

#include <rsf.h>

/*
                set 1     set 2     set 3     set 4   
          a     0.0       0.2      -0.15      0.85
          b     0.0      -0.26      0.28      0.04
          c     0.0       0.23      0.26     -0.04
          d     0.16      0.22      0.24      0.85
          e     0.0       0.0       0.0       0.0
          f     0.0       1.6       0.44      1.6
probability     0.01      0.07      0.07      0.85
*/

static const float  a1[] = {0.00, 0.20,-0.15, 0.75};
static const float  a2[] = {0.00, 0.20,-0.15, 0.85};
static const float  b [] = {0.00,-0.26, 0.28, 0.04};
static const float  c [] = {0.00, 0.23, 0.26,-0.04};
static const float  d [] = {0.16, 0.22, 0.24, 0.85};
static const float  e [] = {0.00, 0.00, 0.00, 0.00};
static const float  f [] = {0.00, 1.60, 0.44, 1.60};
static const double p1[] = {0.10, 0.18, 0.26};
static const double p2[] = {0.01, 0.07, 0.15};

int main(int argc,char *argv[])
{
    int i, j, k, n, nbuf=BUFSIZ, seed;
    float x=0., y=0., x2, y2;
    double r;
    const double *p;
    const float *a;
    bool angle;
    sf_complex *xy;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");

    if (!sf_getint("n",&n)) n=1000;
    /* number of points */
    sf_putint(out,"n1",n);

    if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */
    init_genrand((unsigned long) seed);

    if (!sf_getbool("angle",&angle)) angle=true;
    /* if y, use more angular fern */
    p = angle? p1: p2;
    a = angle? a1: a2;

    sf_setformat(out,"native_complex");
    
    nbuf /= sizeof(sf_complex);
    xy = sf_complexalloc(nbuf);
    
    for (i=j=0; i<n; i++, j++) {
	if (j==nbuf) {
	    sf_complexwrite(xy,nbuf,out);
	    j=0;
	}
	r = genrand_real1();

	if (r < p[0])
            k = 0;
	else if (r < p[1])
            k = 1;
	else if (r < p[2])
            k = 2;
	else
            k = 3;

         x2 = a[k] * x + b[k] * y + e[k];
         y2 = c[k] * x + d[k] * y + f[k];
         x = x2;
         y = y2;
	 xy[j] = sf_cmplx(x,y);
    }
    sf_complexwrite(xy,j,out);

    exit(0);
}
