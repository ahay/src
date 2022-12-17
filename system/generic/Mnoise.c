/* Add random noise to the data.

July 2011 program of the month:
http://ahay.org/blog/2011/07/03/program-of-the-month-sfnoise/
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2021 Colorado School of Mines
  
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

#include <math.h>
#include <time.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
  float   mean, var, range, a, b;
  size_t  nbuf, nsiz;
  int     seed;
  size_t  i;
  bool    normal, rep;
  sf_file in, out;

  bool isreal;      /* true for float; false for complex */
  float      *datR; /*    real data */
  sf_complex *datC; /* complex data */
  float      re,im;

  sf_init (argc,argv);
  in  = sf_input ( "in");
  out = sf_output("out");

  /* ------------------------------------------------------------ */
  if (SF_FLOAT == sf_gettype(in))
    isreal = true;
  else
    isreal = false;

  /* ------------------------------------------------------------ */
  /* random seed */
  if (!sf_getint("seed",&seed)) seed = time(NULL);
  init_genrand((unsigned long) seed);

  /* ------------------------------------------------------------ */
  if (!sf_getbool("type",&normal)) normal=true;
  /* noise distribution, y: normal, n: uniform */

  if (!sf_getfloat("var",&var)) {
	/* noise variance */

  	if (!sf_getfloat("range",&range)) {
  	    /* noise range (default=1) */
  	    a = 1.;
  	} else {
  	    a = normal ? 2.*range/9. : 2.*range;
  	}

  } else {
	   a = normal ? sqrtf(var) : sqrtf(12*var);
  }

  /* ------------------------------------------------------------ */
  if (!sf_getfloat("mean",&mean)) mean=0;
  /* noise mean */
  b = normal ? mean : mean - 0.5*a;

  /* ------------------------------------------------------------ */
  if (!sf_getbool("rep",&rep)) rep=false;
  /* if y, replace data with noise */

  /* ------------------------------------------------------------ */
  if(isreal) {
    nbuf = BUFSIZ/sizeof(float);
    datR = sf_floatalloc(nbuf);
    datC = NULL;
  } else {
    nbuf = BUFSIZ/sizeof(sf_complex);
    datR = NULL;
    datC = sf_complexalloc(nbuf);
  }

  /* ------------------------------------------------------------ */
  for (nsiz = sf_filesize(in); nsiz > 0; nsiz -= nbuf) {
  	if (nbuf > nsiz) nbuf = nsiz;

    if (rep) {

      if(isreal) { /* generate real noise */

        if (normal) {
          for (i=0; i < nbuf; i++)
            datR[i] = a * sf_randn_one_bm() + b;
        } else {
          for (i=0; i < nbuf; i++)
            datR[i] = a * genrand_real1() + b;
        }

      } else { /* generate complex noise */

        if (normal) {
          for (i=0; i < nbuf; i++) {
            re = a * sf_randn_one_bm() + b;
            im = a * sf_randn_one_bm() + b;
            #ifdef SF_HAS_COMPLEX_H
              datC[i]   = re+im*I;
            #else
              datC[i].r = re;
              datC[i].i = im;
            #endif
          }
        } else {
          for (i=0; i < nbuf; i++) {
            re = a * genrand_real1() + b;
            im = a * genrand_real1() + b;
            #ifdef SF_HAS_COMPLEX_H
              datC[i]   = re+im*I;
            #else
              datC[i].r = re;
              datC[i].i = im;
            #endif
          }
        }

      }

    } else {

      if(isreal) sf_floatread  (datR,nbuf,in);
      else       sf_complexread(datC,nbuf,in);

      if(isreal) { /* add real noise */

        if (normal) {
          for (i=0; i < nbuf; i++)
            datR[i] += a * sf_randn_one_bm() + b;
        } else {
          for (i=0; i < nbuf; i++)
            datR[i] += a * genrand_real1() + b;
        }

      } else { /* add complex noise */

        if (normal) {
          for (i=0; i < nbuf; i++) {
            re = a * sf_randn_one_bm() + b;
            im = a * sf_randn_one_bm() + b;
            #ifdef SF_HAS_COMPLEX_H
              datC[i]   += re+im*I;
            #else
              datC[i].r += re;
              datC[i].i += im;
            #endif
          }
        } else {
          for (i=0; i < nbuf; i++) {
            re = a * genrand_real1() + b;
            im = a * genrand_real1() + b;
            #ifdef SF_HAS_COMPLEX_H
              datC[i]   += re+im*I;
            #else
              datC[i].r += re;
              datC[i].i += im;
            #endif
          }
        }

      }

    }

    if(isreal) sf_floatwrite  (datR,nbuf,out);
    else       sf_complexwrite(datC,nbuf,out);
  }

  if(isreal) free(datR);
  else       free(datC);

  exit(0);
}
