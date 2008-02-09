/* shift 1d or 2d array */
/*
  Copyright (C) 2006 Colorado School of Mines
  
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

#include "fft2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;
    int  ompchunk; 
    float del1,del2;

    sf_file Fi=NULL; /* input   */
    sf_file Fo=NULL; /* output  */

    float      **rr=NULL;  /* data */
    sf_complex **cc=NULL;

    /* cube axes */
    sf_axis a1,a2;
    int     i1,i2;
/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */

    if(! sf_getfloat("del1",&del1))   del1=0.;   /* delay on axis 1 */
    if(! sf_getfloat("del2",&del2))   del2=0.;   /* delay on axis 2 */

    Fi = sf_input ("in" );
    Fo = sf_output("out");

    /* input axes */
    a1 = sf_iaxa(Fi,1);
    a2 = sf_iaxa(Fi,2);
    if(verb) sf_raxa(a1);
    if(verb) sf_raxa(a2);

    fft2_init(sf_n(a1),sf_n(a2));

/*------------------------------------------------------------*/

    /* allocate arrays */
    rr = sf_floatalloc2   (sf_n(a1),sf_n(a2)); /*    data */
    cc = sf_complexalloc2 (sf_n(a1),sf_n(a2));

/*------------------------------------------------------------*/
    sf_floatread(rr[0],sf_n(a1)*sf_n(a2),Fi); /* read data */
    for(    i2=0;i2<sf_n(a2);i2++) {
	for(i1=0;i1<sf_n(a1);i1++) {
	    cc[i2][i1] = rr[i2][i1];
	}
    }
    fft2(false,(kiss_fft_cpx**) cc);

/*------------------------------------------------------------*/

    /* init shift */
    sft2_init(-del1,sf_d(a1),
	      -del2,sf_d(a2));
    
    /* apply shift */
    sft2(cc);

/*------------------------------------------------------------*/

    fft2( true,(kiss_fft_cpx**) cc);
    for(    i2=0;i2<sf_n(a2);i2++) {
	for(i1=0;i1<sf_n(a1);i1++) {
	    rr[i2][i1] = cc[i2][i1];
	}
    }
    sf_floatwrite(rr[0],sf_n(a1)*sf_n(a2),Fo); /* write data */
/*------------------------------------------------------------*/

    exit (0);
}
