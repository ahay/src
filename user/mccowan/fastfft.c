/*
  Copyright (C) 2008 Doug McCowan
   
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

#include <rsf.h>

/* FFT tables structure */

typedef struct {
    int       number_points, power;
    int       *unscr;
    float     *table;
} fftTables;

static fftTables *ffttables;

void fft(float *real, float *imag, bool transform)
/*<*************************************************************
 *
 *  fft : Fast fourier transform. Input data must be
 *          a power of two.
 *
 *  INPUT
 *  real       :   array of real data values
 *  imag       :   array of imaginary data values
 *  logn       :   nth power of 2
 *  transform  :   +1, -1 for forward or inverse fft
 *
 *  OUTPUT
 *  real       :   array of real data values
 *  imag       :   array of imaginary values
 *
 *
 *  received from : Tom Botts, PetroCorp, Oklahoma City, Oklahoma
 *
 ***************************************************************>*/     
 {
     int   npoints;
     float d, e, x1, y1, x2, y2, r1, r2;
     float x, y, *uu, *vv, *xx, *yy;
     int   i, ikj, index, j, jj, jm;
     int   k, kkk, l, m, mm, mmm, npt;
     int   sines;
     
     npoints = ffttables->number_points; 
     sines = ( ( transform ? npoints: npoints*3 ) ) >> 2;
     npt = npoints;
     l = npoints >> 1;
     m = l;
     
     /* The following 2 "for" loops take sum and differences of the data
      * which can be factored out of the fft process and saves having to
      * perform complex multiplications in the main loop.
      */
     for( i=0; i<ffttables->power; ++i ) {
	 for( j = 0; j < m; ++j ) {
	     k = j + m;
	     xx = real + j;
	     yy = real + k;
	     x = *xx + *yy;
	     *yy = *xx - *yy;
	     *xx = x;
	     xx = imag + j;
	     yy = imag + k;
	     x = *xx + *yy;
	     *yy = *xx - *yy;
	     *xx = x;
	 }
	 m >>= 1;
     }
     
     for( mmm=2; mmm <= ffttables->power; ++mmm )
     {
	 j = 1 << ( ffttables->power - mmm );
	 kkk = 1 << ( mmm - 1 );
	 k = j << 1;
	 jm = k;
	 index = 2;
	 for( mm = 2; mm <= kkk; ++mm )
	 {
	     ikj = ffttables->unscr[ index ];
	     x1 = ffttables->table[ ikj ];
	     y1 = ffttables->table[ ikj + sines ];
	     ikj = k + j;
	     r1 = x1 + y1;
	     r2 = x1 - y1;
	     for( npt = 0; npt < j; ++npt )
	     {
		 i = ikj + npt;
		 jj = k + npt;
		 xx = real + i;
		 yy = imag + i;
		 uu = real + jj;
		 vv = imag + jj;
		 
		 x2 = *uu;
		 y2 = *vv;
		 d = ( *xx + *yy )*y1;
		 e = r1*(*xx) - d;
		 d = d + r2*(*yy);
		 *xx = x2 - e;
		 *yy = y2 - d;
		 *uu = x2 + e;
		 *vv = y2 + d;
	     }
	     k = jm*mm;
	     index += 2;
	 }
     }
     
     /* fft is completed. Unscramble the results */
     /*  normalize on positive transforms */
     if(! transform )
     {
	 x = 1./npoints;
	 for( i = 0; i < npoints; ++i )
	 {
	     real[ i ] = real[ i ] * x;
	     imag[ i ] = imag[ i ] * x;
	 }
     }
     for( l = 1; l < npoints; ++l )
     {
	 k = ffttables->unscr[ l ];
	 if( k > l )
	 {
	     xx = real + k;
	     yy = imag + k;
	     uu = real + l;
	     vv = imag + l;
	     x = *xx;
	     y = *yy;
	     *xx = *uu;
	     *yy = *vv;
	     *uu = x;
	     *vv = y;
	 }
     }
 } /* end of fft */

void fft_close(void)
/*< deallocate tables >*/
{
    free (ffttables->unscr);
    free (ffttables->table);
    free (ffttables);
}

int fft_init(int n1)
/*<***************************************************************************
 *
 *  FftTables : precompute tables for Fast fourier transform. 
 *
 *    Written for Uneeda Migration Company  by  E. Coli             11/3/98
 *
 ***************************************************************************>*/ 
{
    float cosa, cosb, sina, sinb;
    float r1, r2, x;
    int   i, j, npow, base;
    int   k, n, npp, npt;
    int   q1, q2, q3, q4, q5, q6, q7;
    
    base=1;
    for (npow=0; base < n1; npow++) {
	base *= 2;
    }
    
    ffttables  = (fftTables *) sf_alloc(1,sizeof(fftTables ));
    ffttables->number_points = base;
    ffttables->power = npow;
    ffttables->table = sf_floatalloc(base + 3*base/4 + 1);
    ffttables->unscr = sf_intalloc(base);
    
    /* Now build a table of values from 0 to npoints where the low
     * order logn bits are reversed. This is an unscramble table and
     * also points to the correct trig values to use in the main loop
     * of the actual fft. */
    for(i=0; i<ffttables->number_points; ++i) {
	k = i;
	n = 0;
	for(j=0; j<ffttables->power; ++j) {
	    n = n << 1;
	    n = n | (k & 1 ? 1:0);
	    k = k >> 1;
	}
	ffttables->unscr[i] = n;
    }
    
    /* Now build a cos wave of 1 and 3/4 cycles. 
     * The Hz is 1/npoints. */
    npp = ffttables->number_points/4;
    npt = npp;
    x = 2*SF_PI/(npp*4);
    cosb = cosf(x);
    sinb = sinf(x);
    cosa = 1.0;
    sina = 0.0;
    
    q2 = npp*2;
    q4 = npp*4;
    q6 = npp*6;
    
    /* set the table with values not computed */
    ffttables->table[ 0 ]  =  1.0;
    ffttables->table[ q2 ] = -1.0;
    ffttables->table[ q4 ] =  1.0;
    ffttables->table[ q6 ] = -1.0;
    
    q3 = q2 + npp;
    q5 = q4 + npp;
    q7 = q6 + npp;
    
    ffttables->table[ npt ] = 0;
    ffttables->table[ q3 ]  = 0;
    ffttables->table[ q5 ]  = 0;
    ffttables->table[ q7 ]  = 0;
    
    --q2;
    --q4;
    --q6;
    q3 = q2 + 2;
    q5 = q4 + 2;
    q7 = q6 + 2;

    /* build 1 and 3/4 cycles */
    r1 = cosb + sinb;
    r2 = cosb - sinb;
    for( q1 = 1; q1 < npp; ++q1 ) {
	
	/* Compute the following 
	 * (cosa+isina)*(cosb+isinb) */
	x = ( cosa + sina )*sinb;
	cosa = r1*cosa - x;
	sina = x + r2*sina;
	
	ffttables->table[ q1 ] = cosa;
	ffttables->table[ q4 ] = cosa;
	ffttables->table[ q5 ] = cosa;
	cosa = -cosa;
	ffttables->table[ q2 ] = cosa;
	ffttables->table[ q3 ] = cosa;
	ffttables->table[ q6 ] = cosa;
	ffttables->table[ q7 ] = cosa;
	cosa = -cosa;
	
	++q3;
	++q5;
	++q7;
	--q2;
	--q4;
	--q6;
    }
    
    return base;

} /* end fft_init */

#ifdef TEST
int main(void) 
{
    int i, base=16;
    float a[16]={4.,1.,8.,-4.,2.,-9.,1.,-2.,6.,-6.,-3.,2.,-7.,5.,-2.,8.};
    float b[16]={7.,-2.,-5.,-2.,9.,2.,4.,-8.,-1.,6.,3.,7.,4.,-2.,-9.,1.};
    
/* allocate and precompute fft tables */
    
    fft_init(base);
    
/* print em out */
    
    printf("a=");
    for(i=0; i<base; i++) {
	if(!(i%10)) printf("\n");
	printf("%10.2e", a[i]);
    }
    printf("\n");
    
    printf("b=");
    for(i=0; i<base; i++) {
	if(!(i%10)) printf("\n");
	printf("%10.2e", b[i]);
    }
    printf("\n");
    
/* FFT them - direction */ 
    
    fft( a, b, true);
    
/* print em out */
    
    printf("a=");
    for(i=0; i<base; i++) {
	if(!(i%10)) printf("\n");
	printf("%10.2e", a[i]);
    }
    printf("\n");
    
    printf("b=");
    for(i=0; i<base; i++) {
	if(!(i%10)) printf("\n");
	printf("%10.2e", b[i]);
    }
    printf("\n");
    
/* FFT them + direction */ 
    
    fft( a, b, false);
    
/* print em out */
    
    printf("a=");
    for(i=0; i<base; i++) {
	if(!(i%10)) printf("\n");
	printf("%10.2e", a[i]);
    }
    printf("\n");
    
    printf("b=");
    for(i=0; i<base; i++) {
	if(!(i%10)) printf("\n");
	printf("%10.2e", b[i]);
    }
    printf("\n");
    
    exit(0);
}
#endif
