/*multiply, for complex Matrix */

/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include <string.h>
#include <unistd.h>
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char * argv[])
{

    sf_file inA, outC, inB;
    int an1, an2, bn1, bn2;
    int im, in, ik, m, n, k, nth;
    
    sf_complex **a, **b, **c;
    sf_axis aax1, aax2, bax1, bax2, cax1, cax2;
    
    /* init RSF */
    sf_init (argc, argv);
    inA = sf_input("in");
    inB = sf_input("B");
    outC= sf_output("out");

    if(SF_COMPLEX != sf_gettype(inA)) sf_error("Need complex input!");
    if(SF_COMPLEX != sf_gettype(inB)) sf_error("Need complex input!");
    
    if(!sf_histint(inA, "n1", &an1)) sf_error("No n1 in input");
    if(!sf_histint(inB, "n1", &bn1)) sf_error("No n1 in input");
    if(!sf_histint(inA, "n2", &an2)) an2 = 1;
    if(!sf_histint(inB, "n2", &bn2)) bn2 = 1;

    if(sf_leftsize(inA,2)>1) sf_error("Input should be a matrix!");
    if(sf_leftsize(inB,2)>1) sf_error("Input should be a matrix!");
    
    if(an2 != bn1) sf_error("Input do not match!");

    aax1 = sf_iaxa(inA, 1);
    aax2 = sf_iaxa(inA, 2);
    bax1 = sf_iaxa(inB, 1);
    bax2 = sf_iaxa(inB, 2);
    
    cax1 = aax1;
    cax2 = bax2;
    
    a = sf_complexalloc2(an1, an2);
    b = sf_complexalloc2(bn1, bn2);
    
    sf_complexread(a[0], an1*an2, inA);
    sf_complexread(b[0], bn1*bn2, inB);

    m = an1;
    n = an2;
    k = bn2;
    
#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
      /* omp_set_num_threads(nth); */
    }
    sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    c = sf_complexalloc2(m,k);
#ifdef _OPENMP
#pragma omp parallel for private(im,ik) shared(c,a,b,an1,an2,bn2)
#endif
    for(im=0; im< an1; im++) {
	for (ik=0; ik<bn2; ik++ ) {
	    c[ik][im]=sf_cmplx(0.0,0.0);
	    for (in=0; in<an2; in++) {
#ifdef SF_HAS_COMPLEX_H
		c[ik][im]+=a[in][im]*b[ik][in];
#else
		c[ik][im]+=sf_cmul(a[in][im],b[ik][in]);
#endif
	    }
	}
    }
    sf_oaxa(outC, cax1, 1);
    sf_oaxa(outC, cax2, 2);
    
    sf_complexwrite(c[0], m*k, outC);

    exit(0);
    
}

		
		
		
		
		



    
    
    
 
    




















