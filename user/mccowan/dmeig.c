/*
  Copyright (C) 2007 Doug McCowan
   
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

void dmeig(int n       /* order of matrix */, 
	   float *a    /* [n*n] matrix to be diagonalized, destroyed in calculation */,
	   float *eval /* [n] computed eigenvalues in ascending order */, 
	   float *evec /* [n*n] corresponding eigenvectors, normalized to unit length */) 
/*< World's shortest eigenvalue - eigenvector routine 
  c     works only on symmetric, positive definite matrices 
  c     uses escalator method of Morris (1947)
  c     change root estimator to work on general symmetric matrices 
  c     Written Oct, 1976 by D.W. McCowan (LLab) >*/
{
    int   i, j, k, l, ij, ji, ik, ki, kj, lj, lk, itr;
    float anynum, d, sum, sume, amu, fun, funpr, u, fact;
    
    
/*     initialize by doing i=0 case   */
    
    eval[0] = a[0];
    evec[0] = 1.0f; 
   
/*    main escalator loop    */
   
    for(i=1; i<n; i++) {
	
/*   d = a[i,i]   */
	
	d = a[(n+1)*i]; 
   
/*              t
       compute c  * x   and move old eigenvectors into a 
                     j                                  i-1 
*/
	sume = 0.0f;
	ji = i*n; ij = i; 
	for(j=0; j<i; j++) {
	    sum = 0.0f; 
	    ki = i*n; kj = j*n;
	    for(k=0; k<i; k++) {
		a[kj] = evec[kj];
		sum += a[ki]*evec[kj];
		evec[kj] = 0.0f; 
		kj++; ki++;
	    }
	    evec[ji] = 0.0f; 
	    sume += eval[j]; 
	    a[ij] = sum;
	    ji++; ij += n;
	}
	
/*     Loop for computing each eigenvalue and corresponding eigenvector */
   
	for(j=0; j<=i; j++) {
	    if(j==0) amu = 0.5f*eval[0];
	    else if(j==i) amu = 0.5f*(eval[i-1]+sume+d); 
	    else amu = 0.5f*(eval[j-1]+eval[j]); 
   
/*     root finder iteration loop   */

	    for(itr=0; itr<10; itr++) {   
/*   
      compute secular function and its derivative 
    
                  i-1         2 
      fun = d-mu+ sum( ctxj(k)  / (mu-lambda  ))
                  k=0                       k 
    
                  i-1                          2
      funpr = -1- sum( ctxj(k) / (mu-lambda  )) 
                  k=0                      k
*/   
		fun = d-amu; 
		funpr = -1.0f;
		ik = i;
		for(k=0; k<i; k++) {
		    fun += a[ik]*a[ik]/(amu-eval[k]); 
		    anynum = a[ik]/(amu-eval[k]);
		    funpr -= anynum*anynum;
		    ik += n; 
		}
   
/*     Newton's method    */
   
		amu -= fun/funpr;
	    }     /*   end loop itr over Newton iterations   */
/*   
      Compute new eigenvectors
                 i
      y(l) = u* sum( e  (l) * ctxj(k) / (mu-lambda  ))
                k=0   k                           k 
*/   
	    u = 1.0f/sqrtf(SF_ABS(funpr));
	    ik = i;  
	    for(k=0; k<i; k++) {
		fact = u*a[ik]/(amu-eval[k]);
		lj = j*n; lk = k*n;
		for(l=0; l<i; l++) {  
		    evec[lj] += fact*a[lk];
		    lj++; lk++;
		}
		ik += n;
	    }
	    ij = i+j*n; ji = j+i*n;
	    evec[ij] = u; 
	    a[ji] = amu;
	}                     /*   end loop j over new eigenvalues   */
   
/*     Move new eigenvalues   */
   
	ji = i*n;
	for(j=0; j<=i; j++) {
	    eval[j] = a[ji];
	    ji++;
	}
    }          /*   end escalator loop i   */
   
    return;
} 

#ifdef TEST
#include <stdio.h>

int main(void) {
    int  i, j, k, n=3;
    float a[9] = {4, -1, 2, -1, 3, -2, 2, -2, 2};
    float eval[3], evec[9], b[9], test[3], sum;
    
    void dmeig(int n, float *a, float *eval, float *evec); 
    
    printf(" a=");
    for(i=0; i<n*n; i++) {
        if(!(i%n)) printf("\n");
        printf(" %f", a[i]);
        b[i] = a[i];
    }
    printf("\n");
    
    dmeig(n, a, eval, evec); 
    
    printf("\n eval=");
    for(i=0; i<n; i++) {
        if(!(i%n)) printf("\n");
        printf(" %f", eval[i]);
    }
    printf("\n");
    
    printf("\n evec=");
    for(i=0; i<n*n; i++) {
        if(!(i%n)) printf("\n");
        printf(" %f", evec[i]);
    }
    printf("\n");
    
    for(k=0; k<n; k++) {
        printf("\n Test vector for k= %d is", k);
        for(i=0; i<n; i++) {
            sum = 0.0f;
            for(j=0; j<n; j++) 
                sum += b[i+j*n]*evec[j+k*n];
            test[i] = sum/eval[k];
            if(!(i%n)) printf("\n");
            printf(" %f", test[i]);
        }
        printf("\n");
    }
}
#endif
