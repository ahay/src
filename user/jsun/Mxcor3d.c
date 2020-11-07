/* OpenMP time- or freq-domain reversed cross-correlation on the fourth axes, read entire cube into memory */

/*
  Copyright (C) 2007 Colorado School of Mines
  Copyright (C) 2016 The University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#define rCOR(a,b) (a*b)

#ifdef SF_HAS_COMPLEX_H
#define cCOR(a,b) (crealf(conjf(a)*b))
#else
#define cCOR(a,b) (crealf(sf_cmul(conjf(a),b)))
#endif

int main(int argc, char* argv[])
{
    bool verb,rflg,rev;

    sf_file Fs,Fr,Fi;    /* I/O files -> source, receiver, image */
    sf_axis a1,a2,a3,a4,aa; /* cube axes */
    int     i1,i2,i3,i4;
    int     n1,n2,n3,n4;

    float      ****r_us=NULL,****r_ur=NULL,***ii=NULL;
    sf_complex ****c_us=NULL,****c_ur=NULL;

    float scale;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif
    
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("rev",&rev))   rev=false;  /* reverse the fourth axis of uu */

    Fs = sf_input ("in" );
    Fr = sf_input ("uu" );
    Fi = sf_output("out");

    rflg = (sf_gettype(Fs) == SF_COMPLEX)?false:true;
    if(rflg) {
	sf_warning("real input");
    } else {
	sf_warning("complex input");
	sf_settype(Fi,SF_FLOAT);
    }

    /* read axes */
    a1=sf_iaxa(Fs,1); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fs,2); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fs,3); if(verb) sf_raxa(a3);
    a4=sf_iaxa(Fs,4); if(verb) sf_raxa(a4);

    aa=sf_maxa(1,0,1); /*length,origin,sampling*/ 
    sf_setlabel(aa,""); 
    sf_setunit (aa,""); 

    n1 = sf_n(a1);
    n2 = sf_n(a2);
    n3 = sf_n(a3);
    n4 = sf_n(a4);

    sf_oaxa(Fi,a1,1);
    sf_oaxa(Fi,a2,2);
    sf_oaxa(Fi,a3,3);
    sf_oaxa(Fi,aa,4);
    ii=sf_floatalloc3(n1,n2,n3); 
    scale = 1./n4;

    for        (i3=0; i3<n3; i3++) {
        for    (i2=0; i2<n2; i2++) {
            for(i1=0; i1<n1; i1++) {
                ii[i3][i2][i1]=0.;
            }
        }
    }

    if(rflg) {
	r_us = sf_floatalloc4  (n1,n2,n3,n4);
	r_ur = sf_floatalloc4  (n1,n2,n3,n4);
    } else {
	c_us = sf_complexalloc4(n1,n2,n3,n4);
	c_ur = sf_complexalloc4(n1,n2,n3,n4);
    }

    if(rflg) {
        sf_floatread  (r_us[0][0][0],n1*n2*n3*n4,Fs);
        sf_floatread  (r_ur[0][0][0],n1*n2*n3*n4,Fr);
    } else {
        sf_complexread(c_us[0][0][0],n1*n2*n3*n4,Fs);
        sf_complexread(c_ur[0][0][0],n1*n2*n3*n4,Fr);
    }

    if(rflg) {

        for(i4=0; i4<n4; i4++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
            private(i1,i2,i3)			\
            shared( n1,n2,n3,i4,ii,r_us,r_ur)
#endif
            for        (i3=0; i3<n3; i3++) {
                for    (i2=0; i2<n2; i2++) {
                    for(i1=0; i1<n1; i1++) {
                        if (rev) ii[i3][i2][i1] += rCOR(r_us[i4][i3][i2][i1],r_ur[n4-1-i4][i3][i2][i1]);
                        else     ii[i3][i2][i1] += rCOR(r_us[i4][i3][i2][i1],r_ur[     i4][i3][i2][i1]);
                    }
                }
            }
        } /* i4 */

    } else {

        for(i4=0; i4<n4; i4++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
            private(i1,i2,i3)			\
            shared( n1,n2,n3,i4,ii,c_us,c_ur)
#endif
            for        (i3=0; i3<n3; i3++) {
                for    (i2=0; i2<n2; i2++) {
                    for(i1=0; i1<n1; i1++) {
                        if (rev) ii[i3][i2][i1] += cCOR(c_us[i4][i3][i2][i1],c_ur[n4-1-i4][i3][i2][i1]);
                        else     ii[i3][i2][i1] += cCOR(c_us[i4][i3][i2][i1],c_ur[     i4][i3][i2][i1]);
                    }
                }
            }
        } /* i4 */

    }

    if(verb) fprintf(stderr,"\n");    

    for        (i3=0; i3<n3; i3++) {
        for    (i2=0; i2<n2; i2++) {
            for(i1=0; i1<n1; i1++) {
                ii[i3][i2][i1] *=scale;
            }
        }
    }
    sf_floatwrite(ii[0][0],n3*n2*n1,Fi);

    /*------------------------------------------------------------*/
    free(**ii); free(*ii); free(ii);

    if(rflg) {
        free(***r_us); free(**r_us); free(*r_us); free(r_us);
        free(***r_ur); free(**r_ur); free(*r_ur); free(r_ur);
    } else {
        free(***c_us); free(**c_us); free(*c_us); free(c_us);
        free(***c_ur); free(**c_ur); free(*c_ur); free(c_ur);
    }
    /*------------------------------------------------------------*/

    exit (0);
}
