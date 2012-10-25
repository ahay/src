#include <rsf.h>
#include "eigv.h"

/*------------------------------------------------------------*/
void eigval2( float **a,
	      float  *e)
/*<test>*/
{
/*< return eigenvalues of 2x2 matrix >*/
    float a11,a12,a21,a22;
    float d,cp,cm;

    a11 = a[0][0];
    a12 = a[0][1];
    a21 = a[1][0];
    a22 = a[1][1];

    cp = (a11+a22)/2.;
    cm = (a11-a22)/2.;

    d = sqrt( cm*cm + a12*a21 );
    e[0] = cp + d;
    e[1] = cp - d;
}

/*------------------------------------------------------------*/
void eigvec2( float **a,
	      float   e,
	      float  *u)
/*<test>*/
{
/*< return eigenvector for given eigenvalue of 2x2 matrix >*/
    float /* a11,a12, */ a21,a22;
    float uu;
    
    /*a11 = a[0][0];
      a12 = a[0][1]; */
    a21 = a[1][0];
    a22 = a[1][1];

    if(a21!=0) {
	u[1] = 1.;
	u[0] = -(a22-e)/a21;
    } else {
	u[1] = 1;
	u[0] = 0;
    }

    uu = sqrt(u[0]*u[0]+u[1]*u[1]);    
    if(uu!=0) {
	u[0] /= uu;
	u[1] /= uu;
    }
}
