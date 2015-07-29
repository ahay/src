#include "vector.h"
#include "alloc.h"

#include "_defs.h"
/*^*/

#include "file.h"
/*^*/

#include "point.h"
/*^*/

#ifndef _sf_vector_h

typedef struct{
    double dx,dz;
}vc2d;
/*^*/

typedef struct{
    double dx,dy,dz;
}vc3d;
/*^*/

#endif

/*------------------------------------------------------------*/
double det3(double *m)
/*< compute 3D determinant >*/
{
  double pp,mm;

  pp=m[0]*m[4]*m[8]+m[3]*m[2]*m[7]+m[1]*m[6]*m[5];
  mm=m[2]*m[4]*m[6]+m[1]*m[3]*m[8]+m[0]*m[5]*m[7];

  return pp-mm;
}

/*------------------------------------------------------------*/
double det2(double *m)
/*< compute 2D determinant >*/
{
  return m[0]*m[3]-m[1]*m[2];
}

/*------------------------------------------------------------*/
double jac3d(pt3d *C, pt3d *T, pt3d *P, pt3d *Q)
/*< 3D jacobian >*/
{
  double m[9];
  m[0]=C->x-P->x;  m[1]=C->x-Q->x;  m[2]=C->x-T->x;
  m[3]=C->y-P->y;  m[4]=C->y-Q->y;  m[5]=C->y-T->y;
  m[6]=C->z-P->z;  m[7]=C->z-Q->z;  m[8]=C->z-T->z;
  return det3(m);
}

/*------------------------------------------------------------*/
vc3d vec3d(pt3d* O, pt3d* A)
/*< build 3D vector >*/
{
    vc3d V;

    V.dx = A->x - O->x;
    V.dy = A->y - O->y;
    V.dz = A->z - O->z;

    return V;
}

vc2d vec2d(pt2d* O, pt2d* A)
/*< build 2D vector >*/
{
    vc2d V;

    V.dx = A->x - O->x;
    V.dz = A->z - O->z;

    return V;
}

/*------------------------------------------------------------*/
vc3d axa3d( int n)
/*< build 3D unit vector >*/
{
    vc3d V;
    V.dx=V.dy=V.dz=0.;

    if(n<1 || n>3) n=1;
    
    switch(n) {
	case 3:
	    V.dy=1;
	    break;		
	case 2:
	    V.dx=1;
	    break;	    
	case 1:
	    V.dz=1;
	    break;
    }
    return V;
}

/*------------------------------------------------------------*/
double scp3d(vc3d* U, vc3d* V)
/*< scalar product of 3D vectors >*/
{
    return U->dx*V->dx + U->dy*V->dy + U->dz*V->dz;
}

double scp2d(vc2d* U, vc2d* V)
/*< scalar product of 2D vectors >*/
{
    return U->dx*V->dx + U->dz*V->dz;
}

/*------------------------------------------------------------*/
vc3d vcp3d(vc3d* U, vc3d* V)
/*< vector product of 3D vectors >*/
{
    vc3d W;

    W.dx=(U->dy*V->dz) - (V->dy*U->dz);
    W.dy=(U->dz*V->dx) - (V->dz*U->dx);
    W.dz=(U->dx*V->dy) - (V->dx*U->dy);

    return W;
}

/*------------------------------------------------------------*/
double len3d(vc3d* V)
/*< 3D vector length >*/
{
    double l;
    l = sqrtf( V->dx*V->dx + 
	       V->dy*V->dy + 
	       V->dz*V->dz);
    return l;
}

double len2d(vc2d* V)
/*< 2D vector length >*/
{
    double l;
    l = sqrtf( V->dx*V->dx +
               V->dz*V->dz);
    return l;
}

/*------------------------------------------------------------*/
vc3d nor3d(vc3d* V)
/*< normalize 3D vector >*/
{
    vc3d W;
    double l;
    l = len3d(V);

    W.dx = V->dx / l;
    W.dy = V->dy / l;
    W.dz = V->dz / l;

    return W;
}

vc2d nor2d(vc2d* V)
/*< normalize 2D vector >*/
{
    vc2d W;
    double l;
    l = len2d(V);

    W.dx = V->dx / l;
    W.dz = V->dz / l;

    return W;
}

/*------------------------------------------------------------*/
double ang3d(vc3d* U, vc3d* V)
/*< angle between 3D vectors >*/
{
    double c,a;

    c = U->dx * V->dx +
	U->dy * V->dy +
	U->dz * V->dz;
    c/= len3d(U);
    c/= len3d(V);

    c = SF_SIG(c) * SF_MIN( 1., SF_ABS(c));

    a = 180*acosf(c)/SF_PI;

    return a;
}

double ang2d(vc2d* U, vc2d* V)
/*< angle between 2D vectors >*/
{
    double c,a;

    c = U->dx * V->dx +
        U->dz * V->dz;
    c/= len2d(U);
    c/= len2d(V);

    c = SF_SIG(c) * SF_MIN( 1., SF_ABS(c));

    a = 180*acosf(c)/SF_PI;

    return a;
}


/*------------------------------------------------------------*/
pt3d tip3d(pt3d* O, vc3d* V)
/*< tip of a 3D vector >*/
{
    pt3d A;

    A.x = O->x + V->dx;
    A.y = O->y + V->dy;
    A.z = O->z + V->dz;
    A.v = 0;

    return A;    
}

pt2d tip2d(pt2d* O, vc2d* V)
/*< tip of a 3D vector >*/
{
    pt2d A;

    A.x = O->x + V->dx;
    A.z = O->z + V->dz;
    A.v = 0;

    return A;
}

/*------------------------------------------------------------*/
vc3d scl3d(vc3d* V, float s)
/*< scale a 3D vector >*/
{
     vc3d W;

     W.dx = V->dx * s;
     W.dy = V->dy * s;
     W.dz = V->dz * s;

     return W;
}

vc2d scl2d(vc2d* V, float s)
/*< scale a 2D vector >*/
{
     vc2d W;

     W.dx = V->dx * s;
     W.dz = V->dz * s;

     return W;
}

/*------------------------------------------------------------*/
vc3d rot3d(vc3d *nn,
	   vc3d *aa,
	   float phi)
/*< rotate a 3D vector >*/
{
    vc3d qq;
    float Qxx,Qxy,Qxz, Qyx,Qyy,Qyz, Qzx,Qzy,Qzz;
    float    nx,ny,nz;
    float    ax,ay,az;
    float    qx,qy,qz;
    float sif,cof;

    nx=nn->dx; ny=nn->dy; nz=nn->dz;
    ax=aa->dx; ay=aa->dy; az=aa->dz;

    sif=sinf(phi);
    cof=cosf(phi);

    Qxx=nx*nx+(ny*ny+nz*nz)*cof;
    Qxy=nx*ny*(1-cof)-nz*sif;
    Qxz=nx*nz*(1-cof)+ny*sif;
    
    Qyx=ny*nx*(1-cof)+nz*sif;
    Qyy=ny*ny+(nz*nz+nx*nx)*cof;
    Qyz=ny*nz*(1-cof)-nx*sif;
    
    Qzx=nz*nx*(1-cof)-ny*sif;
    Qzy=nz*ny*(1-cof)+nx*sif;
    Qzz=nz*nz+(nx*nx+ny*ny)*cof;

    qx = Qxx*ax + Qxy*ay + Qxz*az;
    qy = Qyx*ax + Qyy*ay + Qyz*az;
    qz = Qzx*ax + Qzy*ay + Qzz*az;

    qq.dx=qx;
    qq.dy=qy;
    qq.dz=qz;
    
    return qq;
}

/*------------------------------------------------------------*/
void vc3dread1(sf_file F,
               vc3d   *v,
               size_t  n1)
/*< input vector3d 1-D vector >*/
{
    int i1;
    float **w;
    w=sf_floatalloc2(3,n1);

    sf_floatread(w[0],3*n1,F);
    for( i1=0; i1<n1; i1++) {
        v[i1].dx = w[i1][0];
        v[i1].dy = w[i1][1];
        v[i1].dz = w[i1][2];
    }

    free(*w); free(w);
}



