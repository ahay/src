#include "point.h"

#include "file.h"
/*^*/

#ifndef _sf_point_h

/* distance */
#define DST2d(A,B)   ( (B.x-A.x)*(B.x-A.x)+(B.z-A.z)*(B.z-A.z) ) 
/*^*/
/* vector product magnitude (Jacobian) */
#define JAC2d(O,A,B) ( (A.x-O.x)*(B.z-O.z)-(B.x-O.x)*(A.z-O.z) )
/*^*/

typedef struct{
    double x,z;   /* location */
    float  v;     /* value    */
}pt2d;
/*^*/

typedef struct{
    double x,y,z; /* location */
    float  v;     /* value    */
}pt3d;
/*^*/

#endif


/*------------------------------------------------------------*/

void printpt2d(pt2d P)
/*< print point info  >*/
{
    fprintf(stderr,"P info: x=%f z=%f v=%f \n",P.x,P.z,P.v);
}

void printpt3d(pt3d P)
/*< print point info  >*/
{
    fprintf(stderr,"P info: x=%f y=%f z=%f v=%f \n",P.x,P.y,P.z,P.v);
}

/*------------------------------------------------------------*/

void writept2d(sf_file F, pt2d *v, int n, int k)
/*< output point2d vector >*/
{
    int i;
    float w[k];

    for( i=0; i<n; i++) {
	w[0]          = v[i].x; 
	w[1]          = v[i].z; 
	if(k==3) w[2] = v[i].v;
	sf_floatwrite(w,k,F);
    }
}

void writept3d(sf_file F, pt3d *v, int n, int k)
/*< output point3d vector >*/
{
    int i;
    float w[k];

    for( i=0; i<n; i++) {
	w[0]          = v[i].x;
	w[1]          = v[i].y;
	w[2]          = v[i].z;
	if(k==4) w[3] = v[i].v;
	sf_floatwrite(w,k,F);
    }
}

/*------------------------------------------------------------*/

void readpt2d(sf_file F, pt2d *v, int n, int k)
/*< input point2d vector >*/
{
    int i;
    float w[k];

    for( i=0; i<n; i++) {
	sf_floatread(w,k,F);
	v[i].x = w[0]; 
	v[i].z = w[1];
	if(k==3) v[i].v = w[2];
    }
}

void readpt3d(sf_file F, pt3d *v, int n, int k)
/*< input point3d vector >*/
{
    int i;
    float w[k];

    for( i=0; i<n; i++) {
	sf_floatread(w,k,F);
	v[i].x = w[0];
	v[i].y = w[1]; 
	v[i].z = w[2];
	if(k==4) v[i].v = w[3];
    }
}

/*------------------------------------------------------------*/
