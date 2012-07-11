#include "point.h"
#include "alloc.h"

#include "file.h"
/*^*/

#ifndef _sf_point_h

#define DST2d(A,B) ( (B.x-A.x)*(B.x-A.x) + \
                     (B.z-A.z)*(B.z-A.z) ) 
#define DST3d(A,B) ( (B.x-A.x)*(B.x-A.x) + \
                     (B.y-A.y)*(B.y-A.y) + \
                     (B.z-A.z)*(B.z-A.z) ) 
/*^*/

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
/*< print point2d info  >*/
{
    fprintf(stderr,
	    "P info: x=%f z=%f v=%f \n",
	    P.x,P.z,P.v);
}

void printpt3d(pt3d P)
/*< print point3d info  >*/
{
    fprintf(stderr,
	    "P info: x=%f y=%f z=%f v=%f \n",
	    P.x,P.y,P.z,P.v);
}

/*------------------------------------------------------------*/

void pt2dwrite1(sf_file F, 
		pt2d   *v, 
		size_t  n1, 
		int     k)
/*< output point2d 1-D vector >*/
{
    int i1;
    float **w;
    w=sf_floatalloc2(k,n1);

    for( i1=0; i1<(int)n1; i1++) {
	;        w[i1][0] = v[i1].x; 
	;        w[i1][1] = v[i1].z; 
	if(k==3) w[i1][2] = v[i1].v;
    }
    sf_floatwrite(w[0],k*n1,F);

    free( *w); free(w);
}

void pt2dwrite2(sf_file F, 
		pt2d  **v, 
		size_t  n1, 
		size_t  n2,
		int     k)
/*< output point2d 2-D vector >*/
{
    int i1,i2;
    float ***w;
    w=sf_floatalloc3(k,n1,n2);

    for( i2=0; i2<(int)n2; i2++) {
	for( i1=0; i1<(int)n1; i1++) {
	    ;        w[i2][i1][0] = v[i2][i1].x; 
	    ;        w[i2][i1][1] = v[i2][i1].z; 
	    if(k==3) w[i2][i1][2] = v[i2][i1].v;
	}
    }
    sf_floatwrite(w[0][0],k*n1*n2,F);
    
    free(**w); free(*w); free(w);
}

/*------------------------------------------------------------*/

void pt3dwrite1(sf_file F, 
		pt3d   *v, 
		size_t  n1, 
		int     k)
/*< output point3d 1-D vector >*/
{
    int i1;
    float **w;
    w=sf_floatalloc2(k,n1);

    for( i1=0; i1<(int)n1; i1++) {
	;        w[i1][0] = v[i1].x;
	;        w[i1][1] = v[i1].y;
	;        w[i1][2] = v[i1].z;
	if(k==4) w[i1][3] = v[i1].v;
    }
    sf_floatwrite(w[0],k*n1,F);

    free(*w); free(w);
}

void pt3dwrite2(sf_file F, 
		pt3d  **v, 
		size_t  n1, 
		size_t  n2,
		int     k)
/*< output point3d 2-D vector >*/
{
    int i1,i2;
    float ***w;
    w=sf_floatalloc3(k,n1,n2);

    for( i2=0; i2<(int)n2; i2++) {
	for( i1=0; i1<(int)n1; i1++) {
	    ;        w[i2][i1][0] = v[i2][i1].x;
	    ;        w[i2][i1][1] = v[i2][i1].y;
	    ;        w[i2][i1][2] = v[i2][i1].z;
	    if(k==4) w[i2][i1][3] = v[i2][i1].v;
	}
    }
    sf_floatwrite(w[0][0],k*n1*n2,F);

    free(**w); free(*w); free(w);
}

/*------------------------------------------------------------*/

void pt2dread1(sf_file F, 
	       pt2d   *v, 
	       size_t  n1, 
	       int     k)
/*< input point2d 1-D vector >*/
{
    int i1;
    float **w;
    w=sf_floatalloc2(k,n1);

    sf_floatread(w[0],k*n1,F);
    for( i1=0; i1<(int)n1; i1++) {
	;        v[i1].x = w[i1][0]; 
	;        v[i1].z = w[i1][1];
	if(k==3) v[i1].v = w[i1][2];
    }

    free(*w); free(w);
}

void pt2dread2(sf_file F, 
	       pt2d  **v, 
	       size_t  n1,
	       size_t  n2,
	       int     k)
/*< input point2d 2-D vector >*/
{
    int i1,i2;
    float ***w;
    w=sf_floatalloc3(k,n1,n2);

    sf_floatread(w[0][0],k*n1*n2,F);
    for( i2=0; i2<(int)n2; i2++) {
	for( i1=0; i1<(int)n1; i1++) {
	    ;        v[i2][i1].x = w[i2][i1][0]; 
	    ;        v[i2][i1].z = w[i2][i1][1];
	    if(k==3) v[i2][i1].v = w[i2][i1][2];
	}
    }

    free(**w); free(*w); free(w);
}

/*------------------------------------------------------------*/

void pt3dread1(sf_file F, 
	       pt3d   *v, 
	       size_t  n1, 
	       int     k)
/*< input point3d 1-D vector >*/
{
    int i1;
    float **w;
    w=sf_floatalloc2(k,n1);

    sf_floatread(w[0],k*n1,F);
    for( i1=0; i1<(int)n1; i1++) {
	;        v[i1].x = w[i1][0];
	;        v[i1].y = w[i1][1]; 
	;        v[i1].z = w[i1][2];
	if(k==4) v[i1].v = w[i1][3];
    }

    free(*w); free(w);
}

void pt3dread2(sf_file F, 
	       pt3d  **v, 
	       size_t  n1, 
	       size_t  n2,
	       int     k)
/*< input point3d 2-D vector >*/
{
    int i1,i2;
    float ***w;
    w=sf_floatalloc3(k,n1,n2);

    sf_floatread(w[0][0],k*n1*n2,F);
    for( i2=0; i2<(int)n2; i2++) {
	for( i1=0; i1<(int)n1; i1++) {
	    ;        v[i2][i1].x = w[i2][i1][0];
	    ;        v[i2][i1].y = w[i2][i1][1]; 
	    ;        v[i2][i1].z = w[i2][i1][2];
	    if(k==4) v[i2][i1].v = w[i2][i1][3];
	}
    }

    free(**w); free(*w); free(w);
}

/*------------------------------------------------------------*/

pt2d* pt2dalloc1( size_t n1)
/*< alloc point2d 1-D vector >*/
{
    pt2d *ptr;
    ptr =  (pt2d*) sf_alloc(n1,sizeof(pt2d));
    return ptr;
}

void pt2dfree1 (pt2d* tofree) 
/*< free allocated storage >*/
{
	free (tofree);
	tofree = NULL;
	return;
}

pt2d** pt2dalloc2( size_t n1, size_t n2)
/*< alloc point2d 2-D vector >*/
{
    size_t i2;
    pt2d **ptr;

    ptr = (pt2d**) sf_alloc(n2,sizeof(pt2d*));
    ptr[0] = pt2dalloc1(n1*n2);
    for (i2=1; i2<n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

pt2d*** pt2dalloc3( size_t n1,
		    size_t n2,
		    size_t n3)
/*< alloc point2d 3-D vector >*/
{
    size_t i3;
    pt2d ***ptr;

    ptr = (pt2d***) sf_alloc(n3,sizeof(pt2d**));
    ptr[0] = pt2dalloc2(n1,n2*n3);
    for (i3=1; i3<n3; i3++) {
	ptr[i3] = ptr[0]+i3*n2;
    }
    return ptr;
}

/*------------------------------------------------------------*/

pt3d* pt3dalloc1( size_t n1)
/*< alloc point3d 1-D vector >*/
{
    pt3d *ptr;
    ptr =  (pt3d*) sf_alloc(n1,sizeof(pt3d));
    return ptr;
}

pt3d** pt3dalloc2( size_t n1,
		   size_t n2)
/*< alloc point3d 2-D vector >*/
{
    size_t i2;
    pt3d **ptr;

    ptr = (pt3d**) sf_alloc(n2,sizeof(pt3d*));
    ptr[0] = pt3dalloc1(n1*n2);
    for (i2=1; i2<n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

pt3d*** pt3dalloc3( size_t n1,
		    size_t n2,
		    size_t n3)
/*< alloc point3d 3-D vector >*/
{
    size_t i3;
    pt3d ***ptr;

    ptr = (pt3d***) sf_alloc(n3,sizeof(pt3d**));
    ptr[0] = pt3dalloc2(n1,n2*n3);
    for (i3=1; i3<n3; i3++) {
	ptr[i3] = ptr[0]+i3*n2;
    }
    return ptr;
}
