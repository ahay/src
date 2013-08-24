#include <rsf.h>
#include "cemdutil1.h"

#define DEFAULT_THRESHOLD 0.05
#define DEFAULT_TOLERANCE 0.05
#define DEFAULT_NBPHASES 4
#define MAX_ITERATIONS 1000
#define LIM_GMP 30000
#define NBSYM 2
#define SQUARE(A) (A*A)
#define CUBE(A) (A*A*A)
#define GREATER(A,B) ((A) >= (B) ? (A) : (B))
#define SMALLER(A,B) ((A) <  (B) ? (A) : (B))
#ifdef C99_OK
#include <complex.h>
#define COMPLEX_T double complex
#define CREAL creal
#define CIMAG cimag
#define CABS cabs
#else
#define COMPLEX_T complex_data_t
#define CREAL emd_creal
#define CIMAG emd_cimag
#define CABS emd_cabs
#endif
/*^*/

#ifndef _cemdutil1_h
//emd_complex.h
#ifndef EMD_COMPLEX_H
#define EMD_COMPLEX_H

typedef struct {
  double r;
  double i;
} complex_data_t;
/*^*/
#endif

//cextr.h
#ifndef CEXTR_H
#define CEXTR_H

/* structure used to store the local extrema of the signal */
typedef struct {
  int n_min;
  int n_max;
  int *ind_min;
  int *ind_max;
  double *x_min;
  double *ry_min;
  double *iy_min;
  double *x_max;
  double *ry_max;
  double *iy_max;
} extrema_t;
/*^*/
#endif

//interpolation.h
#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#endif

//cio.h
#ifndef EMD_IO_H
#define EMD_IO_H

typedef struct {
    float threshold;
    float tolerance;
} stop_t;
/*^*/

/* structure used to store an IMF and the associated number of iterations */
typedef struct i {
  int nb_iterations;
  COMPLEX_T *pointer;
  struct i *next;
} imf_t;
/*^*/

/* structure of the IMF list */
typedef struct {
  imf_t *first;
  imf_t *last;
  int m; 
  int n; 
} imf_list_t;
/*^*/
#endif

//clocal_mean.h
#ifndef CLOCAL_MEAN_H
#define CLOCAL_MEAN_H

/* structure used to store envelopes and temporary data */
typedef struct {
  int n;
  double *re_min;
  double *ie_min;
  double *re_max;
  double *ie_max;
  double *tmp1;
  double *tmp2;
} envelope_t;
/*^*/
#endif

#endif

//cextr.c
/************************************************************************/
/*                                                                      */
/* INITIALIZATION OF EXTREMA STRUCTURE                                  */
/*                                                                      */
/************************************************************************/

extrema_t init_extr(int n)
/*<initialization for extremas extraction>*/
{
    extrema_t ex;
    ex.ind_min=(int *)malloc(n*sizeof(int));
    ex.ind_max=(int *)malloc(n*sizeof(int));
    ex.x_min=(double *)malloc(n*sizeof(double));
    ex.x_max=(double *)malloc(n*sizeof(double));
    ex.ry_min=(double *)malloc(n*sizeof(double));
    ex.ry_max=(double *)malloc(n*sizeof(double));
    ex.iy_min=(double *)malloc(n*sizeof(double));
    ex.iy_max=(double *)malloc(n*sizeof(double));
    return ex;
}

/************************************************************************/
/*                                                                      */
/* DETECTION OF LOCAL EXTREMA                                           */
/*                                                                      */
/************************************************************************/

void extr(double x[], COMPLEX_T z[], double phi,int n,extrema_t *ex)
/*< extract extremas >*/
{
    int cour;
    double val,valp,valn;
    #ifdef C99_OK
    COMPLEX_T eiphi;
    #endif
    ex->n_min=0;
    ex->n_max=0;
    
    #ifdef C99_OK
    eiphi = cexp(I*phi);
    #endif
    #ifdef C99_OK
    val = CREAL(eiphi*z[0]);
    #else
    val = crealeiphi(phi,z[0]);
    #endif

    #ifdef C99_OK
    valn = CREAL(eiphi*z[1]);
    #else
    valn = crealeiphi(phi,z[1]);
    #endif

    /* search for extrema in direction phi*/
    for(cour=1;cour<(n-1);cour++) {
        valp = val;
        val = valn;
        #ifdef C99_OK
        valn = CREAL(eiphi*z[cour+1]);
        #else
        valn = crealeiphi(phi,z[cour+1]);
        #endif

        if (val<valp && val<valn) /* local minimum */ {
            ex->x_min[ex->n_min+NBSYM]=x[cour];
            ex->ry_min[ex->n_min+NBSYM]=CREAL(z[cour]);
            ex->iy_min[ex->n_min+NBSYM]=CIMAG(z[cour]);
            ex->ind_min[ex->n_min+NBSYM]=cour;
            ex->n_min++;
        }
        if (val>valp && val>valn) /* local maximum */ {
            ex->x_max[ex->n_max+NBSYM]=x[cour];
            ex->ry_max[ex->n_max+NBSYM]=CREAL(z[cour]);
            ex->iy_max[ex->n_max+NBSYM]=CIMAG(z[cour]);
            ex->ind_max[ex->n_max+NBSYM]=cour;
            ex->n_max++;
        }
    }
}

/************************************************************************/
/*                                                                      */
/* EXTRAPOLATION OF EXTREMA TO LIMIT BORDER EFFECTS                     */
/*                                                                      */
/************************************************************************/

void boundary_conditions(double x[],COMPLEX_T z[],double phi,int n,extrema_t *ex) 
/*< setting boundary conditions >*/
{
    int cour,nbsym;
    #ifdef C99_OK
    COMPLEX_T eiphi;
    #endif
    #ifdef C99_OK
    eiphi = cexp(I*phi);
    #endif
    
    nbsym = NBSYM;
    /* reduce the number of symmetrized points if there is not enough extrema */
    while(ex->n_min < nbsym+1 && ex->n_max < nbsym+1) nbsym--;
    if (nbsym < NBSYM) {
        for(cour=0;cour<ex->n_max;cour++) {
            ex->ind_max[nbsym+cour] = ex->ind_max[NBSYM+cour];
            ex->x_max[nbsym+cour] = ex->x_max[NBSYM+cour];
            ex->ry_max[nbsym+cour] = ex->ry_max[NBSYM+cour];
            ex->iy_max[nbsym+cour] = ex->iy_max[NBSYM+cour];
        }
        for(cour=0;cour<ex->n_min;cour++) {
            ex->ind_min[nbsym+cour] = ex->ind_min[NBSYM+cour];
            ex->x_min[nbsym+cour] = ex->x_min[NBSYM+cour];
            ex->ry_min[nbsym+cour] = ex->ry_min[NBSYM+cour];
            ex->iy_min[nbsym+cour] = ex->iy_min[NBSYM+cour];
        }
    }
    
    /* select the symmetrized points and the axis of symmetry at the beginning of the signal*/
    if (ex->x_max[nbsym] < ex->x_min[nbsym]) { /* first = max */
        #ifdef C99_OK
        if (CREAL(eiphi*z[0]) > CREAL(eiphi*z[ex->ind_min[nbsym]])) { /* the edge is not a min */
        #else
        if (crealeiphi(phi,z[0]) > crealeiphi(phi,z[ex->ind_min[nbsym]])) { /* the edge is not a min */
        #endif

            if (2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*ex->x_max[nbsym]-ex->x_max[2*nbsym-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-cour];
                    ex->x_min[cour] = 2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
                }
            }
        } else { /* edge is a min -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-2-cour];
                ex->ry_min[cour] = ex->ry_min[2*nbsym-2-cour];
                ex->iy_min[cour] = ex->iy_min[2*nbsym-2-cour];
            }
            ex->x_min[nbsym-1] = x[0];
            ex->ry_min[nbsym-1] = CREAL(z[0]);
            ex->iy_min[nbsym-1] = CIMAG(z[0]);
        }
    } else { /* first = min */
        #ifdef C99_OK
        if (CREAL(eiphi*z[0]) < CREAL(eiphi*z[ex->ind_max[nbsym]])) { /* the edge is not a max */
        #else
        if (crealeiphi(phi,z[0]) < crealeiphi(phi,z[ex->ind_max[nbsym]])) { /* the edge is not a max */
        #endif

            if (2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*ex->x_min[nbsym]-ex->x_min[2*nbsym-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-cour];
                }
            }
        } else { /* edge is a max -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-2-cour];
                ex->ry_max[cour] = ex->ry_max[2*nbsym-2-cour];
                ex->iy_max[cour] = ex->iy_max[2*nbsym-2-cour];
            }
            ex->x_max[nbsym-1] = x[0];
            ex->ry_max[nbsym-1] = CREAL(z[0]);
            ex->iy_max[nbsym-1] = CIMAG(z[0]);
        }
    }
    
    
    (ex->n_min) += nbsym-1;
    (ex->n_max) += nbsym-1;
    
    
    /* select the symmetrized points and the axis of symmetry at the end of the signal*/
    if (ex->x_max[ex->n_max] < ex->x_min[ex->n_min]) { /* last is a min */
        #ifdef C99_OK
        if (CREAL(eiphi*z[n-1]) < CREAL(eiphi*z[ex->ind_max[ex->n_max]])) { /* the edge is not a max */
        #else
        if (crealeiphi(phi,z[n-1]) < crealeiphi(phi,z[ex->ind_max[ex->n_max]])) { /* the edge is not a max */
        #endif

            if (2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_min[ex->n_min-1-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-1-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-1-cour];
                }
            }
        } else { /* edge is a max -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_max[ex->n_max+2+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                ex->ry_max[ex->n_max+2+cour] = ex->ry_max[ex->n_max-cour];
                ex->iy_max[ex->n_max+2+cour] = ex->iy_max[ex->n_max-cour];
            }
            ex->x_max[ex->n_max+1] = x[n-1];
            ex->ry_max[ex->n_max+1] = CREAL(z[n-1]);
            ex->iy_max[ex->n_max+1] = CIMAG(z[n-1]);
        }
    } else {  /* last is a max */
        #ifdef C99_OK
        if (CREAL(eiphi*z[n-1]) > CREAL(eiphi*z[ex->ind_min[ex->n_min]])) { /* the edge is not a min */
        #else
        if (crealeiphi(phi,z[n-1]) > crealeiphi(phi,z[ex->ind_min[ex->n_min]])) { /* the edge is not a min */
        #endif

            if (2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_max[ex->n_max-1-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-1-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-1-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
                }
            }
        } else { /* edge is a min -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_min[ex->n_min+2+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                ex->ry_min[ex->n_min+2+cour] = ex->ry_min[ex->n_min-cour];
                ex->iy_min[ex->n_min+2+cour] = ex->iy_min[ex->n_min-cour];
            }
            ex->x_min[ex->n_min+1] = x[n-1];
            ex->ry_min[ex->n_min+1] = CREAL(z[n-1]);
            ex->iy_min[ex->n_min+1] = CIMAG(z[n-1]);
        }
    }
    
    
    
    (ex->n_min) = ex->n_min + nbsym + 1;
    (ex->n_max) = ex->n_max + nbsym + 1;
}


/************************************************************************/
/*                                                                      */
/* FREE ALLOCATED MEMORY                                                */
/*                                                                      */
/************************************************************************/

void free_extr(extrema_t ex) 
/*<free allocated extrema struct>*/
{
    free(ex.x_max);
    free(ex.x_min);
    free(ex.ind_max);
    free(ex.ind_min);
    free(ex.ry_max);
    free(ex.ry_min);
    free(ex.iy_max);
    free(ex.iy_min);
}

//interpolation.c
/*************************************************************************/
/*                                                                       */
/* INTERPOLATION                                                         */
/*                                                                       */
/* interpolates the sequence (xs,ys) at instants in x using cubic spline */
/*                                                                       */
/*************************************************************************/

void interpolation(double y[],double xs[],double ys[],int n,double x[], int nx,double *ys2, double *temp) 
/*< interpolates the sequence (xs,ys) at instants in x using cubic spline >*/
{
  int i,j,jfin,cur,prev;
  double p,sig,a,b,c,d,e,f,g,a0,a1,a2,a3,xc;

  /* Compute second derivatives at the knots */
  ys2[0]=temp[0]=0.0;
  for (i=1;i<n-1;i++) {
    sig=(xs[i]-xs[i-1])/(xs[i+1]-xs[i-1]);
    p=sig*ys2[i-1]+2.0;
    ys2[i]=(sig-1.0)/p;
    temp[i]=(ys[i+1]-ys[i])/(xs[i+1]-xs[i])-(ys[i]-ys[i-1])/(xs[i]-xs[i-1]);
    temp[i]=(6.0*temp[i]/(xs[i+1]-xs[i-1])-sig*temp[i-1])/p;
  }
  ys2[n-1]=0.0;
  for (j=n-2;j>=0;j--) ys2[j]=ys2[j]*ys2[j+1]+temp[j];

  /* Compute the spline coefficients */
  cur=0;
  j=0;
  jfin=n-1;
  while (xs[j+2]<x[0]) j++;
  while (xs[jfin]>x[nx-1]) jfin--;
  for (;j<=jfin;j++) {
    /* Compute the coefficients of the polynomial between two knots */
    a=xs[j];
    b=xs[j+1];
    c=b-a;
    d=ys[j];
    e=ys[j+1];
    f=ys2[j];
    g=ys2[j+1];
    a0=(b*d-a*e+CUBE(b)*f/6-CUBE(a)*g/6)/c+c*(a*g-b*f)/6;
    a1=(e-d-SQUARE(b)*f/2+SQUARE(a)*g/2)/c+c*(f-g)/6;
    a2=(b*f-a*g)/(2*c);
    a3=(g-f)/(6*c);


    prev=cur;
    while ((cur<nx) && ((j==jfin) || (x[cur]<xs[j+1]))) cur++;

    /* Compute the value of the spline at the sampling times x[i] */
    for (i=prev;i<cur;i++) {
      xc=x[i];
      y[i]=a0+a1*xc+a2*SQUARE(xc)+a3*CUBE(xc);
    }
  }
}

//cio.c
/************************************************************************/
/*                                                                      */
/* INITIALIZATION OF THE LIST                                           */
/*                                                                      */
/************************************************************************/

imf_list_t init_imf_list(int n) 
/*< define imf list struct >*/
{
  imf_list_t list;
  list.first=NULL;
  list.last=NULL;
  list.n=n;
  list.m=0;
  return list;
}


/************************************************************************/
/*                                                                      */
/* ADD AN IMF TO THE LIST                                               */
/*                                                                      */
/************************************************************************/

void add_imf(imf_list_t *list,COMPLEX_T *p,int nb_it) 
/*< adding imf to imf list >*/
{
  COMPLEX_T *v=(COMPLEX_T *)malloc(list->n*sizeof(COMPLEX_T));
  int i;
  imf_t *mode=(imf_t *)malloc(sizeof(imf_t));
  for (i=0;i<list->n;i++) v[i]=p[i];
  mode->pointer=v;
  mode->nb_iterations=nb_it;
  mode->next=NULL;
  if (!list->first) {
    list->first=mode;
  } else {
    (list->last)->next=mode;
  }
  list->last=mode;
  list->m++;
}


/************************************************************************/
/*                                                                      */
/* FREE MEMORY ALLOCATED FOR THE LIST                                   */
/*                                                                      */
/************************************************************************/

void free_imf_list(imf_list_t list) 
/*<free allocated imf list struct >*/
{
  imf_t *current=list.first, *previous;
  while (current) {
    previous=current;
    current=current->next;
    free(previous->pointer);
    free(previous);
  }
}

//clocal_mean.c
/********************************************************/
/* ALLOCATE MEMORY FOR THE ENVELOPES AND TEMPORARY DATA */
/********************************************************/

envelope_t init_local_mean(int n) 
/*< initialization for local mean of bivariate emd >*/
{
  envelope_t env;
  env.re_min = (double*)malloc(n*sizeof(double));
  env.ie_min = (double*)malloc(n*sizeof(double));
  env.re_max = (double*)malloc(n*sizeof(double));
  env.ie_max = (double*)malloc(n*sizeof(double));
  env.tmp1 = (double*)malloc(n*sizeof(double));
  env.tmp2 = (double*)malloc(n*sizeof(double));
  return env;
}

/*************************/
/* FREE ALLOCATED MEMORY */
/*************************/

void free_local_mean(envelope_t env) 
/*<free allocated local mean struct >*/
{
  free(env.re_min);
  free(env.ie_min);
  free(env.re_max);
  free(env.ie_max);
  free(env.tmp1);
  free(env.tmp2);
}

/***************************************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES AND THE AMPLITUDE OF THE CURRENT IMF */
/***************************************************************************/

int mean_and_amplitude(double *x,COMPLEX_T *z,COMPLEX_T *m,double *a,int n,int nbphases,extrema_t *ex,envelope_t *env) 
/*< compute the mean of the envelopes and the amplitude of the current imf >*/
{
  int i,k;
  double phi;
  
  #ifdef C99_OK
  for (i=0;i<n;i++) m[i]=0;
  #else
  for (i=0;i<n;i++) {
    m[i].r = 0;
    m[i].i = 0;
  }
  #endif
  for (i=0;i<n;i++) a[i]=0;
  
  for(k=0;k<nbphases;k++) {
    
    phi = k*M_PI/nbphases;
    /* detect maxima and minima in direction phi=k*M_PI/nbphases*/
    extr(x,z,phi,n,ex);
    if (ex->n_max+ex->n_min <3){ /* not enough extrema in a direction -> stop */
      return 1;
    }

    /* add extra points at the edges */
    boundary_conditions(x,z,phi,n,ex);
    
    /* interpolation - upper envelope */
    interpolation(env->re_max,ex->x_max,ex->ry_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    interpolation(env->ie_max,ex->x_max,ex->iy_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    
    /* interpolation - lower envelope */
    interpolation(env->re_min,ex->x_min,ex->ry_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    interpolation(env->ie_min,ex->x_min,ex->iy_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    if ((ex->n_min > LIM_GMP)||(ex->n_min > LIM_GMP)) {
      sf_warning("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* compute the mean and amplitude*/
    #ifdef C99_OK
    for (i=0;i<n;i++) m[i]+=(env->re_max[i]+env->re_min[i]+I*(env->ie_max[i]+env->ie_min[i]))/(2*nbphases);
    for (i=0;i<n;i++) a[i]+=CABS(env->re_max[i]-env->re_min[i]+I*(env->ie_max[i]-env->ie_min[i]))/(2*nbphases);
    #else
    for (i=0;i<n;i++) {
      m[i].r+=(env->re_max[i]+env->re_min[i])/(2*nbphases);
      m[i].i+=(env->ie_max[i]+env->ie_min[i])/(2*nbphases);
    }
    for (i=0;i<n;i++) a[i]+=sqrt((env->re_max[i]-env->re_min[i])*(env->re_max[i]-env->re_min[i])+(env->ie_max[i]-env->ie_min[i])*(env->ie_max[i]-env->ie_min[i]))/(2*nbphases);
    #endif
  }
 
  return 0;
}

/*********************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES OF THE CURRENT IMF */
/*********************************************************/

int mean(double *x,COMPLEX_T *z,COMPLEX_T *m,int n,int nbphases,extrema_t *ex,envelope_t *env) 
/*< compute the mean of the envelopes of the current imf >*/
{
  int i,k;
  double phi;
  
  #ifdef C99_OK
  for (i=0;i<n;i++) m[i]=0;
  #else
  for (i=0;i<n;i++) {
    m[i].r = 0;
    m[i].i = 0;
  }
  #endif
  
  for(k=0;k<nbphases;k++) {
    
    phi = k*M_PI/nbphases;
    /* detect maxima and minima in direction phi=k*M_PI/nbphases*/
    extr(x,z,phi,n,ex);
    if (ex->n_max+ex->n_min <3){ /* not enough extrema in a direction -> stop */
      return 1;
    }
    
    boundary_conditions(x,z,phi,n,ex);
    
    /* interpolation - upper envelope */
    if (ex->n_max < LIM_GMP) {
      interpolation(env->re_max,ex->x_max,ex->ry_max,ex->n_max,x,n,env->tmp1,env->tmp2);
      interpolation(env->ie_max,ex->x_max,ex->iy_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    }
    else {
      sf_warning("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* interpolation - lower envelope */
    if (ex->n_min < LIM_GMP) {
      interpolation(env->re_min,ex->x_min,ex->ry_min,ex->n_min,x,n,env->tmp1,env->tmp2);
      interpolation(env->ie_min,ex->x_min,ex->iy_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    }
    else {
      sf_warning("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* compute the mean*/
    #ifdef C99_OK
    for (i=0;i<n;i++) m[i]+=(env->re_max[i]+env->re_min[i]+I*(env->ie_max[i]+env->ie_min[i]))/(2*nbphases);
    #else
    for (i=0;i<n;i++) {
      m[i].r+=(env->re_max[i]+env->re_min[i])/(2*nbphases);
      m[i].i+=(env->ie_max[i]+env->ie_min[i])/(2*nbphases);
    }
    #endif
    
  }
  return 0;
}

//emd_complex.c

double CREAL(complex_data_t c)
/*< creal >*/
{
  return c.r;
}

double CIMAG(complex_data_t c)
/*< cimag >*/
{
  return c.i;
}

double CABS(complex_data_t c)
/*< cabs >*/
{  
  return sqrt((c.r)*(c.r)+(c.i)*(c.i));
}

double crealeiphi(double phi,complex_data_t c)
/*< crealeiphi >*/
{
  return cos(phi)*c.r-sin(phi)*c.i;
}

//cemdc.c

/************************************************************************/
/* STOP TEST FOR THE SIFTING LOOP                                    */
/************************************************************************/

int stop_sifting(COMPLEX_T *m, double *a,extrema_t *ex,stop_t *sp,int n, int counter, int max_iterations) 
/*< decide if stop sifting >*/
{
  int i,count;
  double tol,eps;
  tol = sp->tolerance*n;
  eps = sp->threshold;
  count = 0;
  if (counter >= MAX_ITERATIONS) return 1;
  for (i=0;i<n;i++) {
    if (CABS(m[i]) > eps*a[i]) if (++count>tol) return 0;
  }
  return 1;
}
