#include <rsf.h>
#include "emdutil.h"

#define DEFAULT_THRESHOLD 0.05
#define DEFAULT_TOLERANCE 0.05
#define MAX_ITERATIONS 1000
#define LIM_GMP 30000
#define NBSYM 2
#define SQUARE(A) (A*A)
#define CUBE(A) (A*A*A)
#define GREATER(A,B) ((A) >= (B) ? (A) : (B))
#define SMALLER(A,B) ((A) <  (B) ? (A) : (B))
/*^*/

#ifndef _emdutil_h

#ifndef EXTR_H
#define EXTR_H

/* structure used to store the local extrema of the signal */
typedef struct {
  int n_min;
  int n_max;
  double *x_min;
  double *y_min;
  double *x_max;
  double *y_max;
} extrema_t;
/*^*/
#endif

#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#endif

#ifndef EMD_IO_H
#define EMD_IO_H

typedef struct {
    double threshold;
    double tolerance;
} stop_t;
/*^*/

/* structure used to store an IMF and the associated number of iterations */
typedef struct i {
  int nb_iterations;
  double *pointer;
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

#ifndef LOCAL_MEAN_H
#define LOCAL_MEAN_H

/* structure used to store envelopes and temporary data */
typedef struct {
  int n;
  double *e_min;
  double *e_max;
  double *tmp1;
  double *tmp2;
} envelope_t;
/*^*/
#endif
#endif

/************************************************************************/
/*                                                                      */
/* INITIALIZATION OF EXTREMA STRUCTURE                                  */
/*                                                                      */
/************************************************************************/

extrema_t init_extr(int n) 
/*<initialization for extremas extraction>*/
{
    extrema_t ex;
    ex.x_min=(double *)sf_alloc(n,sizeof(double));
    ex.x_max=(double *)sf_alloc(n,sizeof(double));
    ex.y_min=(double *)sf_alloc(n,sizeof(double));
    ex.y_max=(double *)sf_alloc(n,sizeof(double));
    ex.n_min=n;
    ex.n_max=n;
    return ex;
}


/************************************************************************/
/*                                                                      */
/* DETECTION OF LOCAL EXTREMA                                           */
/*                                                                      */
/************************************************************************/

void extr(double x[],double y[],int n,extrema_t *ex)
/*< extract extremas >*/
 {
    int cour;
    ex->n_min=0;
    ex->n_max=0;
    
  /* search for extrema */
    for(cour=1;cour<(n-1);cour++) {
        if (y[cour]<=y[cour-1] && y[cour]<=y[cour+1]) /* local minimum */ {
            ex->x_min[ex->n_min+NBSYM]=x[cour];
            ex->y_min[ex->n_min+NBSYM]=y[cour];
            ex->n_min++;
        }
        if (y[cour]>=y[cour-1] && y[cour]>=y[cour+1]) /* local maximum */ {
            ex->x_max[ex->n_max+NBSYM]=x[cour];
            ex->y_max[ex->n_max+NBSYM]=y[cour];
            ex->n_max++;
        }
    }
}

/************************************************************************/
/*                                                                      */
/* EXTRAPOLATION OF EXTREMA TO LIMIT BORDER EFFECTS                     */
/*                                                                      */
/************************************************************************/

void boundary_conditions(double x[],double y[],int n,extrema_t *ex) 
/*< setting boundary conditions >*/
{
    int cour,nbsym;
    
    nbsym = NBSYM;
    /* reduce the number of symmetrized points if there is not enough extrema */
    while(ex->n_min < nbsym+1 && ex->n_max < nbsym+1) nbsym--;
    if (nbsym < NBSYM) {
        for(cour=0;cour<ex->n_max;cour++) {
            ex->x_max[nbsym+cour] = ex->x_max[NBSYM+cour];
            ex->y_max[nbsym+cour] = ex->y_max[NBSYM+cour];
        }
        for(cour=0;cour<ex->n_min;cour++) {
            ex->x_min[nbsym+cour] = ex->x_min[NBSYM+cour];
            ex->y_min[nbsym+cour] = ex->y_min[NBSYM+cour];
        }
    }
    
    /* select the symmetrized points and the axis of symmetry at the beginning of the signal*/
    if (ex->x_max[nbsym] < ex->x_min[nbsym]) { /* first = max */
        if (y[0] > ex->y_min[nbsym]) { /* the edge is not a min */
            if (2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                    ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                    ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*ex->x_max[nbsym]-ex->x_max[2*nbsym-cour];
                    ex->y_max[cour] = ex->y_max[2*nbsym-cour];
                    ex->x_min[cour] = 2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1-cour];
                    ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
                }
            }
        } else { /* edge is a min -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-2-cour];
                ex->y_min[cour] = ex->y_min[2*nbsym-2-cour];
            }
            ex->x_min[nbsym-1] = x[0];
            ex->y_min[nbsym-1] = y[0];
        }
    } else { /* first = min */
        
        if (y[0] < ex->y_max[nbsym]) { /* the edge is not a max */
            if (2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                    ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                    ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1-cour];
                    ex->y_max[cour] = ex->y_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*ex->x_min[nbsym]-ex->x_min[2*nbsym-cour];
                    ex->y_min[cour] = ex->y_min[2*nbsym-cour];
                }
            }
        } else { /* edge is a max -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                ex->y_min[cour] = ex->y_min[2*nbsym-1-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-2-cour];
                ex->y_max[cour] = ex->y_max[2*nbsym-2-cour];
            }
            ex->x_max[nbsym-1] = x[0];
            ex->y_max[nbsym-1] = y[0];
        }
    }
    
    
    (ex->n_min) += nbsym-1;
    (ex->n_max) += nbsym-1;
    
    /* select the symmetrized points and the axis of symmetry at the end of the signal*/
    if (ex->x_max[ex->n_max] < ex->x_min[ex->n_min]) { /* last is a min */
        if (y[n-1] < ex->y_max[ex->n_max]) { /* the edge is not a max */
            if (2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                    ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                    ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-cour];
                    ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_min[ex->n_min-1-cour];
                    ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-1-cour];
                }
            }
        } else { /* edge is a max -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_max[ex->n_max+2+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                ex->y_max[ex->n_max+2+cour] = ex->y_max[ex->n_max-cour];
            }
            ex->x_max[ex->n_max+1] = x[n-1];
            ex->y_max[ex->n_max+1] = y[n-1];
        }
    } else {  /* last is a max */
        if (y[n-1] > ex->y_min[ex->n_min]) { /* the edge is not a min */
            if (2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                    ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                    ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_max[ex->n_max-1-cour];
                    ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-1-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-cour];
                    ex->y_min[ex->n_min+1+cour] = ex->y_min[ex->n_min-cour];
                }
            }
        } else { /* edge is a min -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                ex->y_max[ex->n_max+1+cour] = ex->y_max[ex->n_max-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_min[ex->n_min+2+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                ex->y_min[ex->n_min+2+cour] = ex->y_min[ex->n_min-cour];
            }
            ex->x_min[ex->n_min+1] = x[n-1];
            ex->y_min[ex->n_min+1] = y[n-1];
        }
    }
    
    (ex->n_min) = ex->n_min + nbsym + 1;
    (ex->n_max) = ex->n_max + nbsym + 1;
}

void free_extr(extrema_t ex) 
/*<free allocated extrema struct>*/
{
    free(ex.x_max);
    free(ex.x_min);
    free(ex.y_max);
    free(ex.y_min);
}

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

void add_imf(imf_list_t *list,double *p,int nb_it) 
/*< adding imf to imf list >*/
{
    double *v=(double *)sf_alloc(list->n,sizeof(double));
  int i;
  imf_t *mode=(imf_t *)sf_alloc(1,sizeof(imf_t));
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
/*< initialization for local mean of bivariate emd >*/
{
  imf_t *current=list.first, *previous;
  while (current) {
    previous=current;
    current=current->next;
    free(previous->pointer);
    free(previous);
  }
}

/********************************************************/
/* ALLOCATE MEMORY FOR THE ENVELOPES AND TEMPORARY DATA */
/********************************************************/

envelope_t init_local_mean(int n) 
/*<free allocated local mean struct >*/
{
  envelope_t env;
  env.e_min = (double*)sf_alloc(n,sizeof(double));
  env.e_max = (double*)sf_alloc(n,sizeof(double));
  env.tmp1 = (double*)sf_alloc(n,sizeof(double));
  env.tmp2 = (double*)sf_alloc(n,sizeof(double));
  env.n = n;
  return env;
}

/*************************/
/* FREE ALLOCATED MEMORY */
/*************************/

void free_local_mean(envelope_t env) 
/*< free_local_mean >*/
{
  free(env.e_min);
  free(env.e_max);
  free(env.tmp1);
  free(env.tmp2);
}

/*********************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES OF THE CURRENT IMF */
/*********************************************************/

int mean(double *x,double *z,double *m,int n,extrema_t *ex,envelope_t *env) 
/*< compute the mean of the envelopes and the amplitude of the current imf >*/
{
  int i;
  /* detect maxima and minima */
  extr(x,z,n,ex);
  /* if not enough extrema -> stop */
  if (ex->n_min+ex->n_max <7)
    return 1;
  /* add extra points at the edges */
  boundary_conditions(x,z,n,ex);
  /* interpolation - upper envelope */
  interpolation(env->e_max,ex->x_max,ex->y_max,ex->n_max,x,n,env->tmp1,env->tmp2);
  /* interpolation - lower envelope */
  interpolation(env->e_min,ex->x_min,ex->y_min,ex->n_min,x,n,env->tmp1,env->tmp2);
  if ((ex->n_min > LIM_GMP)||(ex->n_min > LIM_GMP)) {
    sf_warning("Too many extrema, interpolation may be unreliable\n");
  }
  /* compute the mean */
  for (i=0;i<n;i++) m[i]=(env->e_max[i]+env->e_min[i])/2;
  return 0;
}

/***************************************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES AND THE AMPLITUDE OF THE CURRENT IMF */
/***************************************************************************/

int mean_and_amplitude(double *x,double *z,double *m,double *a,int n,extrema_t *ex,envelope_t *env) 
/*< compute the mean of the envelopes of the current imf >*/
{
  int i;
  /* detect maxima and minima */
  extr(x,z,n,ex);
  /* if not enough extrema -> stop */
  if (ex->n_min+ex->n_max <7)
    return 1;
  /* add extra points at the edges */
  boundary_conditions(x,z,n,ex);
  /* interpolation - upper envelope */
  interpolation(env->e_max,ex->x_max,ex->y_max,ex->n_max,x,n,env->tmp1,env->tmp2);
  /* interpolation - lower envelope */
  interpolation(env->e_min,ex->x_min,ex->y_min,ex->n_min,x,n,env->tmp1,env->tmp2);
  /* compute the mean */
  for (i=0;i<n;i++) m[i]=(env->e_max[i]+env->e_min[i])/2;
  /* compute the amplitude */
  for (i=0;i<n;i++) a[i]=(env->e_max[i]-env->e_min[i])/2;
  return 0;
}

/************************************************************************/
/* ABSOLUTE VALUE                                                       */
/************************************************************************/

double emd_fabs(double x) 
/*< absolute value  >*/
{
  if (x <0) return -x;
  else return x;
}


/************************************************************************/
/* STOP TEST FOR THE SIFTING LOOP                                       */
/************************************************************************/

int stop_sifting(double *m, double *a,extrema_t *ex,stop_t *sp,int n, int counter, int max_iterations)
/*< decide if stop sifting >*/
{
  int i,count;
  double tol,eps;
  tol = sp->tolerance*n;
  eps = sp->threshold;
  count = 0;
  if (counter >= MAX_ITERATIONS) return 1;
  for (i=0;i<ex->n_min;i++) if (ex->y_min[i] > 0) return 0;
  for (i=0;i<ex->n_max;i++) if (ex->y_max[i] < 0) return 0;
  for (i=0;i<n;i++) {
    if (emd_fabs(m[i]) > eps*emd_fabs(a[i])) if (++count>tol) return 0;
  }
  return 1;
}
