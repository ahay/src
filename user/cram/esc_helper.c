/* Auxiliary functions for escape solvers */
/*
  Copyright (C) 2012 University of Texas at Austin

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

/* ================================= Quaternions ================================================ */

void sf_quat_rotmat (float *q, float *M)
/*< Build rotation matrix[3][3] from quaternion q[4] >*/
{
    M[0] = 1.0 - 2.0*q[2]*q[2] - 2.0*q[3]*q[3];
    M[1] = 2.0*q[1]*q[2] - 2.0*q[3]*q[0];
    M[2] = 2.0*q[1]*q[3] + 2.0*q[2]*q[0];
    M[3] = 2.0*q[1]*q[2] + 2.0*q[3]*q[0];
    M[4] = 1.0 - 2.0*q[1]*q[1] - 2.0*q[3]*q[3];
    M[5] = 2.0*q[3]*q[2] - 2.0*q[1]*q[0];
    M[6] = 2.0*q[1]*q[3] - 2.0*q[2]*q[0];
    M[7] = 2.0*q[3]*q[2] + 2.0*q[1]*q[0];
    M[8] = 1.0 - 2.0*q[1]*q[1] - 2.0*q[2]*q[2];
}

void sf_quat_norm (float *q)
/*< Normalize quaternion q[4] >*/
{
    float f = sqrtf (q[0]*q[0] + q[1]*q[1] + 
               q[2]*q[2] + q[3]*q[3]);
    q[0] /= f;
    q[1] /= f;
    q[2] /= f;
    q[3] /= f;
    return;
}

void sf_quat_vecrot (float *vf, float *vt, float *q)
/*< Find quaternion q[4] from rotation of vector vf[3] to vt[3] >*/
{
    float s, invs, d, c[3];

    /* Dot product */
    d = vf[0]*vt[0] + vf[1]*vt[1] + vf[2]*vt[2];
    
    if (d >= 1.0) { /* Same vectors */
        q[0] = 1.0;
        q[1] = 0.0;
        q[2] = 0.0;
        q[3] = 0.0;
        return;
    } else if (d <= -1.0) { /* Opposite vectors */
        /* Cross product with the vertical direction */
        c[0] = 0.0;
        c[1] = vf[2];
        c[2] = -vf[1];
        q[0] = 0.0;
    } else {
        s = sqrtf ((1.0 + d)*2.0);
        invs = 1.0/s;
        c[0] = (vf[1]*vt[2] - vf[2]*vt[1])*invs;
        c[1] = (vf[2]*vt[0] - vf[0]*vt[2])*invs;
        c[2] = (vf[0]*vt[1] - vf[1]*vt[0])*invs;
        q[0] = 0.5*s;
    }
    q[1] = c[0];
    q[2] = c[1];
    q[3] = c[2];
    /* Normalize */
    sf_quat_norm (q);
    return;
}

/* ================================= LU matrix decomposition ================================================ */

/* LU decomposition with Doolittle method, original code can
   be found at http://mymathlib.com */

int sf_ludlt_decomposition (float *A, int *pivot, int n)
/*< Find LU decomposition of matrix A[nxn], return
    -1 for singular matrix, otherwise - 0 >*/
{
    int i, j, k;
    float *p_k, *p_row, *p_col = NULL;
    float max;

    /* For each row and column, k = 0, ..., n-1, */
    for (k = 0, p_k = A; k < n; p_k += n, k++) {
        /* Find the pivot row */
        pivot[k] = k;
        max = fabsf (*(p_k + k));
        for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n) {
            if (max < fabsf (*(p_row + k))) {
                max = fabsf (*(p_row + k));
                pivot[k] = j;
                p_col = p_row;
            }
        }
        /* and if the pivot row differs from the current row, then
           interchange the two rows.*/
        if (pivot[k] != k) {
            for (j = 0; j < n; j++) {
                max = *(p_k + j);
                *(p_k + j) = *(p_col + j);
                *(p_col + j) = max;
            }
        }
        /* and if the matrix is singular, return error */
        if (*(p_k + k) == 0.0)
            return -1;
        /* otherwise find the lower triangular matrix
           elements for column k.*/
        for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++) {
            *(p_row + k) /= *(p_k + k);
        }
        /* update remaining matrix */
        for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++) {
            for (j = k+1; j < n; j++) {
                *(p_row + j) -= *(p_row + k) * *(p_k + j);
            }
        }
    }
    return 0;
}

/* Matrix vector multiplication for a LU-decomposed system with Doolittle method,
   original code can be found at http://mymathlib.com */

int sf_ludtl_solve (float *A, float *b, int *pivot, float *x, int n)
/*< Find vecotr x[n] for RHS b[n] and prveiously
    LU-processed matrix A[nxn], return -1 for improper
    matrix, otherwise - 0 >*/
{
    int i, k;
    float *p_k;
    float dum;
    
    /* Solve the linear equation Lx = B for x, where L is a lower
       triangular matrix with an implied 1 along the diagonal.*/
    for (k = 0, p_k = A; k < n; p_k += n, k++) {
        if (pivot[k] != k) {
            dum = b[k];
            b[k] = b[pivot[k]];
            b[pivot[k]] = dum;
        }
        x[k] = b[k];
        for (i = 0; i < k; i++) {
            x[k] -= x[i] * *(p_k + i);
        }
    }
    /* Solve the linear equation Ux = y, where y is the solution
       obtained above of Lx = B and U is an upper triangular matrix. */
    for (k = n-1, p_k = A + n*(n-1); k >= 0; k--, p_k -= n) {
        if (pivot[k] != k) {
            dum = b[k];
            b[k] = b[pivot[k]];
            b[pivot[k]] = dum;
        }
        for (i = k + 1; i < n; i++) {
            x[k] -= x[i] * *(p_k + i);
        }
        if (*(p_k + k) == 0.0)
            return -1;
        x[k] /= *(p_k + k);
    }
    return 0;
}

/* Thin plate spline function */
static double sf_tps_base_func (double r) {
    if (r == 0.0)
        return 0.0;
    else
        return r*r*log (r);
}

/* ================================= TPS interpolation ================================================ */

/* Thin plate spline interpolation, original code by Jarno Elonen */

void sf_tps_build_matrix (float **L, float *x, float *y, int nd,
                          float eps)
/*< Build thin plate spline matrix L[nd+3][nd+3] for coordinate
    locations x[nd] and y[nd]; eps is bending energy >*/
{
    int i, j;
    double a = 0.0, elen;

    if (nd < 3)
        return;

    /* Fill K (p x p, upper left of L) and calculate
       mean edge length from control points.
       K is symmetrical, thus only half of the 
       coefficients has to be calculated. */
    for (i = 0; i < nd; i++) {
        for (j = i + 1; j < nd; j++) {
            elen = sqrt ((x[i] - x[j])*(x[i] - x[j]) + 
                         (y[i] - y[j])*(y[i] - y[j]));
            L[i][j] = sf_tps_base_func (elen);
            L[j][i] = L[i][j];
            a += elen*2;
        }
    }
    a /= (double)(nd*nd);

    /* Fill the rest of L */
    for (i = 0; i < nd; i++) {
        /* Diagonal: reqularization parameters (lambda * a^2) */
        L[i][i] = eps*(a*a);
        /* P (p x 3, upper right) */
        L[i][nd] = 1.0;
        L[i][nd + 1] = x[i];
        L[i][nd + 2] = y[i];
        /* P transposed (3 x p, bottom left) */
        L[nd][i] = 1.0;
        L[nd + 1][i] = x[i];
        L[nd + 2][i] = y[i];
    }
    /* O (3 x 3, lower right) */
    for (i = nd; i < nd + 3; i++) {
        for (j = nd; j < nd + 3; j++) {
            L[i][j] = 0.0;
        }
    }
}

float sf_tps_compute_point (float *w, float *x, float *y, int nd,
                            float xc, float yc)
/*< Compute point at xc,yc using thin-plate spline coefficients w[nd+3]
    and control nodes at x[nd], y[nd] >*/
{
    int i;
    float h;

    /* Do interpolation */
    h = w[nd] + w[nd + 1]*xc + w[nd + 2]*yc;
    for (i = 0; i < nd; i++) {
        h += w[i]*sf_tps_base_func (sqrtf ((x[i] - xc)*(x[i] - xc) + 
                                           (y[i] - yc)*(y[i] - yc)));
            }
    return h;
}

#ifdef PTHREADS

/* ================================= Pthreads work queue ================================================ */
               
#include <pthread.h>
#include <semaphore.h>
/*^*/

/* Original code by Johan Hanssen Seferidis, Licence: LGPL */

/* Individual job */
typedef struct thpool_job_t {
        void*  (*function)(void* arg); /* function pointer         */
        void*                     arg; /* function's argument      */
        struct thpool_job_t*     next; /* pointer to next job      */
        struct thpool_job_t*     prev; /* pointer to previous job  */
} thpool_job_t;

/* Job queue as doubly linked list */
typedef struct thpool_jobqueue {
        thpool_job_t    *head;      /* pointer to head of queue */
        thpool_job_t    *tail;      /* pointer to tail of queue */
        int              jobsN;     /* amount of jobs in queue  */
        sem_t           *queueSem;  /* semaphore(this is probably just holding the same as jobsN) */
} thpool_jobqueue;

typedef struct thpool_t *sf_thpool;
/* abstract data type */
/*^*/

/* The thread pool */
typedef struct thpool_t {
        pthread_t*       threads;   /* pointer to threads' ID   */
        int              threadsN;  /* amount of threads        */
        thpool_jobqueue* jobqueue;  /* pointer to the job queue */
        pthread_mutex_t  mutex;     /* used to serialize queue access */
        bool             alive;
} thpool_t;

/* Remove job from queue */
static bool thpool_jobqueue_removelast (thpool_t* tp_p) {
    thpool_job_t *oldLastJob;
    oldLastJob = tp_p->jobqueue->tail;

    /* fix jobs' pointers */
    switch(tp_p->jobqueue->jobsN){        
        case 0:     /* if there are no jobs in queue */
            return false;
            break;
        case 1:     /* if there is only one job in queue */
            tp_p->jobqueue->tail = NULL;
            tp_p->jobqueue->head = NULL;
            break;
        default:     /* if there are more than one jobs in queue */
            oldLastJob->prev->next = NULL;
            tp_p->jobqueue->tail = oldLastJob->prev;
    }

    (tp_p->jobqueue->jobsN)--;

    int sval;
    sem_getvalue (tp_p->jobqueue->queueSem, &sval);
    return true;
}

/* Get first element from queue */
static thpool_job_t* thpool_jobqueue_peek (thpool_t* tp_p) {
    return tp_p->jobqueue->tail;
}

/* What each individual thread is doing */
static void* thpool_thread_do (void *ud) {
    thpool_t* tp_p = (thpool_t*)ud;

    while (tp_p->alive) {
        if (sem_wait (tp_p->jobqueue->queueSem)) { /* Waiting until there is work in the queue */
            perror ("thpool_thread_do(): Waiting for semaphore");
            exit (-1);
        }

        /* Read job from queue and execute it */
        void*(*func_buff)(void* arg);
        void* arg_buff;
        thpool_job_t* job_p;

        pthread_mutex_lock (&tp_p->mutex);

        job_p = thpool_jobqueue_peek (tp_p);
        func_buff = job_p->function;
        arg_buff  = job_p->arg;
        thpool_jobqueue_removelast (tp_p);

        pthread_mutex_unlock (&tp_p->mutex);

        func_buff (arg_buff); /* run job */
    }
    return NULL;
}

/* Initialize queue */
static int thpool_jobqueue_init (thpool_t* tp_p) {
    tp_p->jobqueue = (thpool_jobqueue*)malloc (sizeof(thpool_jobqueue));
    if (tp_p->jobqueue == NULL)
        return -1;
    tp_p->jobqueue->tail = NULL;
    tp_p->jobqueue->head = NULL;
    tp_p->jobqueue->jobsN = 0;
    return 0;
}

sf_thpool sf_thpool_init (int threadsN)
/*< Initialize thread pool >*/
{
    thpool_t* tp_p;

    if (!threadsN || threadsN < 1)
        threadsN = 1;

    /* Make new thread pool */
    tp_p = (thpool_t*)malloc (sizeof(thpool_t));
    if (tp_p == NULL) {
        fprintf (stderr, "thpool_init(): Could not allocate memory for thread pool\n");
        return NULL;
    }
    tp_p->threads = (pthread_t*)malloc (threadsN*sizeof(pthread_t));
    if (tp_p->threads == NULL) {
        fprintf (stderr, "thpool_init(): Could not allocate memory for thread IDs\n");
        return NULL;
    }
    tp_p->threadsN = threadsN;

    /* Initialise the job queue */
    if (thpool_jobqueue_init (tp_p) == -1) {
        fprintf (stderr, "thpool_init(): Could not allocate memory for job queue\n");
        return NULL;
    }

    /* Initialise semaphore */
    tp_p->jobqueue->queueSem = (sem_t*)malloc (sizeof(sem_t));
    sem_init (tp_p->jobqueue->queueSem, 0, 0); /* no shared, initial value */

    if (pthread_mutex_init (&tp_p->mutex, NULL) != 0) {
        fprintf (stderr, "sem_init(): Could not create a semaphore\n");
        return NULL;
    }
    tp_p->alive = true;

    /* Make threads in pool */
    int t;
    for (t = 0; t < threadsN; t++) {
        pthread_create (&(tp_p->threads[t]), NULL, thpool_thread_do,
                        (void *)tp_p);
    }

    return tp_p;
}

/* Remove and deallocate all jobs in queue */
static void thpool_jobqueue_empty (thpool_t* tp_p) {
    
    thpool_job_t* curjob;
    curjob=tp_p->jobqueue->tail;

    while (tp_p->jobqueue->jobsN) {
        tp_p->jobqueue->tail = curjob->prev;
        free (curjob);
        curjob = tp_p->jobqueue->tail;
        tp_p->jobqueue->jobsN--;
    }

    /* Fix head and tail */
    tp_p->jobqueue->tail = NULL;
    tp_p->jobqueue->head = NULL;
}

int sf_thpool_jobsn (sf_thpool tp_p)
/*< Return number of jobs in the queue >*/
{
    return tp_p->jobqueue->jobsN;
}

void sf_thpool_destroy (sf_thpool tp_p)
/*< Destroy the thread pool >*/
{
    int t;

    /* End each thread's infinite loop */
    tp_p->alive = false; 

    /* Awake idle threads waiting at semaphore */
    for (t = 0; t < (tp_p->threadsN); t++){
        if (sem_post (tp_p->jobqueue->queueSem)) {
            fprintf (stderr, "thpool_destroy(): Could not bypass sem_wait()\n");
        }
    }

    /* Kill semaphore */
    if (sem_destroy (tp_p->jobqueue->queueSem) != 0) {
        fprintf (stderr, "thpool_destroy(): Could not destroy semaphore\n");
    }

    /* Wait for threads to finish */
    for (t = 0; t < (tp_p->threadsN); t++) {
        pthread_join (tp_p->threads[t], NULL);
    }

    thpool_jobqueue_empty (tp_p);

    pthread_mutex_destroy (&tp_p->mutex);
    free (tp_p->threads);
    free (tp_p->jobqueue->queueSem);
    free (tp_p->jobqueue);
    free (tp_p);
}

/* Add job to queue */
static void thpool_jobqueue_add (thpool_t* tp_p, thpool_job_t* newjob_p) {
    newjob_p->next = NULL;
    newjob_p->prev = NULL;

    thpool_job_t *oldFirstJob;
    oldFirstJob = tp_p->jobqueue->head;

    /* fix jobs' pointers */
    switch(tp_p->jobqueue->jobsN) {
        case 0:     /* if there are no jobs in queue */
            tp_p->jobqueue->tail = newjob_p;
            tp_p->jobqueue->head = newjob_p;
            break;
        default:     /* if there are already jobs in queue */
            oldFirstJob->prev = newjob_p;
            newjob_p->next = oldFirstJob;
            tp_p->jobqueue->head = newjob_p;
    }

    (tp_p->jobqueue->jobsN)++;     /* increment amount of jobs in queue */
    sem_post (tp_p->jobqueue->queueSem);

    int sval;
    sem_getvalue (tp_p->jobqueue->queueSem, &sval);
}

bool sf_thpool_add_work (sf_thpool tp_p, void *(*function_p)(void*), void* arg_p)
/*< Add work to the thread pool >*/
{
    thpool_job_t* newJob;

    newJob = (thpool_job_t*)malloc (sizeof(thpool_job_t));
    if (newJob == NULL) {
        fprintf (stderr, "thpool_add_work(): Could not allocate memory for new job\n");
        exit (-1);
    }

    /* add function and argument */
    newJob->function = function_p;
    newJob->arg = arg_p;

    /* add job to queue */
    pthread_mutex_lock (&tp_p->mutex);
    thpool_jobqueue_add (tp_p, newJob);
    pthread_mutex_unlock (&tp_p->mutex);

    return true;
}

#endif /* PTHREADS */

