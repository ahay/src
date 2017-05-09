/* Testing Posix threads */
/*
  Copyright (C) 2009 Dough McCowan
   
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
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#define MAX_QUEUE_SIZE 37   /* max number of work orders in queue */
#define MAX_THREAD_FOLD 8   /* max number of migrations per thread call */

typedef struct mig_work {
    int id; int ntodo; float* data;
    float xdist2[MAX_THREAD_FOLD]; float ydist2[MAX_THREAD_FOLD];
    float shot_dist[MAX_THREAD_FOLD]; float rcvr_dist[MAX_THREAD_FOLD]; 
    int offsets[MAX_THREAD_FOLD];
    float* lineap[MAX_THREAD_FOLD]; float* cdpap[MAX_THREAD_FOLD]; 
    float* vrms[MAX_THREAD_FOLD];
    float* image[MAX_THREAD_FOLD]; float* image_fold[MAX_THREAD_FOLD];
    int migrated_gathers[MAX_THREAD_FOLD]; float* gathers[MAX_THREAD_FOLD]; 
    float* gathers_fold[MAX_THREAD_FOLD]; struct mig_work *next;
} mig_work_order_t;

typedef struct tpool {
    int   id;
    int   *working;                   /* me working flags */
    int in_nsamps;                    /* never changes */
    int *hitarray;                    /* never changes */
    int digi;                         /* never changes */
    float stretch_mute;               /* never changes */
    int ntimes;                       /* never changes */
    int amps;                         /* never changes */
    float dxbar;                      /* never changes */
    int   cur_queue_size;
    int   max_queue_size;
    mig_work_order_t *queue_head;
    mig_work_order_t *queue_tail;
    pthread_mutex_t queue_lock;
    pthread_cond_t queue_empty;
    pthread_cond_t queue_not_empty;
    pthread_cond_t queue_not_full;
    int   shutdown;
} *tpool_t;

#include <rsf.h>

int main(int argc, char **argv)
{
    int         h, i, j, k, n;
    int         nkernels, id, kcalls=0, fixed_up=0, nfilled=0, nread=0;
    float       stretch_mute=0.0, dxbar=0.0;
    int         amps=0;
    int         input_number_samples=0;
    int         digi=0, ntimes=0, onetwo;
    float       a,b=1.0;
    float       *temp1=NULL, *temp2=NULL, *temp;

    unsigned int delay=1000;
    tpool_t tpool;
    pthread_t *peers;
    pthread_attr_t attr;
    mig_work_order_t *work_order=NULL;

    void   post_work_order(tpool_t tpool, mig_work_order_t *work_order);
    void   peer_mig(tpool_t tpool);

    sf_init(argc,argv);

    if (!sf_getint("kernels",&nkernels)) sf_error("Need kernels=");
    /* Number of kernel threads to create */

    if (!sf_getint("times",&ntimes)) sf_error("Need times=");
    /* Number of SQRT loops to execute */

    tpool = (tpool_t) sf_alloc(1,sizeof(struct tpool));
    tpool->cur_queue_size = 0;
    tpool->max_queue_size = MAX_QUEUE_SIZE;
    tpool->queue_head = NULL;
    tpool->queue_tail = NULL;
    tpool->shutdown = 0;
    pthread_mutex_init(&(tpool->queue_lock), NULL);
    pthread_cond_init(&(tpool->queue_empty), NULL);
    pthread_cond_init(&(tpool->queue_not_empty), NULL);
    pthread_cond_init(&(tpool->queue_not_full), NULL);

    peers = (pthread_t*) sf_alloc(nkernels,sizeof(pthread_t));
    pthread_attr_init(&attr);        /* set default thread attributes */
#ifdef SGI
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_BOUND_NP);  /* allow SGI single process */
    pthread_setconcurrency(nkernels);     /* allow SGI multi kernels */
#else
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);  /* allow multi process */
#endif
    temp = temp1; onetwo=1;

/* kickoff the kernels into a wait state */

    tpool->working = sf_intalloc(nkernels);
    tpool->in_nsamps = input_number_samples;    /* never changes */
    tpool->hitarray = sf_intalloc(nkernels);
    tpool->digi = digi;                         /* never changes */
    tpool->stretch_mute = stretch_mute;         /* never changes */
    tpool->ntimes = ntimes;                     /* never changes */
    tpool->amps = amps;                         /* never changes */
    tpool->dxbar = dxbar;                       /* never changes */
    for(k=0; k<nkernels; k++) tpool->working[k] = 0;

    for(id=0; id<nkernels; id++) {
	tpool->id = id;
	n=pthread_create(&(peers[id]), &attr, (void*(*)(void*))peer_mig,
			 (void*)tpool);
	if(n) { sf_warning("\n peer %d kickoff error %d",id,n); exit(0); }
	sf_warning("peer %d kicked off",id);
	usleep(delay);    /* give threads time to retrieve parameters */
    }

    nread = 0;
    for(h=0; h<4; h++) {              /* loop over input traces */
	nread++;

	for(k=0; k<nkernels; k++) tpool->hitarray[k] = 0;
	fixed_up = 0;
	for(i=0; i<3; i++) {            /* loop over lines */
	    for(j=0; j<73; j++) {         /* loop over xlines */

/* any aperture branch-arounds would be here */

		if (fixed_up == 0) {

		    sf_warning("\npreparing trace %d", nread); 
		    for(k=0; k<10000; k++) a=sqrt(b);
     
		    fixed_up = 1;
		    nfilled = 0;
		    kcalls = 0;
		    tpool->cur_queue_size = 0;      /* just to */
		    tpool->queue_head = NULL;         /* be */
		    tpool->queue_tail = NULL;        /* sure */
		}

/* fill out the new work order for this batch */

		k=kcalls%MAX_THREAD_FOLD;

		work_order = (mig_work_order_t*) sf_alloc(1,sizeof(mig_work_order_t));

		if(!k) {
		    work_order->id = (nread-1)*1000+kcalls; 
		    work_order->ntodo = 0; 
		    work_order->next = NULL; 
		} 

/* fill out this entry in the work order */

		work_order->data = temp;         /* all the same for dis work order */
		work_order->xdist2[k] = 0.0;
		work_order->ydist2[k] = 0.0;
		work_order->shot_dist[k] = 0.0;
		work_order->rcvr_dist[k] = 0.0;
		work_order->offsets[k] = 0;
		work_order->lineap[k] = NULL;
		work_order->cdpap[k] = NULL;
		work_order->vrms[k] = NULL;
		work_order->image[k] = NULL;
		work_order->image_fold[k] = NULL;
		work_order->migrated_gathers[k] = 0;
		work_order->gathers[k] = NULL;
		work_order->gathers_fold[k] = NULL;
		(work_order->ntodo)++;                /* one more in dis batch */
		nfilled++;

/* send it off for migration */

		if(nfilled==MAX_THREAD_FOLD) {
		    sf_warning(" main posting filled work order %d with %d vectors for trace %d" ,
			       work_order->id,nfilled,nread);  
		    post_work_order(tpool, work_order);
		    nfilled = 0;
		}
 
		kcalls++;

	    } /* end loop j */
	} /* end loop i */

/* send last unfilled work order off for migration */

	if(nfilled) {
	    sf_warning(" main posting partially-filled work order %d with %d vectors for trace %d" ,
		       work_order->id,nfilled,nread);  
	    post_work_order(tpool, work_order);
	    nfilled = 0;
	}
 
/* wait for queue to empty before migrating next trace */

	pthread_mutex_lock(&(tpool->queue_lock));    /* lock queue */
	while(tpool->cur_queue_size>0) {             /* check flag */
	    sf_warning(" Waiting in main on queue_size= %d",tpool->cur_queue_size);
	    pthread_cond_wait(&(tpool->queue_empty), &(tpool->queue_lock));
	} 
	pthread_mutex_unlock(&(tpool->queue_lock));  /* unlock queue */

/* switch temp array */

	onetwo = 3-onetwo;
	if(onetwo == 1) temp = temp1;
	else            temp = temp2;

    } /* end h loop */

/* wait for kernels to finish last batch of work orders */

    for(k=0; k<nkernels; k++) {while(tpool->working[k]>0) usleep(delay);} 

/* kill kernel threads */

    sf_warning("\n Killing off kernels");
    pthread_mutex_lock(&(tpool->queue_lock));    /* lock queue */
    tpool->shutdown = 1;                         /* request shutdown */
    tpool->cur_queue_size=1;                     /* one in queue */
    pthread_mutex_unlock(&(tpool->queue_lock));  /* unlock queue */
    pthread_cond_broadcast(&(tpool->queue_not_empty)); /* broadcast queue not empty */
    for(id=0; id<nkernels; id++) {
	pthread_join(peers[id], NULL);               /* wait for kernel to flush */
    }


} /* end */


void peer_mig(tpool_t tpool)
{
    int  i,j,me,ntimes;
    int  *working;
    float a,b;
    mig_work_order_t *my_work_order;

/* Routine to start off a peer thread */

    me = tpool->id;
    ntimes = tpool->ntimes;
    working = tpool->working;
    sf_warning(" Hello from thread %d I'm alive and well",me);

/* loop forever */

    for(;;) {

/* sleep while queue empty */

	pthread_mutex_lock(&(tpool->queue_lock));    /* lock queue */
	while( (tpool->cur_queue_size==0) && (!tpool->shutdown) ) { 
	    sf_warning(" Thread %d waiting for work",me);
	    pthread_cond_wait(&(tpool->queue_not_empty), 
			      &(tpool->queue_lock));  /* wait for work */
	}

/* quit if asked to shutdown */

	if(tpool->shutdown) {
	    sf_warning(" Thread %d shutting down",me);
	    pthread_mutex_unlock(&(tpool->queue_lock));  /* unlock queue */
	    pthread_exit(NULL);                          /* quit */
	}

/* otherwise get a work order from the head of the list */

	if(tpool->cur_queue_size==0) continue;      /* just to be careful */
	my_work_order = tpool->queue_head;
	tpool->cur_queue_size--;
	if(tpool->cur_queue_size==0) tpool->queue_head = NULL;
	else                         tpool->queue_head = my_work_order->next;
      
/* now migrate da data */

	if(tpool->cur_queue_size==tpool->max_queue_size-1)  /* signal main thread */
	    pthread_cond_signal(&(tpool->queue_not_full));    /* queue not full */
	if(tpool->cur_queue_size==0)                        /* signal main thread */
	    pthread_cond_signal(&(tpool->queue_empty));       /* queue empty */
	working[me] = 1;                                    /* me now working */
	pthread_mutex_unlock(&(tpool->queue_lock));         /* unlock queue */

	sf_warning(" Thread %d migrating work order %d ntodo= %d",
		   me,my_work_order->id,my_work_order->ntodo);             
/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   call migration kernel routine here 
   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

	for(i=0; i<ntimes; i++) {
	    b=(float)i;
	    for(j=0; j<ntimes; j++)  a=sqrt(b);
	}
	working[me] = 0;                                    /* me finish working */

	sf_warning(" Thread %d finished migrating work order %d", me,my_work_order->id);                   
	free(my_work_order); 
    }
}


void post_work_order(tpool_t tpool, mig_work_order_t *work_order)
{

/* wait for queue size to decrease */

    pthread_mutex_lock(&(tpool->queue_lock));    /* lock queue */
    while(tpool->cur_queue_size>=tpool->max_queue_size) {  /* check queue length */
	pthread_cond_wait(&(tpool->queue_not_full), &(tpool->queue_lock));
    } 
          
/* now add work order to queue linked list */

    if(tpool->cur_queue_size==0) {      /* if queue empty */
	tpool->queue_head = work_order;   /* link it in as queue_head */
	tpool->queue_tail = work_order;   /* and also as queue_tail */
	pthread_cond_broadcast(&(tpool->queue_not_empty)); /* broadcast queue not empty */
    } else {                                    /* if queue not empty */
	(tpool->queue_tail)->next = work_order;   /* link it in as queue_tail */
	tpool->queue_tail = work_order;
    }
    tpool->cur_queue_size++;        /* one more queue entry now */
    pthread_mutex_unlock(&(tpool->queue_lock));  /* unlock queue */

    return;
}


