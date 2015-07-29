/* Differential source Fast marching main interface. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
/*^*/

#include "fastm.h"
#include "neighbors.h"

/*#define AT(x,y,z) ((z)+n3*(y)+n23*(x))*/
#define AT(z,y,x) ((z)+n1*(y)+n12*(x))

#define INSERT nm--; mask[i] = FMM_FRONT; pt=ttime+i; sf_pqueue_insert (pt)

#define FMM_OUT '\0'
#define FMM_IN 'i'
#define FMM_FRONT 'f'

static double tp;
static float rd1,rd2,rd3,d1,d2,d3;
static int nm, n12, show=0, n, n1, n2, n3;

void updateds (int p1, int p2, int p3, float* tj, float* dtj, unsigned char* mj, float t, float tx, float ty, float tz, 
	       float s, float y);
void update2ds (int p1, int p2, int p3, float* tj, float* dtj, unsigned char* mj, float t, float tx, float ty, float tz,
	       float s, float y);


void fastm_init (int n3,int n2,int n1) 
/*< Initialize data dimensions >*/
{
    int maxband;
    
    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;

    sf_pqueue_init (10*maxband);
}

void fastm (float* time                /* time */, 
		float* v                   /* slowness */, 
		int* in                    /* in/front/out flag */, 
		bool* plane                /* if plane source */,
		int   n3,  int n2,  int n1 /* dimensions */,
		float o3,float o2,float o1 /* origin */,
		float d3,float d2,float d1 /* sampling */,
		float s3,float s2,float s1 /* source */,
		int   b3,  int b2,  int b1 /* box around the source */,
		int order                  /* accuracy order (1,2,3) */)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], d[3], *p;
    int n[3], b[3], npoints, i;
    
    n[0] = n1; xs[0] = s1-o1; b[0] = b1; d[0] = d1;
    n[1] = n2; xs[1] = s2-o2; b[1] = b2; d[1] = d2;
    n[2] = n3; xs[2] = s3-o3; b[2] = b3; d[2] = d3;

    sf_pqueue_start();
    neighbors_init (in, d, n, order, time);

    for (npoints =  nearsource (xs, b, d, v, plane);
	 npoints > 0;
	 npoints -= neighbours(i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

	/* sf_warning("npoints=%d",npoints); */

	p = sf_pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p - time;

	in[i] = SF_IN;
    }
}

void fastm_close (void)
/*< Free allocated storage >*/
{
    sf_pqueue_close();
}

/* Numerically solving the differential source linear equation */
void fastds (float* time                /* time */, 
		float* v                   /* slowness */, 
		int* in                    /* in/front/out flag */, 
		bool* plane                /* if plane source */,
		int   nn3,  int nn2,  int nn1 /* dimensions */,
		float o3,float o2,float o1 /* origin */,
		float dd3,float dd2,float dd1 /* sampling */,
		float s3,float s2,float s1 /* source */,
		int   b3,  int b2,  int b1 /* box around the source */,
	     int order                  /* accuracy order (1,2,3) */,
	     int sorder                  /* accuracy order (1,2,3) of the source perturbation*/,
	     float dy            /*Source shift*/)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], ***dw, ***tx, ***ty, ***tz, ***dD, *pt, *ptt, *ttime, *dtime, *ddtime;
    int i, i1, i2, i3, j1, j2, nh;
    unsigned char *mask, *pm;

    
    xs[0] = (s1-o1)/dd1; d1 = dd1;  n1 = nn1;
    xs[1] = (s2-o2)/dd2; d2 = dd2;  n2 = nn2;
    xs[2] = (s3-o3)/dd3; d3 = dd3;  n3 = nn3;
    n12 = n2*n1;
    nm=n3*n12;
    n = nm;
    rd1 = 1./d1;
    rd2 = 1./d2;
    rd3 = 1./d3;

    dw  = sf_floatalloc3(n1,n2,n3);
    tx  = sf_floatalloc3(n1,n2,n3);
    ty  = sf_floatalloc3(n1,n2,n3);
    tz  = sf_floatalloc3(n1,n2,n3);

    for(i3 = 0; i3 < n3; i3++){
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,0,i3); j2 = AT(i1,1,i3);
	dw[i3][0][i1] = 0.5*rd2*(v[j2]-v[i]);
      }
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,n2-1,i3); j1 = AT(i1,n2-2,i3);
	dw[i3][n2-1][i1] = 0.5*rd2*(v[i]-v[j1]);
      }
      for(i2 = 1; i2 < n2-1; i2++) {
	for(i1 = 0; i1 < n1; i1++) {
	  i = AT(i1,i2,i3); j1 = AT(i1,i2-1,i3); j2 = AT(i1,i2+1,i3);
	  dw[i3][i2][i1] = 0.25*rd2*(v[j2]-v[j1]);
	}
      }
    }

    for(i3 = 0; i3 < n3; i3++){
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,0,i3); j2 = AT(i1,1,i3);
	ty[i3][0][i1] = rd2*(time[j2]-time[i]);
      }
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,n2-1,i3); j1 = AT(i1,n2-2,i3);
	ty[i3][n2-1][i1] = rd2*(time[i]-time[j1]);
      }
      for(i2 = 1; i2 < n2-1; i2++) {
	for(i1 = 0; i1 < n1; i1++) {
	  i = AT(i1,i2,i3); j1 = AT(i1,i2-1,i3); j2 = AT(i1,i2+1,i3);
	  ty[i3][i2][i1] = 0.5*rd2*(time[j2]-time[j1]);
	}
      }
    }
    
    for(i3 = 0; i3 < n3; i3++){
      for(i2 = 0; i2 < n2; i2++) {
	i = AT(0,i2,i3); j2 = AT(1,i2,i3);
	tz[i3][i2][0] = rd1*(time[j2]-time[i]);
	i = AT(n1-2,i2,i3); j2 = AT(n1-1,i2,i3);
	tz[i3][i2][n1-1] = rd1*(time[j2]-time[i]);
	for(i1 = 1; i1 < n1-1; i1++) {
	  i = AT(i1,i2,i3); j1 = AT(i1-1,i2,i3); j2 = AT(i1+1,i2,i3);
	  tz[i3][i2][i1] = 0.5*rd1*(time[j2]-time[j1]);
	}
      }
    }

    if(n3>1){
      for(i2 = 0; i2 < n2; i2++) {
	for(i1 = 0; i1 < n1; i1++) {
	  i = AT(i1,i2,0); j2 = AT(i1,i2,1);
	  tx[0][i2][i1] = rd3*(time[j2]-time[i]);
	  i = AT(i1,i2,n3-2); j2 = AT(i1,i2,n3-1);
	  tx[n3-1][i2][i1] = rd3*(time[j2]-time[i]);
	}
      }
      for(i3 = 1; i3 < n3-1; i3++){
	for(i2 = 0; i2 < n2; i2++) {
	  for(i1 = 0; i1 < n1; i1++) {
	    i = AT(i1,i2,i3); j1 = AT(i1,i2,i3-1); j2 = AT(i1,i2,i3+1);
	    tx[i3][i2][i1] = 0.5*rd3*(time[j2]-time[j1]);
	  }
	}
      }
    } else {
      for(i3 = 0; i3 < n3; i3++)
	for(i2 = 0; i2 < n2; i2++) 
	  for(i1 = 0; i1 < n1; i1++)
	    tx[i3][i2][i1] = 0.0;
    }


    ttime  = (float *) malloc (n*sizeof(float));
    dtime  = (float *) malloc (n*sizeof(float));
    mask  = (unsigned char *) malloc (n*sizeof(unsigned char));
    for (i=0; i<nm; i++) {
	ttime[i] = SF_HUGE;
	dtime[i] = 0.0;
	mask[i] = FMM_OUT;
    }


    /* nh is an estimate of the maximum front size */
    nh = 0;
    if (n1 > 1) nh += 2*n2*n3;
    if (n2 > 1) nh += 2*n1*n3;
    if (n3 > 1) nh += 2*n1*n2;

    sf_pqueue_init (nh);
    sf_pqueue_start();

    i1 = SF_MAX(SF_MIN(xs[0],n1-2),0);
    i2 = SF_MAX(SF_MIN(xs[1],n2-2),0);
    i3 = SF_MAX(SF_MIN(xs[2],n3-2),0);
    i = AT(i1,i2,i3); INSERT; ptt=dtime+i;
    dtime[i] = SF_SIG(dw[i3][i2][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2][i1]));
    ttime[i] = time[i];

    if(n1>1){
      i = AT(i1+1,i2,i3); INSERT; ptt=dtime+i;
      dtime[i] = SF_SIG(dw[i3][i2][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2][i1+1]));
      ttime[i]= time[i];
    }

    if(n2>1){
      i = AT(i1,i2+1,i3); INSERT; ptt=dtime+i;
      dtime[i] = SF_SIG(dw[i3][i2+1][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2-1)*(xs[1]-i2-1)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2+1][i1]));
      ttime[i]= time[i];
    }

    if(n3>1){
      i = AT(i1,i2,i3+1); INSERT; ptt=dtime+i;
      dtime[i] = SF_SIG(dw[i3+1][i2][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3-1)*(xs[2]-i3-1)*d3*d3)*dw[i3+1][i2][i1]));
      ttime[i]= time[i];
    }

    if(n1>1 && n2>1){
      i = AT(i1+1,i2+1,i3); INSERT; ptt=dtime+i;
      dtime[i] = SF_SIG(dw[i3][i2+1][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2-1)*(xs[1]-i2-1)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2+1][i1+1]));
      ttime[i]= time[i];
    }

    if(n1>1 && n3>1){
      i = AT(i1+1,i2,i3+1); INSERT; ptt=dtime+i;
      dtime[i] = SF_SIG(dw[i3+1][i2][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3-1)*(xs[2]-i3-1)*d3*d3)*dw[i3+1][i2][i1+1]));
      ttime[i]= time[i];
    }

    if(n2>1 && n3>1){
      i = AT(i1,i2+1,i3+1); INSERT; ptt=dtime+i;
      dtime[i] = SF_SIG(dw[i3+1][i2+1][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2+1)*(xs[1]-i2+1)*d2*d2+(xs[2]-i3+1)*(xs[2]-i3+1)*d3*d3)*dw[i3+1][i2+1][i1]));
      ttime[i]= time[i];
    }

    if(n1>1 && n2>1 && n3>1){
      i = AT(i1+1,i2+1,i3+1); INSERT; ptt=dtime+i;
      dtime[i] =  SF_SIG(dw[i3+1][i2+1][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2-1)*(xs[1]-i2-1)*d2*d2+(xs[2]-i3+1)*(xs[2]-i3+1)*d3*d3)*dw[i3+1][i2+1][i1+1]));
      ttime[i]= time[i];
    }


    if(order==1){
      /* start marching */
      while (nm > 0) {  /* "far away" points */
	pt = sf_pqueue_extract ();
	
	i = pt-ttime;
	/*warn("i=%d ttime=%f ttime2=%f rd1=%f rd2=%f rd3=%f nm=%d",i,*pt,*(pt-1),rd1,rd2,rd3,nm);*/
	i1 = i%n1;
	i2 = (i/n1)%n2;
	i3 = i/n12;
	sf_warning("i1=%d i2=%d i3=%d n1=%d n2=%d n3=%d i=%d time=%f,nm=%d",i1,i2,i3,n1,n2,n3,i,time[i],nm);

	*(pm = mask+i) = FMM_IN;
	ptt = dtime+i;
	if (i3 < n3-1 && *(pm+n12) != FMM_IN) updateds (i1,i2,i3+1, pt+n12, ptt+n12, pm+n12,time[i+n12], 
							 tx[i3+1][i2][i1], ty[i3+1][i2][i1], tz[i3+1][i2][i1],
							dw[i3+1][i2][i1],dy); 
	if (i3 > 0    && *(pm-n12) != FMM_IN) updateds (i1,i2,i3-1, pt-n12, ptt-n12, pm-n12, time[i-n12],
							 tx[i3-1][i2][i1], ty[i3-1][i2][i3], tz[i3-1][i2][i1],
							 dw[i3-1][i2][i1],dy);
	if (i2 < n2-1 && *(pm+ n1) != FMM_IN) updateds (i1,i2+1,i3, pt+ n1, ptt+ n1, pm+ n1, time[i+n1],
							 tx[i3][i2+1][i1], ty[i3][i2+1][i1], tz[i3][i2+1][i1],
							 dw[i3][i2+1][i1],dy);
	if (i2 > 0    && *(pm- n1) != FMM_IN) updateds (i1,i2-1,i3, pt- n1, ptt- n1, pm- n1, time[i-n1],
							 tx[i3][i2-1][i1], ty[i3][i2-1][i1], tz[i3][i2-1][i1],
							 dw[i3][i2-1][i1],dy);
	if (i1 < n1-1 && *(pm+  1) != FMM_IN) updateds (i1+1,i2,i3, pt+  1, ptt+  1, pm+  1, time[i+1],
							 tx[i3][i2][i1+1], ty[i3][i2][i1+1], tz[i3][i2][i1+1],
							 dw[i3][i2][i1+1],dy); 
	if (i1 > 0    && *(pm-  1) != FMM_IN) updateds (i1-1,i2,i3, pt-  1, ptt-  1, pm-  1, time[i-1],
							 tx[i3][i2][i1-1], ty[i3][i2][i1-1], tz[i3][i2][i1-1],
							 dw[i3][i2][i1-1],dy);
      }
    } else {

      while (nm > 0) {
	pt = sf_pqueue_extract ();
	i = pt-ttime;
	if(nm<30){
	  show=1;
	  sf_warning("i=%d ttime=%f ttime2=%f rd1=%f rd2=%f rd3=%f nm=%d",i,*pt,*(pt-1),rd1,rd2,rd3,nm);
	}
	i1 = i%n1;
	i2 = (i/n1)%n2;
	i3 = i/n12;
      
	*(pm = mask+i) = FMM_IN;
	ptt = dtime+i;
	if (i3 < n3-1 && *(pm+n12) != FMM_IN) update2ds (i1,i2,i3+1, pt+n12, ptt+n12, pm+n12, time[i+n12], 
							 tx[i3+1][i2][i1], ty[i3+1][i2][i1], tz[i3+1][i2][i1],
							 dw[i3+1][i2][i1],dy); 
	if (i3 > 0    && *(pm-n12) != FMM_IN) update2ds (i1,i2,i3-1, pt-n12, ptt-n12, pm-n12, time[i-n12], 
							 tx[i3-1][i2][i1], ty[i3-1][i2][i1], tz[i3-1][i2][i1],
							 dw[i3-1][i2][i1],dy);
	if (i2 < n2-1 && *(pm+ n1) != FMM_IN) update2ds (i1,i2+1,i3, pt+ n1, ptt+ n1, pm+ n1, time[i+n1],
							 tx[i3][i2+1][i1], ty[i3][i2+1][i1], tz[i3][i2+1][i1],
							 dw[i3][i2+1][i1],dy); 
	if (i2 > 0    && *(pm- n1) != FMM_IN) update2ds (i1,i2-1,i3, pt- n1, ptt- n1, pm- n1, time[i-n1],
							 tx[i3][i2-1][i1], ty[i3][i2-1][i1], tz[i3][i2-1][i1],
							 dw[i3][i2-1][i1],dy);
	if (i1 < n1-1 && *(pm+  1) != FMM_IN) update2ds (i1+1,i2,i3, pt+  1, ptt+  1, pm+  1, time[i+1],
							 tx[i3][i2][i1+1], ty[i3][i2][i1+1], tz[i3][i2][i1+1],
							 dw[i3][i2][i1+1],dy); 
	if (i1 > 0    && *(pm-  1) != FMM_IN) update2ds (i1-1,i2,i3, pt-  1, ptt-  1, pm-  1, time[i-1],
							 tx[i3][i2][i1-1], ty[i3][i2][i1-1], tz[i3][i2][i1-1],
							 dw[i3][i2][i1-1],dy);
      }
    }

    if(sorder==1){
      for (i=0; i<n; i++) {
	/*time[i] = ttime[i];*/
	time[i] -= dtime[i]*dy;
	sf_warning("time=%f time2=%f dtime=%f dy=%f",time[i],time[i],dtime[i],dy);
      }
      return;
    }

     /* end marching */
  sf_pqueue_close ();


  nm = n;
    dD  = sf_floatalloc3(n1,n2,n3);

    for(i3 = 0; i3 < n3; i3++){
      for(i2 = 1; i2 < n2-1; i2++) {
	for(i1 = 0; i1 < n1; i1++) {
	  i = AT(i1,i2,i3); j1 = AT(i1,i2-1,i3); j2 = AT(i1,i2+1,i3);
	  dw[i3][i2][i1] = 0.5*rd2*rd2*(v[j2]-2*v[i]+v[j1]);
	}
      }
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,0,i3); j2 = AT(i1,1,i3);
	dw[i3][0][i1] = dw[i3][1][i1];
      }
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,n2-1,i3); j1 = AT(i1,n2-2,i3);
	dw[i3][n2-1][i1] = dw[i3][n2-2][i1];
      }
    }

    for(i3 = 0; i3 < n3; i3++){
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,0,i3); j2 = AT(i1,1,i3);
	dD[i3][0][i1] = rd2*(dtime[j2]-dtime[i]);
      }
      for(i1 = 0; i1 < n1; i1++) {
	i = AT(i1,n2-1,i3); j1 = AT(i1,n2-2,i3);
	dD[i3][n2-1][i1] = rd2*(dtime[i]-dtime[j1]);
      }
      for(i2 = 1; i2 < n2-1; i2++) {
	for(i1 = 0; i1 < n1; i1++) {
	  i = AT(i1,i2,i3); j1 = AT(i1,i2-1,i3); j2 = AT(i1,i2+1,i3);
	  dD[i3][i2][i1] = 0.5*rd2*(dtime[j2]-dtime[j1]);
	}
      }
    }

    for(i3 = 0; i3 < n3; i3++)
	for(i2 = 0; i2 < n2; i2++) 
	  for(i1 = 0; i1 < n1; i1++)
	    dw[i3][i2][i1] -= dD[i3][i2][i1]*dD[i3][i2][i1];
    
    for(i3 = 0; i3 < n3; i3++){
      for(i2 = 0; i2 < n2; i2++) {
	i = AT(0,i2,i3); j2 = AT(1,i2,i3);
	dD[i3][i2][0] = rd1*(dtime[j2]-dtime[i]);
	i = AT(n1-2,i2,i3); j2 = AT(n1-1,i2,i3);
	dD[i3][i2][n1-1] = rd1*(dtime[j2]-dtime[i]);
	for(i1 = 1; i1 < n1-1; i1++) {
	  i = AT(i1,i2,i3); j1 = AT(i1-1,i2,i3); j2 = AT(i1+1,i2,i3);
	  dD[i3][i2][i1] = 0.5*rd1*(dtime[j2]-dtime[j1]);
	}
      }
    }

    for(i3 = 0; i3 < n3; i3++)
	for(i2 = 0; i2 < n2; i2++) 
	  for(i1 = 0; i1 < n1; i1++)
	    dw[i3][i2][i1] -= dD[i3][i2][i1]*dD[i3][i2][i1];

    if(n3>1){
      for(i2 = 0; i2 < n2; i2++) {
	for(i1 = 0; i1 < n1; i1++) {
	  i = AT(i1,i2,0); j2 = AT(i1,i2,1);
	  dD[0][i2][i1] = rd3*(dtime[j2]-dtime[i]);
	  i = AT(i1,i2,n3-2); j2 = AT(i1,i2,n3-1);
	  dD[n3-1][i2][i1] = rd3*(dtime[j2]-dtime[i]);
	}
      }
      for(i3 = 1; i3 < n3-1; i3++){
	for(i2 = 0; i2 < n2; i2++) {
	  for(i1 = 0; i1 < n1; i1++) {
	    i = AT(i1,i2,i3); j1 = AT(i1,i2,i3-1); j2 = AT(i1,i2,i3+1);
	    dD[i3][i2][i1] = 0.5*rd3*(dtime[j2]-dtime[j1]);
	  }
	}
      }

      for(i3 = 0; i3 < n3; i3++)
	for(i2 = 0; i2 < n2; i2++) 
	  for(i1 = 0; i1 < n1; i1++)
	    dw[i3][i2][i1] -= dD[i3][i2][i1]*dD[i3][i2][i1];
    }

    ddtime  = (float *) malloc (n*sizeof(float));
    for (i=0; i<nm; i++) {
	ttime[i] = SF_HUGE;
	ddtime[i] = 0.0;
	mask[i] = FMM_OUT;
    }
    
    sf_pqueue_init (nh);
    sf_pqueue_start();

    i1 = SF_MAX(SF_MIN(xs[0],n1-2),0);
    i2 = SF_MAX(SF_MIN(xs[1],n2-2),0);
    i3 = SF_MAX(SF_MIN(xs[2],n3-2),0);
    i = AT(i1,i2,i3); INSERT; ptt=dtime+i;
    ddtime[i] = SF_SIG(dw[i3][i2][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2][i1]));
    ttime[i] = time[i];

    if(n1>1){
      i = AT(i1+1,i2,i3); INSERT; ptt=dtime+i;
      ddtime[i] = SF_SIG(dw[i3][i2][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2][i1+1]));
      ttime[i] = time[i];
    }

    if(n2>1){
      i = AT(i1,i2+1,i3); INSERT; ptt=dtime+i;
      ddtime[i] = SF_SIG(dw[i3][i2+1][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2-1)*(xs[1]-i2-1)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2+1][i1]));
      ttime[i] = time[i];
    }

    if(n3>1){
      i = AT(i1,i2,i3+1); INSERT; ptt=dtime+i;
      ddtime[i] = SF_SIG(dw[i3+1][i2][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3-1)*(xs[2]-i3-1)*d3*d3)*dw[i3+1][i2][i1]));
      ttime[i] = time[i]; 
    }

    if(n1>1 && n2>1){
      i = AT(i1+1,i2+1,i3); INSERT; ptt=dtime+i;
      ddtime[i] = SF_SIG(dw[i3][i2+1][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2-1)*(xs[1]-i2-1)*d2*d2+(xs[2]-i3)*(xs[2]-i3)*d3*d3)*dw[i3][i2+1][i1+1]));
      ttime[i] = time[i];
    }

    if(n1>1 && n3>1){
      i = AT(i1+1,i2,i3+1); INSERT; ptt=dtime+i;
      ddtime[i] = SF_SIG(dw[i3+1][i2][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2)*(xs[1]-i2)*d2*d2+(xs[2]-i3-1)*(xs[2]-i3-1)*d3*d3)*dw[i3+1][i2][i1+1]));
      ttime[i] = time[i];
    }

    if(n2>1 && n3>1){
      i = AT(i1,i2+1,i3+1); INSERT; ptt=dtime+i;
      ddtime[i] = SF_SIG(dw[i3+1][i2+1][i1])*sqrt(SF_ABS(((xs[0]-i1)*(xs[0]-i1)*d1*d1+(xs[1]-i2+1)*(xs[1]-i2+1)*d2*d2+(xs[2]-i3+1)*(xs[2]-i3+1)*d3*d3)*dw[i3+1][i2+1][i1]));
      ttime[i] = time[i];
    }

    if(n1>1 && n2>1 && n3>1){
      i = AT(i1+1,i2+1,i3+1); INSERT; ptt=dtime+i;
      ddtime[i] =  SF_SIG(dw[i3+1][i2+1][i1+1])*sqrt(SF_ABS(((xs[0]-i1-1)*(xs[0]-i1-1)*d1*d1+(xs[1]-i2-1)*(xs[1]-i2-1)*d2*d2+(xs[2]-i3+1)*(xs[2]-i3+1)*d3*d3)*dw[i3+1][i2+1][i1+1]));
      ttime[i] = time[i];
    }

    if(order==1){
      /* start marching */
      while (nm > 0) {  /* "far away" points */
	pt = sf_pqueue_extract ();
	
	i = pt-ttime;
	/*sf_warning("i=%d ttime=%f ttime2=%f rd1=%f rd2=%f rd3=%f nm=%d",i,*pt,*(pt-1),rd1,rd2,rd3,nm);*/
	i1 = i%n1;
	i2 = (i/n1)%n2;
	i3 = i/n12;

	*(pm = mask+i) = FMM_IN;
	ptt = ddtime+i;
	if (i3 < n3-1 && *(pm+n12) != FMM_IN) updateds (i1,i2,i3+1, pt+n12, ptt+n12, pm+n12,time[i+n12], 
							 tx[i3+1][i2][i1], ty[i3+1][i2][i1], tz[i3+1][i2][i1],
							dw[i3+1][i2][i1],dy); 
	if (i3 > 0    && *(pm-n12) != FMM_IN) updateds (i1,i2,i3-1, pt-n12, ptt-n12, pm-n12, time[i-n12],
							 tx[i3-1][i2][i1], ty[i3-1][i2][i3], tz[i3-1][i2][i1],
							 dw[i3-1][i2][i1],dy);
	if (i2 < n2-1 && *(pm+ n1) != FMM_IN) updateds (i1,i2+1,i3, pt+ n1, ptt+ n1, pm+ n1, time[i+n1],
							 tx[i3][i2+1][i1], ty[i3][i2+1][i1], tz[i3][i2+1][i1],
							 dw[i3][i2+1][i1],dy);
	if (i2 > 0    && *(pm- n1) != FMM_IN) updateds (i1,i2-1,i3, pt- n1, ptt- n1, pm- n1, time[i-n1],
							 tx[i3][i2-1][i1], ty[i3][i2-1][i1], tz[i3][i2-1][i1],
							 dw[i3][i2-1][i1],dy);
	if (i1 < n1-1 && *(pm+  1) != FMM_IN) updateds (i1+1,i2,i3, pt+  1, ptt+  1, pm+  1, time[i+1],
							 tx[i3][i2][i1+1], ty[i3][i2][i1+1], tz[i3][i2][i1+1],
							 dw[i3][i2][i1+1],dy); 
	if (i1 > 0    && *(pm-  1) != FMM_IN) updateds (i1-1,i2,i3, pt-  1, ptt-  1, pm-  1, time[i-1],
							 tx[i3][i2][i1-1], ty[i3][i2][i1-1], tz[i3][i2][i1-1],
							 dw[i3][i2][i1-1],dy);
      }
    } else {

      while (nm > 0) {
	pt = sf_pqueue_extract ();
	i = pt-ttime;
	if(nm<30){
	  show=1;
	  sf_warning("i=%d ttime=%f ttime2=%f rd1=%f rd2=%f rd3=%f nm=%d",i,*pt,*(pt-1),rd1,rd2,rd3,nm);
	}
	i1 = i%n1;
	i2 = (i/n1)%n2;
	i3 = i/n12;
      
	*(pm = mask+i) = FMM_IN;
	ptt = ddtime+i;
	if (i3 < n3-1 && *(pm+n12) != FMM_IN) update2ds (i1,i2,i3+1, pt+n12, ptt+n12, pm+n12, time[i+n12], 
							 tx[i3+1][i2][i1], ty[i3+1][i2][i1], tz[i3+1][i2][i1],
							 dw[i3+1][i2][i1],dy); 
	if (i3 > 0    && *(pm-n12) != FMM_IN) update2ds (i1,i2,i3-1, pt-n12, ptt-n12, pm-n12, time[i-n12], 
							 tx[i3-1][i2][i1], ty[i3-1][i2][i1], tz[i3-1][i2][i1],
							 dw[i3-1][i2][i1],dy);
	if (i2 < n2-1 && *(pm+ n1) != FMM_IN) update2ds (i1,i2+1,i3, pt+ n1, ptt+ n1, pm+ n1, time[i+n1],
							 tx[i3][i2+1][i1], ty[i3][i2+1][i1], tz[i3][i2+1][i1],
							 dw[i3][i2+1][i1],dy); 
	if (i2 > 0    && *(pm- n1) != FMM_IN) update2ds (i1,i2-1,i3, pt- n1, ptt- n1, pm- n1, time[i-n1],
							 tx[i3][i2-1][i1], ty[i3][i2-1][i1], tz[i3][i2-1][i1],
							 dw[i3][i2-1][i1],dy);
	if (i1 < n1-1 && *(pm+  1) != FMM_IN) update2ds (i1+1,i2,i3, pt+  1, ptt+  1, pm+  1, time[i+1],
							 tx[i3][i2][i1+1], ty[i3][i2][i1+1], tz[i3][i2][i1+1],
							 dw[i3][i2][i1+1],dy); 
	if (i1 > 0    && *(pm-  1) != FMM_IN) update2ds (i1-1,i2,i3, pt-  1, ptt-  1, pm-  1, time[i-1],
							 tx[i3][i2][i1-1], ty[i3][i2][i1-1], tz[i3][i2][i1-1],
							 dw[i3][i2][i1-1],dy);
      }
    }

    if(sorder==2){
      for (i=0; i<n; i++) {
	/*time[i] = ttime[i];*/
	time[i] -= dtime[i]*dy+0.5*ddtime[i]*dy*dy;
	sf_warning("time2=%f ddtime=%f dy=%f",time[i],ddtime[i],dy*dy);
      }
    } else {
      for (i=0; i<n; i++) {
	/*time[i] = ttime[i];*/
	time[i] -= dtime[i]*dtime[i]*dy/(dtime[i]-0.5*ddtime[i]*dy);
	sf_warning("time2=%f dtime==%f ddtime==%f ddtimeS=%f dy=%f",time[i],dtime[i]*dy,dtime[i]*dy-0.5*ddtime[i]*dy*dy,dtime[i]*dtime[i]*dy/(0.5*ddtime[i]*dy+dtime[i]),dy*dy);
      }
    }
   
}

void updateds (int p1, int p2, int p3, float* tj, float* dtj, unsigned char* mj, float t, float tx, float ty, float tz, 
	       float s, float dy)
{
  float b, c, t1=0., t2=0., u, den,tp1,dt1=0.;
  unsigned int k, i;

  b = c = 0; i = k = 0;
  if ((p3 > 0   ) && *(mj-n12) && ((t1 = *(tj-n12)) < *tj)){i |= 0x01; dt1 = *(dtj-n12);}
  if ((p3 < n3-1) && *(mj+n12) && ((t2 = *(tj+n12)) < *tj) && 
      ((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; dt1 = *(dtj+n12);} 
  if (i & 0x01) {
    u = rd3*SF_ABS(tx); b += u*dt1; c += u;
    i ^= 0x01; k |= 0x01;
  }
  if ((p2 > 0   ) && *(mj-n1) && ((t1 = *(tj-n1)) < *tj)) {i |= 0x01; dt1 = *(dtj-n1);}
  if ((p2 < n2-1) && *(mj+n1) && ((t2 = *(tj+n1)) < *tj) && 
      ((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; dt1 = *(dtj+n1);} 
  if (i & 0x01) {
    u = rd2*SF_ABS(ty); b += u*dt1; c += u;
    i ^= 0x01; k |= 0x02;
  }

  if ((p1 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) {i |= 0x01; dt1 = *(dtj-1);}
  if ((p1 < n1-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
      ((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; dt1 = *(dtj+1);}  
  if (i & 0x01) {
    u = rd1*SF_ABS(tz); b += u*dt1; c += u;
    i ^= 0x01; k |= 0x04;
  }

  if (!k) return;

  den = (SF_ABS(c) < 0.00000001 ? SF_SIG(c)*100000000. : 1./c);
  /*tp1 = (b-(sqrt(1+2.*eta)-1)*sqrt(1+2.*eta0)*tr2*(1-rsv*tz*tz))*den;*/
  tp1 = (b-s)*den;
  tp = t -tp1*dy;

  if(SF_ABS(c) < 0.00000001) 
    sf_error("stop p1=%d p2=%d p3=%d t=%f k=%d tp1=%f b=%f c=%f tx=%f ty=%f tz=%f",p1,p2,p3,tp,k,tp1,b,c,tx,ty,tz);

  if (t < *tj) {
    *tj = t; *dtj = tp1;
    if (*mj == FMM_OUT) {
      nm--; 
      *mj = FMM_FRONT; 
      sf_pqueue_insert (tj);
    } 
  }
}

void update2ds (int p1, int p2, int p3, float* tj, float* dtj, unsigned char* mj, float t, float tx, float ty, float tz,
		float s, float dy)
{
  float t1=0., t2=0., u, dt1=0.;
  double den;
  unsigned int k, i;
  double bbb,ccc,ddd,tp1;
  unsigned int jjj;

  i = k = 0;
  ddd= 1.0; jjj = 0;
  bbb = ccc = 0.;

  if ((p3 > 0   ) && *(mj-n12) && ((t1 = *(tj-n12)) < *tj)) {i |= 0x01; dt1 = *(dtj-n12);}
  if ((p3 < n3-1) && *(mj+n12) && ((t2 = *(tj+n12)) < *tj) && 
      ((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; dt1 = *(dtj+n12);} 
  if (i & 0x01) {
    ddd=rd3;
    if ((jjj & 0x01) && (p3 < n3-2) && *(mj+2*n12)) { 
      ddd *= 1.5; dt1 *= 4.; dt1 -= *(dtj+2*n12); dt1 /= 3.;
    } else if ((p3 > 1) && *(mj-2*n12)) {   
      ddd *= 1.5; dt1 *= 4.; dt1 -= *(dtj-2*n12); dt1 /= 3.;
    }
    u = ddd*SF_ABS(tx); bbb += u*dt1; ccc += u;
    i ^= 0x01; k |= 0x01;
  }
  jjj =0;

  if ((p2 > 0   ) && *(mj-n1) && ((t1 = *(tj-n1)) < *tj)) {i |= 0x01; dt1 = *(dtj-n1);}
  if ((p2 < n2-1) && *(mj+n1) && ((t2 = *(tj+n1)) < *tj) && 
      ((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; dt1 = *(dtj+n1);} 
  if (i & 0x01) {
    ddd=rd2;
    if ((jjj & 0x01) && (p2 < n2-2) && *(mj+2*n1)) { 
      ddd *= 1.5; dt1 *= 4.; dt1 -= *(dtj+2*n1); dt1 /= 3.; 
    } else if ((p2 > 1) && *(mj-2*n1)) {   
      ddd *= 1.5; dt1 *= 4.; dt1 -= *(dtj-2*n1); dt1 /= 3.;
    }
    u = ddd*SF_ABS(ty); bbb += u*dt1; ccc += u;
    i ^= 0x01; k |= 0x02;
  }
  jjj =0;

  if ((p1 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) {i |= 0x01; dt1 = *(dtj-1);}
  if ((p1 < n1-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
      ((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; dt1 = *(dtj+1);}  
  if (i & 0x01) {
    ddd=rd1;
    if ((jjj & 0x01) && (p1 < n1-2) && *(mj+2)) { 
      ddd *= 1.5; dt1 *= 4.; dt1 -= *(dtj+2); dt1 /= 3.; 
    } else if ((p1 > 1) && *(mj-2)) {   
      ddd *= 1.5; dt1 *= 4.; dt1 -= *(dtj-2); dt1 /= 3.;
    }
    u = ddd*SF_ABS(tz); bbb += u*dt1; ccc += u;
    i ^= 0x01; k |= 0x04;
  }
  jjj =0;

  if (!k) return;

  den = (SF_ABS(ccc) < 0.00000001 ? SF_SIG(ccc)*100000000. : 1./ccc);
  /*den = (SF_ABS(ccc) < 0.00000001 ? 0.0 : 1./ccc);*/
  tp1 = (bbb-s)*den;
  /*tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(tp1),0.5);*/
  tp = t-tp1*dy;
  if(show)
    sf_warning("p1=%d p2=%d p3=%d t=%f k=%d tp1=%f bbb=%f ccc=%f tx=%f ty=%f tz=%f",p1,p2,p3,tp,k,tp1,bbb,ccc,tx,ty,tz);

  if (t < *tj && t>0) {
    *tj = t; *dtj = tp1;
    if (*mj == FMM_OUT) {
      nm--; 
      *mj = FMM_FRONT; 
      sf_pqueue_insert (tj);
    } 
  }
}

/* 	$Id: fastm.c 4136 2009-05-07 17:20:32Z tariq $	 */
