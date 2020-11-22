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

#include "fastvti.h"

#define AT(x,y,z) ((z)+n3*(y)+n23*(x))
/*#define AT(z,y,x) ((z)+n1*(y)+n12*(x))*/

#define INSERT nm--; mask[j] = FMM_FRONT; pt=ttime+j; sf_pqueue_insert (pt)

#define FMM_OUT '\0'
#define FMM_IN 'i'
#define FMM_FRONT 'f'

static double tp;
static float a, dd[8], rd1,rd2,rd3,d1,d2,d3;
static int nm, n23, n, n1, n2, n3, jj2, jj3;

static void update (int p1, int p2, int p3, float* tk, float* tj, unsigned char* mj, float s, float rsv, float eta);
static void update2 (int p1, int p2, int p3, float* tj, unsigned char* mj, float s, float rsv, float eta);
static void update3 (int p1, int p2, int p3, float* tj, unsigned char* mj, float s, float rsv, float eta);

/* Numerically solving the differential source linear equation */
void fastvti (float* time                /* time */, 
	      float* slow                   /* NMO velocity squared */, 
	      float* slowv                   /* vertical velocity squared */, 
	      float* eta                   /*  eta */, 
	      int* in                    /* in/front/out flag */, 
	      bool* plane                /* if plane source */,
	      int   nn3,  int nn2,  int nn1 /* dimensions */,
	      float o3,float o2,float o1 /* origin */,
	      float dd3,float dd2,float dd1 /* sampling */,
	      float ss3,float ss2,float ss1 /* source */,
	      int   b3,  int b2,  int b1 /* box around the source */,
	      int order                  /* accuracy order (1,2,3) */)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], *pt, *ttime;
    float sl = 0., svl = 0., etal = 0., sh;
    float s2, sv2, sint2, cost2, dist2, ss;
    int i, i1, i2, i3, j1, j2, j3, nh, j;
    unsigned char *mask, *pm;
    
    xs[2] = (ss1-o1)/dd1; d3 = dd1;  n3 = nn1;
    xs[1] = (ss2-o2)/dd2; d2 = dd2;  n2 = nn2;
    xs[0] = (ss3-o3)/dd3; d1 = dd3;  n1 = nn3;
    n23 = n2*n3;
    nm=n1*n23;
    n = nm;
    rd1 = 1./d1;
    rd2 = 1./d2;
    rd3 = 1./d3;

    j1 = SF_MAX(SF_MIN(xs[0],n1-2),0);
    j2 = SF_MAX(SF_MIN(xs[1],n2-2),0);
    j3 = SF_MAX(SF_MIN(xs[2],n3-2),0);

    /* should probably check the return value of malloc */
    ttime  = (float *) malloc (n*sizeof(float));
    mask  = (unsigned char *) malloc (n*sizeof(unsigned char));

    slow[0] = 1./slow[1];
    slowv[0] = slowv[1];
    eta[0] = eta[1];
    ttime[0] = SF_HUGE;
    mask[1] = FMM_OUT;

    for (i=1; i<n; i++) {
	ttime[i] = SF_HUGE;
	mask[i] = FMM_OUT;
	slow[i] = 1./slow[i];
	slowv[i] = slowv[i];
    }

    j = AT(j1,j2,j3);
    s2 = slow[j];
    sf_warning("j=%d slow=%f",j,slow[0]);
    sl = sqrt(s2);
    sv2 = 1./slowv[j];
    svl = sqrt(sv2);
    etal = eta[j];
    sh = sl/sqrt(1+2*etal);
    ttime[j] = 0.;
    mask[j] = FMM_IN;
    nm = n-1;

    /* nh is an estimate of the maximum front size */
    nh = 0;
    if (n1 > 1) nh += 2*n2*n3;
    if (n2 > 1) nh += 2*n1*n3;
    if (n3 > 1) nh += 2*n1*n2;

    sf_pqueue_init(nh);
    sf_pqueue_start();

    /* initialize source */
    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
	if (i1 >= 0 && i1 < n1) {
	    j = AT(i1,j2,j3); INSERT;
	    *pt = sh * d1;
	}
    }
    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
	if (i2 >= 0 && i2 < n2) {
	    j = AT(j1,i2,j3); INSERT;
	    *pt = sh * d2;
	}
    }
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    j = AT(j1,j2,i3); INSERT;
	    *pt = svl * d3; 
	}
    }
    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
	if (i2 >= 0 && i2 < n2) {
	    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
		if (i1 >= 0 && i1 < n1) {
		    j = AT(i1,i2,j3); INSERT;
		    *pt = sh * sqrt(d1*d1+d2*d2);
		}
	    }
	}
    }

    dist2 = d1*d1+d3*d3;
    cost2 = d3*d3/dist2;
    sint2 = 1-cost2;
    ss=2*s2*sv2/(cost2*s2+sint2*(sv2+2*etal*sv2) +
		 sqrt(-8*cost2*etal*s2*sint2*sv2+(cost2*s2+(1 +2*etal)*sint2*sv2)
		      *(cost2*s2+(1+2*etal)*sint2*sv2)));
    dist2 = sqrt(ss*dist2);
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
		if (i1 >= 0 && i1 < n1) {
		    j = AT(i1,j2,i3); INSERT;
		    *pt = dist2;
		}
	    }
	}
    }

    dist2 = d2*d2+d3*d3;
    cost2 = d3*d3/dist2;
    sint2 = 1-cost2;
    ss=2*s2*sv2/(cost2*s2+sint2*(sv2+2*etal*sv2) +
		 sqrt(-8*cost2*etal*s2*sint2*sv2+(cost2*s2+(1 +2*etal)*sint2*sv2)
		      *(cost2*s2+(1+2*etal)*sint2*sv2)));
    dist2 = sqrt(ss*dist2);
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
		if (i2 >= 0 && i2 < n2) {
		    j = AT(j1,i2,i3); INSERT;
		    *pt = dist2;
		}
	    }
	}
    }

    dist2 = d1*d1+d2*d2+d3*d3;
    cost2 = d3*d3/dist2;
    sint2 = 1-cost2;
    ss=2*s2*sv2/(cost2*s2+sint2*(sv2+2*etal*sv2) +
		 sqrt(-8*cost2*etal*s2*sint2*sv2+(cost2*s2+(1 +2*etal)*sint2*sv2)
		      *(cost2*s2+(1+2*etal)*sint2*sv2)));
    dist2 = sqrt(ss*dist2);
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
		if (i2 >= 0 && i2 < n2) {
		    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
			if (i1 >= 0 && i1 < n1) {
			    j = AT(i1,i2,i3); INSERT;
			    *pt = dist2;
			}
		    }
		}
	    }
	}
    }
    /* source initialized */
    sf_warning("v=%f ss=%f dist2=%f nm=%d",slow[1],ss,dist2,nm);

    /* precompute some of the finite-difference coefficients */
    dd[0] = 0.;
    dd[1] = d1*d1;
    dd[2] = d2*d2;
  
    d1 = 1./dd[1];
    d2 = 1./dd[2];
    d3 = 1./(d3*d3);

    dd[3] = 1./(d1+d2);

    if(order==1){
	/* start marching */
	while (nm > 0) {  /* "far away" points */
	    pt = sf_pqueue_extract ();
	    i = pt-ttime;
	    i3 = i%n3;
	    i2 = (i/n3)%n2;
	    i1 = i/n23;
	    sf_warning("i1=%d i2=%d i3=%d i=%d time=%f slow=%f slow=%f eta=%f nm=%d",i1,i2,i3,i,ttime[i],slow[i],slowv[i],eta[i],nm);

	    *(pm = mask+i) = FMM_IN;
	    if (i1 < n1-1 && *(pm+n23) != FMM_IN) update (i1+1,i2,i3, pt, pt+n23, pm+n23, slow[i+n23],slowv[i+n23],eta[i+n23]); 
	    if (i1 > 0    && *(pm-n23) != FMM_IN) update (i1-1,i2,i3, pt, pt-n23, pm-n23, slow[i-n23],slowv[i-n23],eta[i-n23]);
	    if (i2 < n2-1 && *(pm+ n3) != FMM_IN) update (i1,i2+1,i3, pt, pt+ n3, pm+ n3, slow[i+ n3],slowv[i+ n3],eta[i+ n3]); 
	    if (i2 > 0    && *(pm- n3) != FMM_IN) update (i1,i2-1,i3, pt, pt- n3, pm- n3, slow[i- n3],slowv[i- n3],eta[i- n3]);
	    if (i3 < n3-1 && *(pm+  1) != FMM_IN) update (i1,i2,i3+1, pt, pt+  1, pm+  1, slow[i+  1],slowv[i+  1],eta[i+  1]); 
	    if (i3 > 0    && *(pm-  1) != FMM_IN) update (i1,i2,i3-1, pt, pt-  1, pm-  1, slow[i-  1],slowv[i-  1],eta[i-  1]);
	}
    } else if (order==2) {
	/* start marching */
	while (nm > 0) {  /* "far away" points */
	    pt = sf_pqueue_extract ();
	    i = pt-ttime;
	    i3 = i%n3;
	    i2 = (i/n3)%n2;
	    i1 = i/n23;
	    sf_warning("i1=%d i2=%d i3=%d i=%d time=%f slow=%f slow=%f eta=%f nm=%d",i1,i2,i3,i,ttime[i],slow[i],slowv[i],eta[i],nm);
      
	    *(pm = mask+i) = FMM_IN;
	    if (i1 < n1-1 && *(pm+n23) != FMM_IN) update2 (i1+1,i2,i3, pt+n23, pm+n23, slow[i+n23],slowv[i+n23],eta[i+n23]); 
	    if (i1 > 0    && *(pm-n23) != FMM_IN) update2 (i1-1,i2,i3, pt-n23, pm-n23, slow[i-n23],slowv[i-n23],eta[i-n23]);
	    if (i2 < n2-1 && *(pm+ n3) != FMM_IN) update2 (i1,i2+1,i3, pt+ n3, pm+ n3, slow[i+ n3],slowv[i+ n3],eta[i+ n3]); 
	    if (i2 > 0    && *(pm- n3) != FMM_IN) update2 (i1,i2-1,i3, pt- n3, pm- n3, slow[i- n3],slowv[i- n3],eta[i- n3]);
	    if (i3 < n3-1 && *(pm+  1) != FMM_IN) update2 (i1,i2,i3+1, pt+  1, pm+  1, slow[i+  1],slowv[i+  1],eta[i+  1]); 
	    if (i3 > 0    && *(pm-  1) != FMM_IN) update2 (i1,i2,i3-1, pt-  1, pm-  1, slow[i-  1],slowv[i-  1],eta[i-  1]);
	}
    } else {
	/* start marching */
	while (nm > 0) {  /* "far away" points */
	    pt = sf_pqueue_extract ();
	    i = pt-ttime;
	    i3 = i%n3;
	    i2 = (i/n3)%n2;
	    i1 = i/n23;
	    sf_warning("i1=%d i2=%d i3=%d i=%d time=%f slow=%f slow=%f eta=%f nm=%d",i1,i2,i3,i,ttime[i],slow[i],slowv[i],eta[i],nm);

	    *(pm = mask+i) = FMM_IN;
	    if (i1 < n1-1 && *(pm+n23) != FMM_IN) update3 (i1+1,i2,i3, pt+n23, pm+n23, slow[i+n23],slowv[i+n23],eta[i+n23]);
	    if (i1 > 0    && *(pm-n23) != FMM_IN) update3 (i1-1,i2,i3, pt-n23, pm-n23, slow[i-n23],slowv[i-n23],eta[i-n23]);
	    if (i2 < n2-1 && *(pm+ n3) != FMM_IN) update3 (i1,i2+1,i3, pt+ n3, pm+ n3, slow[i+ n3],slowv[i+ n3],eta[i+ n3]); 
	    if (i2 > 0    && *(pm- n3) != FMM_IN) update3 (i1,i2-1,i3, pt- n3, pm- n3, slow[i- n3],slowv[i- n3],eta[i- n3]);
	    if (i3 < n3-1 && *(pm+  1) != FMM_IN) update3 (i1,i2,i3+1, pt+  1, pm+  1, slow[i+  1],slowv[i+  1],eta[i+  1]); 
	    if (i3 > 0    && *(pm-  1) != FMM_IN) update3 (i1,i2,i3-1, pt-  1, pm-  1, slow[i-  1],slowv[i-  1],eta[i-  1]);
	}
    }

    for (i=0; i<n; i++) {
	time[i] = ttime[i];
	sf_warning("time=%f dtime=%f",time[i],ttime[i]);
    }
    
    /* end marching */
    sf_pqueue_close ();
    free (mask);
}

static void update (int p1, int p2, int p3, float* tk, float* tj, unsigned char* mj, float s, float rsv, float eta)
{
    float r, b, c, t1=0, t2=0, u, d3r= d3*s*rsv, b2;
    float d3rr,c1,b1,f,g,gg,ggg,gggg,ggggg,ff,f1;
    float A0,A1,A2,m1,m2,m3,tp4,tpA,tp=0.;
    double tp1,tp2,tp3,den=0.;
    unsigned int k, i;
    int getin=0;

    dd[4] = 1./d3r;

    dd[5] = 1./(d1+d3r);
    dd[6] = 1./(d2+d3r);

    dd[7] = 1./(d1+d2+d3r);

    b = c = 0.; i = k = 0; f = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2=0.; tp3=0.; tp4=0.;
    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d1*t1; b += u; c += u*t1; f += d1;
	i ^= 0x01; k |= 0x01;
    }
    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d2*t1; b += u; c += u*t1; f += d2;
	i ^= 0x01; k |= 0x02;
    }
    if(p3==jj3 || getin || p2==jj2)
	sf_warning("p1=%d p2=%d p3=%d n1=%d n2=%d n3=%d",p1,p2,p3,n1,n2,n3);

    if(k){
	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*f-b;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(c-tp*(2.*b-tp*f))*den;
	tp2 = (tp1*(2.*(b-tp*f)-0.5*tp1*f))*den;
	tp3 = -(tp1*tp1*f-tp2*(2.*(b-tp*f)-tp1*f))*den;
	tp4 = -(tp2*f*(2*tp1+0.5*tp2)-tp3*(2.*(b-tp*f)-tp1*f))*den;
    }
    if(p3==jj3 || getin || p2==jj2)
	sf_warning("tp=%f tp1=%f tp2=%f t1=%f t2=%f den=%f",tp,tp1,tp2,t1,t2,den);

    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.; b2 = 0.; f1 = 0.;
    d3rr = -d3*c*rsv;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;}
    if(p3==jj3 || getin || p2==jj2)
	sf_warning("t1=%f t2=%f",t1,t2);
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2;
	u = d3rr*t1; b1 = b+u; c1= c+u*t1; g= b+t1*f; gg = f; ggg=t1*t1*b;
	gggg = t1*t1*f; ggggg = t1*b; ff =f+d3r;
	u = d3r*t1; b += u; c += u*t1; f += d3rr; b2 = b; f1 = f;
	i ^= 0x01; k |= 0x04;

	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*ff-b2;
	/*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,ff,b,d3r,d1,d2,d3);*/
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*d3*(gggg+4.*ggggg);
	m2 = f1-m1+tp*d3*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*d3*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
					  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
								  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den),SF_ABS(tp2));
	tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(-(tp1*tp1*(f1-m1+6*rsv*d3*tp*
						   (g-tp*gg))-2.*tp2*m3+tp1*tp2*ff)*den),SF_ABS(tp3));
	tp4 = SF_SIG(tp4)*SF_MIN(SF_ABS(-(2.*tp1*(tp1*tp1*rsv*d3*(g-2.*tp*gg)+tp2*(f1-m1+6.*rsv*tp*d3*(g-tp*gg)))+
					  tp2*tp2*0.5*ff-2.*tp3*m3+tp1*tp3*ff)*den),SF_ABS(tp4));
	/*tp1 = -(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
	  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den;
	  tp2 = (2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
	  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den;*/
    }

    if (!k) return;
    /*tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(tp1),.1);
      tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS(tp2),1.);
      tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(tp3),1.);*/
    A0=tp; A1=A0+tp1*eta; A2=A1+tp2*eta*eta; 
    /*sf_warning("tp=%f tp1=%f tp2=%f tp3=%f tp4=%f A1=%f A2=%f A3=%f A4=%f den=%f",tp,tp1,tp2,tp3,tp4,A1,A2,A3,A4,den);*/
    if(p3==jj3 || getin || p2==jj2)
	sf_warning("tp=%f tp1=%f tp2=%f A1=%f A2=%f t1=%f t2=%f den=%f",tp,tp1,tp2,A1,A2,t1,t2,den);
    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tpA = (tp*tp1-(tp*tp2-tp1*tp1)*eta)*den;
    } else
	tpA = tp;

    /*if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001){
      den = A0+A2-2*A1;
      den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
      tpA = (A2*A0-A1*A1)*den;
      } else
      tpA = tp;*/

    /*if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001 || SF_ABS(tp3)>0.00000001 || SF_ABS(tp4)>0.00000001){
      den1 = (tp2-eta*tp3)*(eta*tp1*tp3*tp3*tp3+eta*eta*tp2*tp2*tp3*tp4+ 
      tp2*tp2*tp2*(-tp3+eta*tp4)+tp2*tp3*(tp1*tp3-eta*eta*tp3*tp3-2*eta*tp1*tp4));
      den = (SF_ABS(den1)<0.00000001 ? SF_SIG(den1)*100000000. : 1./den1);
      tpB = (tp*den1+eta*(eta*tp2*tp2*tp2*tp2*tp2*(-tp3 +eta*tp4)+ 
      tp1*tp2*tp2*tp2*(tp2-2*eta*tp3)*(-tp3+eta*tp4)+ 
      tp1*tp1*tp3*(tp2-eta*tp3)*(tp2*tp3+eta*tp3*tp3-2*eta*tp2*tp4)))*den;
      } else
      tpB=tp;*/

    /*if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001 || SF_ABS(tp3)>0.00000001 || SF_ABS(tp4)>0.00000001){
      den = A1+A3-2*A2;
      den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
      S1 = (A3*A1-A2*A2)*den;
      den = A2+A4-2*A3;
      den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
      S2 = (A4*A2-A3*A3)*den;
      den = tpA+S2-2*S1;
      den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
      tpB = (S2*tpA-S1*S1)*den;
      } else
      tpB=tp;*/

    tp = tpA;
    if(p3==jj3  || p3==jj3+1 || p3==jj3-1 || getin || p2==jj2 || p2==jj2+1  || p2==jj2-1 || p2==n2-1 
       || p3==n3-1){
	sf_warning("tp=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,den,k,p1,p2,p3);
	return;
    }

    if(tp<0 || tp >10){
	sf_warning("tp=%f tp1=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,tp1,den,k,p1,p2,p3);
	exit(0);
    }
    /*sf_warning("tp=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,den,k,p1,p2,p3);*/
    /*  tp = r + sqrt (fabs(tp*tp + r*r - c*a)); */
    if (tp < *tj && tp> *tk) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}

void updatelll (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float r, b, c, t1=0., t2=0., u, d3r= d3*s*rsv, b2;
    float d3rr,c1,b1,f,g,gg,ggg,gggg,ggggg,ff,f1;
    double tp1,tp2,den;
    unsigned int k, i;

    dd[4] = 1./d3r;

    dd[5] = 1./(d1+d3r);
    dd[6] = 1./(d2+d3r);

    dd[7] = 1./(d1+d2+d3r);

    b = c = 0.; i = k = 0; f = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2=0.;
    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d1*t1; b += u; c += u*t1; f += d1;
	i ^= 0x01; k |= 0x01;
    }
    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d2*t1; b += u; c += u*t1; f += d2;
	i ^= 0x01; k |= 0x02;
    }

    if(k){
	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*f-b;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(c-tp*(2.*b-tp*f))*den;
	tp2 = (tp1*(2.*(b-tp*f)-0.5*tp1*f))*den;
    }

    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.; b2 = 0.; f1 = 0.;
    d3rr = -d3*c*rsv;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;}  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2;
	u = d3rr*t1; b1 = b+u; c1= c+u*t1; g= b+t1*f; gg = f; ggg=t1*t1*b;
	gggg = t1*t1*f; ggggg = t1*b; ff =f+d3r;
	u = d3r*t1; b += u; c += u*t1; f += d3rr; b2 = b; f1 = f;
	i ^= 0x01; k |= 0x04;

	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*ff-b2;
	/*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,ff,b,d3r,d1,d2,d3);*/
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
					  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
								  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den),SF_ABS(tp2));
    }
    if (!k) return;
    /*tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(tp1),.1);
      tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS(tp2),1.);
      tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(tp3),1.);*/
    /*sf_warning("tp=%f tp1=%f tp2=%f A1=%f A2=%f A3=%f den=%f",tp,tp1,tp2,A1,A2,A3,den);*/
    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tp = (tp*tp1-(tp*tp2-tp1*tp1)*eta)*den;
    }

    /*sf_warning("tp=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,den,k,p1,p2,p3);*/
    /*  tp = r + sqrt (fabs(tp*tp + r*r - c*a)); */
    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}

void updateold (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float r, b, c, t1=0., t2=0., u, d3r= d3*s*rsv, b2;
    float d3rr,c1,b1,f,g,gg,ggg,gggg,ggggg,ff,f1;
    double tp1,tp2,den;
    unsigned int k, i;

    dd[4] = 1./d3r;

    dd[5] = 1./(d1+d3r);
    dd[6] = 1./(d2+d3r);

    dd[7] = 1./(d1+d2+d3r);

    b = c = 0.; i = k = 0; f = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2=0.;
    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d1*t1; b += u; c += u*t1; f += d1;
	i ^= 0x01; k |= 0x01;
    }
    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d2*t1; b += u; c += u*t1; f += d2;
	i ^= 0x01; k |= 0x02;
    }

    if(k){
	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*f-b;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(c-tp*(2.*b-tp*f))*den;
	tp2 = (tp1*(2.*(b-tp*f)-0.5*tp1*f))*den;
    }

    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.; b2 = 0.; f1 = 0.;
    d3rr = -d3*c*rsv;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;}  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2;
	u = d3rr*t1; b1 = b+u; c1= c+u*t1; g= b+t1*f; gg = f; ggg=t1*t1*b;
	gggg = t1*t1*f; ggggg = t1*b; ff =f+d3r;
	u = d3r*t1; b += u; c += u*t1; f += d3rr; b2 = b; f1 = f;
	i ^= 0x01; k |= 0x04;

	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*ff-b2;
	/*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,ff,b,d3r,d1,d2,d3);*/
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
					  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
								  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den),SF_ABS(tp2));
    }
    if (!k) return;

    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tp += tp1*tp1*eta*den;
    }

    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}

void updateh (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float r, b, c, t1=0., t2=0., u, d3r= d3*s*rsv, b2;
    float d3rr,c1,b1,f,g,gg,ggg,gggg,ggggg,ff,f1;
    float m1,m2,m3,tpA;
    double tp1,tp2,tp3=0.,tp4=0.,den;
    unsigned int k, i;

    dd[4] = 1./d3r;

    dd[5] = 1./(d1+d3r);
    dd[6] = 1./(d2+d3r);

    dd[7] = 1./(d1+d2+d3r);

    b = c = 0.; i = k = 0; f = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2=0.;
    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d1*t1; b += u; c += u*t1; f += d1;
	i ^= 0x01; k |= 0x01;
    }
    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d2*t1; b += u; c += u*t1; f += d2;
	i ^= 0x01; k |= 0x02;
    }

    if(k){
	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*f-b;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(c-tp*(2.*b-tp*f))*den;
	tp2 = (tp1*(2.*(b-tp*f)-0.5*tp1*f))*den;
	tp3 = -(tp1*tp1*f-tp2*(2.*(b-tp*f)-tp1*f))*den;
	tp4 = -(tp2*f*(2*tp1+0.5*tp2)-tp3*(2.*(b-tp*f)-tp1*f))*den;
    }

    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.; b2 = 0.; f1 = 0.;
    d3rr = -d3*c*rsv;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;}  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2;
	u = d3rr*t1; b1 = b+u; c1= c+u*t1; g= b+t1*f; gg = f; ggg=t1*t1*b;
	gggg = t1*t1*f; ggggg = t1*b; ff =f+d3r;
	u = d3r*t1; b += u; c += u*t1; f += d3rr; b2 = b; f1 = f;
	i ^= 0x01; k |= 0x04;

	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*ff-b2;
	/*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,ff,b,d3r,d1,d2,d3);*/
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*d3*(gggg+4.*ggggg);
	m2 = f1-m1+tp*d3*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*d3*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
					  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
								  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den),SF_ABS(tp2));
	tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(-(tp1*tp1*(f1-m1+6*rsv*d3*tp*
						   (g-tp*gg))-2.*tp2*m3+tp1*tp2*ff)*den),SF_ABS(tp3));
	tp4 = SF_SIG(tp4)*SF_MIN(SF_ABS(-(2.*tp1*(tp1*tp1*rsv*d3*(g-2.*tp*gg)+tp2*(f1-m1+6.*rsv*tp*d3*(g-tp*gg)))+
					  tp2*tp2*0.5*ff-2.*tp3*m3+tp1*tp3*ff)*den),SF_ABS(tp4));
    }
    if (!k) return;

    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001 && SF_ABS(tp3)>0.0000001 && SF_ABS(tp4)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tpA = tp+tp1*tp1*eta*den;
	den = (tp2-eta*tp3)*(eta*tp1*tp3*tp3*tp3+eta*eta*tp2*tp2*tp3*tp4+ 
			     tp2*tp2*tp2*(-tp3+eta*tp4)+tp2*tp3*(tp1*tp3-eta*eta*tp3*tp3-2*eta*tp1*tp4));
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tp += eta*(eta*tp2*tp2*tp2*tp2*tp2*(-tp3 +eta*tp4)+ 
		   tp1*tp2*tp2*tp2*(tp2-2*eta*tp3)*(-tp3+eta*tp4)+ 
		   tp1*tp1*tp3*(tp2-eta*tp3)*(tp2*tp3+eta*tp3*tp3-2*eta*tp2*tp4))*den;
	sf_warning("tp=%f tpA=%f",tp,tpA);
    }

    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}



static void update2 (int p1, int p2, int p3, float* tj, unsigned char* mj, float s, float rsv, float eta)
{
    float t1=0., t2=0., u, d3r=d3*s*rsv;
    float d3rr,c1,b1,g,gg,ggg,gggg,ggggg,ff;
    float A0,A1,A2,A3,A4,aaa1,bbb1,m1,m2,m3,S1,S2;
    double tp1,tp2,tp3,den=0.,tp4,tpA;
    unsigned int k, i;
    double aaa,bbb,ccc,ddd,ddd1;
    unsigned int jjj;

    i = k = 0;
    ddd=1.; jjj = 0;
    aaa = bbb = ccc = 0.; bbb1 = 0.; aaa1 = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2 = 0.; tp3 = 0.; tp4 = 0.;

    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d1;
	if ((jjj & 0x01) && (p1 < n1-2) && *(mj+2*n23)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n23); t1 /= 3.;
	} else if ((p1 > 1) && *(mj-2*n23)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n23); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x01;
    }
    jjj =0;

    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d2;
	if ((jjj & 0x01) && (p2 < n2-2) && *(mj+2*n3)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n3); t1 /= 3.; 
	} else if ((p2 > 1) && *(mj-2*n3)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n3); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x02;
    }
    jjj =0;
    if(p2==jj2)
	sf_warning("p1=%d p2=%d p3=%d n1=%d n2=%d n3=%d",p1,p2,p3,n1,n2,n3);

    if(k){
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*(ccc-s),0.)))/aaa;
	den = tp*aaa-bbb;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(ccc-tp*(2.*bbb-tp*aaa))*den;
	tp2 = (tp1*(2.*(bbb-tp*aaa)-0.5*tp1*aaa))*den;
	tp3 = -(tp1*tp1*aaa-tp2*(2.*(bbb-tp*aaa)-tp1*aaa))*den;
	tp4 = -(tp2*aaa*(2*tp1+0.5*tp2)-tp3*(2.*(bbb-tp*aaa)-tp1*aaa))*den;
    }
    if(p2==jj2)
	sf_warning("tp=%f tp1=%f tp2=%f v=%f vv=%f eta=%f t1=%f den=%f",tp,tp1,tp2,1/sqrt(s),sqrt(rsv),eta,t1,den);

    ddd1=d3;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; }  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2; 
	ddd=d3r;
	if ((jjj & 0x01) && (p3 < n3-2) && *(mj+2)) { 
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj+2); t1 /= 3.; 
	} else if ((p3 > 1) && *(mj-2)) {   
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj-2); t1 /= 3.;
	}
	d3rr = -ddd1*ccc*rsv;
	u = d3rr*t1; b1 = bbb+u; c1 = ccc+u*t1;
	g = bbb+t1*aaa; gg = aaa; ggg = bbb*t1*t1;
	gggg = aaa*t1*t1; ggggg = t1*bbb; ff = aaa+d3rr;
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1; bbb1 = bbb; aaa1 = aaa;
	i ^= 0x01; k |= 0x04;

	ccc += -s;
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*ccc,0.)))/aaa;
	den = tp*aaa1-bbb1;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*ddd1*(gggg+4.*ggggg);
	m2 = ff-m1+tp*ddd1*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*ddd1*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*ff+tp*tp*tp*ddd1*rsv*
					  (2.*g-tp*gg)-tp*tp*m1+rsv*ddd1*tp*2.*ggg)*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*m3-0.5*tp1*tp1*aaa1)*den),SF_ABS(tp2));
	tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(-(tp1*tp1*(ff-m1+6.*rsv*ddd1*tp*
						   (g-tp*gg))-2.*tp2*m3+tp1*tp2*aaa1)*den),SF_ABS(tp3));
	tp4 = SF_SIG(tp4)*SF_MIN(SF_ABS(-(2.*tp1*(tp1*tp1*rsv*ddd1*(g-2.*tp*gg)+tp2*(ff-m1+6.*rsv*tp*ddd1*(g-tp*gg)))+
					  tp2*tp2*0.5*aaa1-2.*tp3*m3+tp1*tp3*aaa1)*den),SF_ABS(tp4));
    }
    jjj =0;

    if (!k) return;
    /*den = tp*aaa1-bbb1;*/
    /*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,aaa,bbb,d3r,d1,d2,d3);*/
    /*tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(tp1),.1);
      tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS(tp2),1.);
      tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(tp3),1.);*/
    A0=tp; A1=A0+tp1*eta; A2=A1+tp2*eta*eta; A3=A2+tp3*eta*eta*eta; A4=A3+tp4*eta*eta*eta*eta;
    if(p2==jj2)
	sf_warning("tp=%f tp1=%f tp2=%f A1=%f A2=%f t1=%f den=%f",tp,tp1,tp2,tp4,A1,A2,t1,den);
    if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000000001 ? SF_SIG(den)*1000000000000. : 1./den);
	tpA = (tp*tp1-(tp*tp2-tp1*tp1)*eta)*den;
    } else
	tpA = tp;

    /*if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001){
      den = A0+A2-2*A1;
      den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
      tpA = (A2*A0-A1*A1)*den;
      } else
      tpA = tp;*/

    /*if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001 || SF_ABS(tp3)>0.00000001 || SF_ABS(tp4)>0.00000001){
      den1 = (tp2-eta*tp3)*(eta*tp1*tp3*tp3*tp3+eta*eta*tp2*tp2*tp3*tp4+ 
      tp2*tp2*tp2*(-tp3+eta*tp4)+tp2*tp3*(tp1*tp3-eta*eta*tp3*tp3-2*eta*tp1*tp4));
      den = (SF_ABS(den1)<0.000000000001 ? SF_SIG(den1)*1000000000000. : 1./den1);
      tp = (tp*den1+eta*(eta*tp2*tp2*tp2*tp2*tp2*(-tp3 +eta*tp4)+ 
      tp1*tp2*tp2*tp2*(tp2-2*eta*tp3)*(-tp3+eta*tp4)+ 
      tp1*tp1*tp3*(tp2-eta*tp3)*(tp2*tp3+eta*tp3*tp3-2*eta*tp2*tp4)))*den;
      }*/
    if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001 || SF_ABS(tp3)>0.00000001 || SF_ABS(tp4)>0.00000001){
	den = A1+A3-2*A2;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	S1 = (A3*A1-A2*A2)*den;
	den = A2+A4-2*A3;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	S2 = (A4*A2-A3*A3)*den;
	den = tpA+S2-2*S1;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
    } else
    tp = tpA;
    if(p2==jj2)
	sf_warning("tp=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,den,k,p1,p2,p3);
    if(tp<0 || tp >10){
	sf_warning("tp=%f tp1=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,tp1,den,k,p1,p2,p3);
	exit(0);
    }
 
    /*  tp = r + sqrt (fabs(tp*tp + r*r - c*a)); */
    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}


void update2old (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float t1=0., t2=0., u, d3r=d3*s*rsv;
    float d3rr,c1,b1,g,gg,ggg,gggg,ggggg,ff;
    float aaa1,bbb1,m1,m2,m3;
    double tp1,tp2,den,tpp1,tpp2,delta=0;
    unsigned int k, i;
    double aaa,bbb,ccc,ddd,ddd1;
    unsigned int jjj;

    i = k = 0;
    ddd=1.; jjj = 0;
    aaa = bbb = ccc = 0.; bbb1 = 0.; aaa1 = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2 = 0.;

    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d1;
	if ((jjj & 0x01) && (p1 < n1-2) && *(mj+2*n23)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n23); t1 /= 3.;
	} else if ((p1 > 1) && *(mj-2*n23)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n23); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x01;
    }
    jjj =0;

    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d2;
	if ((jjj & 0x01) && (p2 < n2-2) && *(mj+2*n3)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n3); t1 /= 3.; 
	} else if ((p2 > 1) && *(mj-2*n3)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n3); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x02;
    }
    jjj =0;

    if(k){
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*(ccc-s),0.)))/aaa;
	den = tp*aaa-bbb;
	if(SF_ABS(den)<0.00000001)
	    sf_warning("p1=%d p2=%d p3=%d tp=%f den=%f",p1,p2,p3,tp,den);
	den = (SF_ABS(den)<0.00000001 ? SF_SIG(den)*100000000. : 1./den);
	tp1 = -(ccc-tp*(2.*bbb-tp*aaa))*den;
	tp2 = (tp1*(2.*(bbb-tp*aaa)-0.5*tp1*aaa))*den;

	if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001){
	    den = tp1-tp2*eta;
	    den = (SF_ABS(den)<0.00000000000001 ? SF_SIG(den)*1000000000000000. : 1./den);
	    delta = tp1*tp1*eta*den;
	}
    }

    ddd1=d3;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; }  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2; 
	ddd=d3r;
	if ((jjj & 0x01) && (p3 < n3-2) && *(mj+2)) { 
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj+2); t1 /= 3.; 
	} else if ((p3 > 1) && *(mj-2)) {   
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj-2); t1 /= 3.;
	}
	d3rr = -ddd1*ccc*rsv;
	u = d3rr*t1; b1 = bbb+u; c1 = ccc+u*t1;
	g = bbb+t1*aaa; gg = aaa; ggg = bbb*t1*t1;
	gggg = aaa*t1*t1; ggggg = t1*bbb; ff = aaa+d3rr;
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1; bbb1 = bbb; aaa1 = aaa;
	i ^= 0x01; k |= 0x04;

	ccc += -s;
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*ccc,0.)))/aaa;
	den = tp*aaa1-bbb1;
	if(SF_ABS(den)<0.00000001)
	    sf_warning("p1=%d p2=%d p3=%d tp=%f den=%f",p1,p2,p3,tp,den);
	den = (SF_ABS(den)<0.00000001 ? SF_SIG(den)*100000000. : 1./den);
	m1 = rsv*ddd1*(gggg+4.*ggggg);
	m2 = ff-m1+tp*ddd1*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*ddd1*ggg;
	tpp1 = -(c1-2.*tp*b1+tp*tp*ff+tp*tp*tp*ddd1*rsv*
		 (2.*g-tp*gg)-tp*tp*m1+rsv*ddd1*tp*2.*ggg)*den;
	tpp2 = (2.*tp1*m3-0.5*tp1*tp1*aaa1)*den;
	/*if(SF_ABS(tpp1)>SF_ABS(tp1) || SF_ABS(tpp2)>SF_ABS(tp2))
	  sf_warning("tp1=%f tp1new=%f tp2=%f tp2new=%f p1=%d p2=%d p3=%d",tp1,tpp1,tp2,tpp2,p1,p2,p3);*/
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*ff+tp*tp*tp*ddd1*rsv*
					  (2.*g-tp*gg)-tp*tp*m1+rsv*ddd1*tp*2.*ggg)*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*m3-0.5*tp1*tp1*aaa1)*den),SF_ABS(tp2));
	tp1 = tpp1;
	tp2 = tpp2;
    }
    jjj =0;

    if (!k) return;
  
    if(SF_ABS(tp1)>0.000000000001 || SF_ABS(tp2)>0.000000000001){
	den = tp1-tp2*eta;
	if(SF_ABS(den)<0.000000000000001)
	    sf_warning("p1=%d p2=%d p3=%d tp=%f tp1=%f tp2=%f den=%f",p1,p2,p3,tp,tp1,tp2,den);
	den = (SF_ABS(den)<0.000000000000001 ? SF_SIG(den)*1000000000000000. : 1./den);
	tp += SF_SIG(delta)*SF_MIN(SF_ABS(tp1*tp1*eta*den),SF_ABS(delta));
    }

    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}

void update2h (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float t1=0., t2=0., u, d3r=d3*s*rsv;
    float d3rr,c1,b1,g,gg,ggg,gggg,ggggg,ff;
    float aaa1,bbb1,m1,m2,m3;
    double tp1,tp2,tp3=0.,tp4=0.,den;
    unsigned int k, i;
    double aaa,bbb,ccc,ddd,ddd1;
    unsigned int jjj;

    i = k = 0;
    ddd=1.; jjj = 0;
    aaa = bbb = ccc = 0.; bbb1 = 0.; aaa1 = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2 = 0.;

    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d1;
	if ((jjj & 0x01) && (p1 < n1-2) && *(mj+2*n23)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n23); t1 /= 3.;
	} else if ((p1 > 1) && *(mj-2*n23)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n23); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x01;
    }
    jjj =0;

    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d2;
	if ((jjj & 0x01) && (p2 < n2-2) && *(mj+2*n3)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n3); t1 /= 3.; 
	} else if ((p2 > 1) && *(mj-2*n3)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n3); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x02;
    }
    jjj =0;

    if(k){
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*(ccc-s),0.)))/aaa;
	den = tp*aaa-bbb;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(ccc-tp*(2.*bbb-tp*aaa))*den;
	tp2 = (tp1*(2.*(bbb-tp*aaa)-0.5*tp1*aaa))*den;
	tp3 = -(tp1*tp1*aaa-tp2*(2.*(bbb-tp*aaa)-tp1*aaa))*den;
	tp4 = -(tp2*aaa*(2*tp1+0.5*tp2)-tp3*(2.*(bbb-tp*aaa)-tp1*aaa))*den;
    }

    ddd1=d3;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; }  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2; 
	ddd=d3r;
	if ((jjj & 0x01) && (p3 < n3-2) && *(mj+2)) { 
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj+2); t1 /= 3.; 
	} else if ((p3 > 1) && *(mj-2)) {   
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj-2); t1 /= 3.;
	}
	d3rr = -ddd1*ccc*rsv;
	u = d3rr*t1; b1 = bbb+u; c1 = ccc+u*t1;
	g = bbb+t1*aaa; gg = aaa; ggg = bbb*t1*t1;
	gggg = aaa*t1*t1; ggggg = t1*bbb; ff = aaa+d3rr;
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1; bbb1 = bbb; aaa1 = aaa;
	i ^= 0x01; k |= 0x04;

	ccc += -s;
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*ccc,0.)))/aaa;
	den = tp*aaa1-bbb1;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*ddd1*(gggg+4.*ggggg);
	m2 = ff-m1+tp*ddd1*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*ddd1*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*ff+tp*tp*tp*ddd1*rsv*
					  (2.*g-tp*gg)-tp*tp*m1+rsv*ddd1*tp*2.*ggg)*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*m3-0.5*tp1*tp1*aaa1)*den),SF_ABS(tp2));
	tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(-(tp1*tp1*(ff-m1+6.*rsv*ddd1*tp*
						   (g-tp*gg))-2.*tp2*m3+tp1*tp2*aaa1)*den),SF_ABS(tp3));
	tp4 = SF_SIG(tp4)*SF_MIN(SF_ABS(-(2.*tp1*(tp1*tp1*rsv*ddd1*(g-2.*tp*gg)+tp2*(ff-m1+6.*rsv*tp*ddd1*(g-tp*gg)))+
					  tp2*tp2*0.5*aaa1-2.*tp3*m3+tp1*tp3*aaa1)*den),SF_ABS(tp4));
    }
    jjj =0;

    if (!k) return;

    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001 && SF_ABS(tp3)>0.0000001 && SF_ABS(tp4)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	den = (tp2-eta*tp3)*(eta*tp1*tp3*tp3*tp3+eta*eta*tp2*tp2*tp3*tp4+ 
			     tp2*tp2*tp2*(-tp3+eta*tp4)+tp2*tp3*(tp1*tp3-eta*eta*tp3*tp3-2*eta*tp1*tp4));
	if(SF_ABS(den)<0.0000000000000001)
	    sf_warning("p1=%d p2=%d p3=%d tp=%f tp1=%f tp2=%f tp3=%f tp4=%f den=%f",p1,p2,p3,tp,tp1,tp2,tp3,tp4,den);
	den = (SF_ABS(den)<0.0000000000000001 ? SF_SIG(den)*10000000000000000. : 1./den);
	tp += eta*(eta*tp2*tp2*tp2*tp2*tp2*(-tp3 +eta*tp4)+ 
		   tp1*tp2*tp2*tp2*(tp2-2*eta*tp3)*(-tp3+eta*tp4)+ 
		   tp1*tp1*tp3*(tp2-eta*tp3)*(tp2*tp3+eta*tp3*tp3-2*eta*tp2*tp4))*den;
    }

    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}

void update3old (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float t1=0., t2=0., u, d3r=d3*s*rsv;
    float d3rr,c1,b1,g,gg,ggg,gggg,ggggg,ff;
    float aaa1,bbb1,m1,m2,m3;
    double tp1,tp2,den;
    unsigned int k, i;
    double aaa,bbb,ccc,ddd,ddd1;
    unsigned int jjj;
    int xxx=1;

    i = k = 0;
    ddd=1.; jjj = 0;
    aaa = bbb = ccc = 0.; bbb1 = 0.; aaa1 = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2 = 0.;

    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d1;
	if (jjj & 0x01) {
	    if   ((p1 < n1-2) && *(mj+2*n23)) {
		if ((p1 < n1-3) && *(mj+3*n23) && xxx) {
		    ddd*=(121./36.); t1*=18.; t1-=*(tj+2*n23)*9.; t1+=*(tj+3*n23)*2.; t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n23); t1 /= 3.;
		}
	    }
	} else {
	    if   ((p1 > 1) && *(mj-2*n23)) {
		if ((p1 > 2) && *(mj-3*n23) && xxx) {
		    ddd*=(121./36.); t1*=18.; t1-=*(tj-2*n23)*9.; t1+=*(tj-3*n23)*2.; t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n23); t1 /= 3.;
		}
	    }
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x01;
    }
    jjj =0;

    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d2;
	if (jjj & 0x01) {
	    if   ((p2 < n2-2) && *(mj+2*n3)) {
		if ((p2 < n2-3) && *(mj+3*n3) && xxx) {
		    ddd*=(121./36.);t1*=18.;t1-=*(tj+2*n3)*9.;t1+=*(tj+3*n3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n3); t1 /= 3.;
		}
	    }
	} else {
	    if   ((p2 > 1) && *(mj-2*n3)) {
		if ((p2 > 2) && *(mj-3*n3) && xxx) {
		    ddd*=(121./36.);t1*=18.;t1-=*(tj-2*n3)*9.;t1+=*(tj-3*n3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n3); t1 /= 3.;
		}
	    }
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x02;
    }
    jjj =0;

    if(k){
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*(ccc-s),0.)))/aaa;
	den = tp*aaa-bbb;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(ccc-tp*(2.*bbb-tp*aaa))*den;
	tp2 = (tp1*(2.*(bbb-tp*aaa)-0.5*tp1*aaa))*den;
    }

    ddd1=d3;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; }  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2; 
	ddd=d3r;
	if (jjj & 0x01) {
	    if   ((p3 < n3-2) && *(mj+2)) {
		if ((p3 < n3-3) && *(mj+3) && xxx) {
		    ddd*=(121./36.); ddd1*=(121./36.);t1*=18.;t1-=*(tj+2)*9.;
		    t1+=*(tj+3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj+2); t1 /= 3.;
		}
	    }
	} else {
	    if   ((p3 > 1) && *(mj-2)) {
		if ((p3 > 2) && *(mj-3) && xxx) {
		    ddd*=(121./36.);ddd1*=(121./36.);t1*=18.;t1-=*(tj-2)*9.;
		    t1+=*(tj-3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj-2); t1 /= 3.;
		}
	    }
	}
	d3rr = -ddd1*ccc*rsv;
	u = d3rr*t1; b1 = bbb+u; c1 = ccc+u*t1;
	g = bbb+t1*aaa; gg = aaa; ggg = bbb*t1*t1;
	gggg = aaa*t1*t1; ggggg = t1*bbb; ff = aaa+d3rr;
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1; bbb1 = bbb; aaa1 = aaa;
	i ^= 0x01; k |= 0x04;

	ccc += -s;
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*ccc,0.)))/aaa;
	den = tp*aaa1-bbb1;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*ddd1*(gggg+4.*ggggg);
	m2 = ff-m1+tp*ddd1*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*ddd1*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*ff+tp*tp*tp*ddd1*rsv*
					  (2.*g-tp*gg)-tp*tp*m1+rsv*ddd1*tp*2.*ggg)*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*m3-0.5*tp1*tp1*aaa1)*den),SF_ABS(tp2));
    }
    jjj =0;

    if (!k) return;

    /*sf_warning("tp=%f tp1=%f tp2=%f den=%f",tp,tp1,tp2,den);*/
  
    if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000000001 ? SF_SIG(den)*1000000000000. : 1./den);
	tp += tp1*tp1*eta*den;
    }
    /*sf_warning("tp=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,den,k,p1,p2,p3);*/
    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}


static void update3 (int p1, int p2, int p3, float* tj, unsigned char* mj, float s, float rsv, float eta)
{
    float t1=0., t2=0., u, d3r=d3*s*rsv;
    float d3rr,c1,b1,g,gg,ggg,gggg,ggggg,ff;
    float aaa1,bbb1,m1,m2,m3;
    double tp1,tp2,den=0.;
    unsigned int k, i;
    double aaa,bbb,ccc,ddd,ddd1;
    unsigned int jjj;
    int xxx=1;

    i = k = 0;
    ddd=1.; jjj = 0;
    aaa = bbb = ccc = 0.; bbb1 = 0.; aaa1 = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2 = 0.;

    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d1;
	if (jjj & 0x01) {
	    if   ((p1 < n1-2) && *(mj+2*n23)) {
		if ((p1 > 0) && *(mj-n23) && ((t2 = *(tj-n23)) < *tj)) {
		    ddd*=0.25; t1*=6.; t1-=*(tj+2*n23); t1+=t2*2.; t1/=3.;
		} else if ((p1 < n1-3) && *(mj+3*n23) && xxx) {
		    ddd*=(121./36.); t1*=18.; t1-=*(tj+2*n23)*9.; t1+=*(tj+3*n23)*2.; t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n23); t1 /= 3.;
		}
	    }
	} else {
	    if   ((p1 > 1) && *(mj-2*n23)) {
		if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj)) {
		    ddd*=0.25; t1*=6.; t1-=*(tj-2*n23); t1+=t2*2.; t1/=3.;
		} else if ((p1 > 2) && *(mj-3*n23) && xxx) {
		    ddd*=(121./36.); t1*=18.; t1-=*(tj-2*n23)*9.; t1+=*(tj-3*n23)*2.; t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n23); t1 /= 3.;
		}
	    }
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x01;
    }
    jjj =0;

    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d2;
	if (jjj & 0x01) {
	    if   ((p2 < n2-2) && *(mj+2*n3)) {
		if ((p2 > 0) && *(mj-n3) && ((t2 = *(tj-n3)) < *tj)) {
		    ddd*=0.25; t1*=6.; t1-=*(tj+2*n3); t1+=t2*2.; t1/=3.;
		} else if ((p2 < n2-3) && *(mj+3*n3) && xxx) {
		    ddd*=(121./36.);t1*=18.;t1-=*(tj+2*n3)*9.;t1+=*(tj+3*n3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n3); t1 /= 3.;
		}
	    }
	} else {
	    if   ((p2 > 1) && *(mj-2*n3)) {
		if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj)) {
		    ddd*=0.25; t1*=6.; t1-=*(tj-2*n3); t1+=t2*2.; t1/=3.;
		} else if ((p2 > 2) && *(mj-3*n3) && xxx) {
		    ddd*=(121./36.);t1*=18.;t1-=*(tj-2*n3)*9.;t1+=*(tj-3*n3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n3); t1 /= 3.;
		}
	    }
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x02;
    }
    jjj =0;

    if(k){
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*(ccc-s),0.)))/aaa;
	den = tp*aaa-bbb;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(ccc-tp*(2.*bbb-tp*aaa))*den;
	tp2 = (tp1*(2.*(bbb-tp*aaa)-0.5*tp1*aaa))*den;
    }

    ddd1=d3;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; }  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2; 
	ddd=d3r;
	if (jjj & 0x01) {
	    if   ((p3 < n3-2) && *(mj+2)) {
		if ((p3 > 0) && *(mj-1) && ((t2 = *(tj-1)) < *tj)) {
		    ddd*=0.25; t1*=6.; t1-=*(tj+2); t1+=t2*2.; t1/=3.;
		} else if ((p3 < n3-3) && *(mj+3) && xxx) {
		    ddd*=(121./36.); ddd1*=(121./36.);t1*=18.;t1-=*(tj+2)*9.;
		    t1+=*(tj+3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj+2); t1 /= 3.;
		}
	    }
	} else {
	    if   ((p3 > 1) && *(mj-2)) {
		if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj)) {
		    ddd*=0.25; t1*=6.; t1-=*(tj-2); t1+=t2*2.; t1/=3.;
		} else if ((p3 > 2) && *(mj-3) && xxx) {
		    ddd*=(121./36.);ddd1*=(121./36.);t1*=18.;t1-=*(tj-2)*9.;
		    t1+=*(tj-3)*2.;t1/=11.;
		} else {
		    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj-2); t1 /= 3.;
		}
	    }
	}
	d3rr = -ddd1*ccc*rsv;
	u = d3rr*t1; b1 = bbb+u; c1 = ccc+u*t1;
	g = bbb+t1*aaa; gg = aaa; ggg = bbb*t1*t1;
	gggg = aaa*t1*t1; ggggg = t1*bbb; ff = aaa+d3rr;
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1; bbb1 = bbb; aaa1 = aaa;
	i ^= 0x01; k |= 0x04;

	ccc += -s;
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*ccc,0.)))/aaa;
	den = tp*aaa1-bbb1;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*ddd1*(gggg+4.*ggggg);
	m2 = ff-m1+tp*ddd1*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*ddd1*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*ff+tp*tp*tp*ddd1*rsv*
					  (2.*g-tp*gg)-tp*tp*m1+rsv*ddd1*tp*2.*ggg)*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*m3-0.5*tp1*tp1*aaa1)*den),SF_ABS(tp2));
    }
    jjj =0;

    if (!k) return;

    /*sf_warning("tp=%f tp1=%f tp2=%f den=%f",tp,tp1,tp2,den);*/
  
    if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000000001 ? SF_SIG(den)*1000000000000. : 1./den);
	tp += tp1*tp1*eta*den;
    }
    if(tp<0 || tp >10){
	sf_warning("tp=%f tp1=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,tp1,den,k,p1,p2,p3);
	exit(0);
    }
    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}
