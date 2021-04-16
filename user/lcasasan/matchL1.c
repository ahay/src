/* Matched filter L1 estimation */
/*
  Copyright (C) 2010 Politecnico di Milano
  
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
#include "helix_icaf.h"

#include "matchL1.h"


static float *mm;
static sf_filter aa;

static void fit(bool adj, bool add, int nm, int nd, float *m, float *d)
/* L1 fitting criterion */
{

    sf_adjnull(adj, add, nm, nd, m, d);
    sf_copy_lop(adj, true, nd, nd, m, d);

    helix_icaf_init(mm, aa, 1);

    helix_icaf_lop(adj,true,aa->nh,nd,m+nd,d);
}


void find_matchL1(int nd   /* data size */,
	      float* mm1       /* match [nd] */,
		  float* dd        /* data [nd] */,
	      sf_filter aa1    /* estimated filter */,
	      int liter        /* number of CG iterations */,
	      int niter        /* number of POCS iterations */,
	      float perc	   /* percentage for sharpening */,
	      float eps        /* dumping parameter */)
/*< find L1 matched filter >*/
{
	float *model;
    int i,iter,nmodel;

    aa=aa1;
	mm=mm1;
	nmodel = nd+aa->nh;

	model = sf_floatalloc(nmodel);

    for (i=0;i<nd;i++) { // this is a slight dump to avoid 0 = aa->flt * 0
    	mm[i]+=eps;
    	dd[i]+=eps;
    }

    for (i=0; i < nmodel; i++)
    	model[i]=0.0;

    sf_sharpen_init(nd,perc,0.5);
    /* initialization of sharpening regularization*/


    sf_sharpen_init(nd,perc,0.5);

    for (iter=0; iter < niter; iter++) {
    	sf_solver(fit,sf_cgstep,nmodel,nd,model,dd,liter,"x0",model,"verb",false,"end");
    	sf_sharpen(model);
    	sf_weight_apply(nd,model);
    /* apply sharpening regularization*/
    }

    //aa1->flt=model+nd; /* this dont work

    for (i=0;i<aa1->nh;i++)
    	aa1->flt[i]=model[nd+i];

    free(model);
}



//void find_matchL1(int nd   /* data size */,
//	      float* m      /* match [nd] */,
//		  float* d        /* data [nd] */,
//	      sf_filter aa    /* estimated filter */,
//	      int liter        /* number of CG iterations */,
//	      int niter        /* number of POCS iterations */,
//	      float perc	   /* percentage for sharpening */,
//	      float eps        /* dumping parameter */)
///*< find L1 matched filter >*/
//{
//	float  *dd,*n,*r;
//	float  *aa_tmp,ea,en;
//    int ia,id,iter;
//
//    aa_tmp=sf_floatalloc(aa->nh);
//
//    /* set the helix filter coefficients to zero */
//    for (ia=0; ia < aa->nh; ia++) {
//    	aa_tmp[ia] = aa->flt[ia]=0.0;
//    }
//
//
//    dd = sf_floatalloc(nd);
//    n = sf_floatalloc(nd);
//    r = sf_floatalloc(nd);
//
//    for (id=0;id<nd;id++) {
//    	// this is a slight dump to avoid 0 = aa->flt * 0
//    	m[id]+=eps;
//    	d[id]+=eps;
//        /* set auxiliary data size vector to zero */
//    	dd[id]=n[id]=r[id]=0.0;
//    }
//
//    sf_sharpen_init(nd,perc);
//     /* initialization of sharpening regularization*/
//
//    helix_icaf_init(m, aa, 1);
//
//    for (iter=0; iter < niter; iter++) {
//		/* Solve | (d - n) - M * aa |_2 */
//		/* -------------------------------------- */
//    	ea=en=0;
//    	for (id=0; id < nd; id++)
//    		r[id] = dd[id] = d[id]-n[id]; /* dd = d - n */
//
//		for (ia=0; ia < aa->nh; ia++)
//			aa_tmp[ia] = (-1) *aa->flt[ia];  /* -aa */
//
//	    helix_icaf_lop(false,true,aa->nh,nd,aa_tmp,r); /* r = dd - M * aa;*/
//
//    	sf_solver(helix_icaf_lop,sf_cgstep,aa->nh,nd,aa->flt,dd,liter,"x0",aa->flt,"verb",false,"end");
//
//    	for (ia=0; ia < aa->nh; ia++) {
//    		aa_tmp[ia]= -aa_tmp[ia] - aa->flt[ia];  /* -da */
//    		ea+=aa->flt[ia]*aa->flt[ia];
//    	}
//
//    	for (id=0; id < nd; id++)
//    		n[id] += r[id];
//
//	    helix_icaf_lop(false,true,aa->nh,nd,aa_tmp,n); /* n += r - M * da;*/
//
//	    for (id=0; id < nd; id++)
//	    	en+= n[id]*n[id];
//
//	    /* apply sharpening regularization*/
//	    sf_sharpen(n);
//	    sf_weight_apply(nd,n);
//
//    }
//    free(dd);
//    free(n);
//    free(r);
//    free(aa_tmp);
//}

