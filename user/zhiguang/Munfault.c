/* Unfault image */
/*
 Copyright (C) 2017 University of Texas at Austin
 
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
#include <rsfpwd.h>

int i1, i2, n1, n2, n12, ig, ng, i, j, k;
int *np, ***fxy, ***sxy;

void sf_fault_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
	int ii, jj;

	sf_adjnull(adj, add, nx, ny, x, y);

	if(adj){
		for(i1=0; i1<n1; i1++)
			x[i1] += y[i1];

		k=n1;
		for(ig=0; ig<ng; ig++){
			for(i1=0; i1<np[ig]; i1++){
				i =fxy[ig][i1][0];
				j =fxy[ig][i1][1];
				ii=i+sxy[ig][i1][0];
				jj=j+sxy[ig][i1][1];

				x[j*n1+i] += y[k];
				x[jj*n1+ii] += -y[k];
				k++;
			}
		}
	}else{
		for(i1=0; i1<n1; i1++)
			y[i1] += x[i1];

		k=n1;
		for(ig=0; ig<ng; ig++){
			for(i1=0; i1<np[ig]; i1++){
				i =fxy[ig][i1][0];
				j =fxy[ig][i1][1];
				ii=i+sxy[ig][i1][0];
				jj=j+sxy[ig][i1][1];

				y[k] += x[j*n1+i]-x[jj*n1+ii];
				k++;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	bool inv, mode;
	int niter, order, ns, nknown;
	float **dd, **ss, **mm, **bb, **ff;
	float *xx, *yy, *p, **pp, eps, lam;
	sf_file in, out, slip, mask=NULL, shift=NULL, dip=NULL;

	sf_init(argc, argv);

	if(!sf_getbool("inv", &inv)) inv=false;
	if(!sf_getbool("mode", &mode)) mode=true;

	in=sf_input("in");
	slip=sf_input("slip");
	out=sf_output("out");
	if(inv){
		shift=sf_input("shift");
		dip=sf_input("dip");
	}else mask=sf_output("mask");

	if(!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if(!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	n12=n1*n2;

	if(!sf_getint("niter", &niter)) niter=100;
	/* number of iterations */
	if(!sf_getint("order", &order)) order=1;
	/* accuracy order */
	if(!sf_getint("ns", &ns)) ns=1;
	/* smoothing radius */
	if(!sf_getfloat("eps", &eps)) eps=0.01;
	/* regularization */

	dd=sf_floatalloc2(n1, n2);
	ss=sf_floatalloc2(n1, n2);
	mm=sf_floatalloc2(n1, n2);

	sf_floatread(dd[0], n12, in);
	sf_floatread(ss[0], n12, slip);
	memset(mm[0], 0, n12*sizeof(float));

	/* search fault */
	ng=0;
	if(mode){
		for (i1=1; i1<n1; i1++){
			for (i2=1; i2<n2; i2++){
				if(ss[i2][i1]>-100.) {
					if(ss[i2-1][i1-1]>-100) mm[i2][i1]=mm[i2-1][i1-1]; 
					else if (ss[i2][i1-1]>-100) mm[i2][i1]=mm[i2][i1-1]; 
					else if (ss[i2+1][i1-1]>-100) mm[i2][i1]=mm[i2+1][i1-1]; 
					else if (ss[i2-1][i1]>-100) mm[i2][i1]=mm[i2-1][i1]; 
					else {
						mm[i2][i1]=ng+1;
						ng++;
					}
				} // slip!=0 
			} // end if i1
		} // end of i2
	}else{ // mode=false
		for (i2=1; i2<n2; i2++){
			for (i1=1; i1<n1; i1++){
				if(ss[i2][i1]>-100.) {
					if(ss[i2-1][i1-1]>-100) mm[i2][i1]=mm[i2-1][i1-1]; 
					else if (ss[i2-1][i1]>-100) mm[i2][i1]=mm[i2-1][i1]; 
					else if (ss[i2-1][i1+1]>-100) mm[i2][i1]=mm[i2-1][i1+1]; 
					else if (ss[i2][i1-1]>-100) mm[i2][i1]=mm[i2][i1-1]; 
					else {
						mm[i2][i1]=ng+1;
						ng++;
					}
				} // slip!=0 
			} // end if i1
		} // end of i2
	} // end of mode

	/* extract fault boundary, position, number of sample */
	np=sf_intalloc(ng);
	bb=sf_floatalloc2(n1, 2*ng);
	fxy=sf_intalloc3(2, n1, ng);
	memset(bb[0], 0, 2*ng*n1*sizeof(float));

	for (ig=0; ig<ng; ig++){
		i=0;
		for (i1=0; i1<n1; i1++){
			for (i2=0; i2<n2; i2++){
				if(mm[i2][i1]==ig+1){
					bb[ig][i]=dd[i2-1][i1];
					bb[ig+ng][i]=dd[i2+1][i1];
					fxy[ig][i][0]=i1;
					fxy[ig][i][1]=i2;
					i++;
				} // one fault
			} // end of i2
		} // end of i1
		np[ig]=i;
	} // end of ig

	if(!inv){
		sf_floatwrite(mm[0], n12, mask);
		sf_putint(out, "n2", 2*ng);
		sf_floatwrite(bb[0], n1*2*ng, out);
		exit(0);
	}

	/* read fault slip calculated using local similarity */
	ff =sf_floatalloc2(n1, ng);
	sxy=sf_intalloc3(2, n1, ng);
	memset(sxy[0][0], 0, 2*ng*n1*sizeof(int));

	sf_floatread(ff[0], ng*n1, shift);
	for(ig=0; ig<ng; ig++){
		for(i1=0; i1<np[ig]; i1++){
			sxy[ig][i1][0]=ff[ig][i1]+0.5;

			i=i1+sxy[ig][i1][0]; 
			if(i<0) i=0;
			if(i>np[ig]-1) i=np[ig]-1;
			sxy[ig][i1][1]=fxy[ig][i][1]-fxy[ig][i1][1];
			//if(ig==3) sf_warning("ig=%d i1=%d np[ig]=%d fz=%d fx=%d, sz=%d sx=%d",
			//ig, i1, np[ig], fxy[ig][i1][0], fxy[ig][i1][1], sxy[ig][i1][0], sxy[ig][i1][1]);
		}
	}

	/* read dip and initialize pwsmooth */
	pp=sf_floatalloc2(n1, n2);
	sf_floatread(pp[0], n12, dip);
	pwsmooth_init(ns, n1, n2, order, eps);

	/* figure out scaling */
	nknown=0;
	for(ig=0; ig<ng; ig++) nknown+=np[ig];
	lam = sqrtf(1.*nknown/n12);
	nknown+=n1;

	/* 1. prepare the yy in Ax=yy */
	xx=sf_floatalloc(n12);
	yy=sf_floatalloc(nknown);
	p=sf_floatalloc(n12);
	memset(yy, 0, nknown*sizeof(float));

	k=n1;
	for (ig=0; ig<ng; ig++){
		for (i1=0; i1<np[ig]; i1++){
			yy[k]=sxy[ig][i1][1];
			k++;
		}
	}

	/* 1. solve shift in z */
	sf_conjgrad_init(n12, n12, nknown, nknown, lam, 10*FLT_EPSILON, true, false);
	pwsmooth_set(pp);
	sf_conjgrad(NULL, sf_fault_lop, pwsmooth_lop, p, xx, yy, niter);

	sf_putint(out, "n3", 2);
	sf_floatwrite(xx, n12, out);

	/* 2. prepare the yy in Ax=yy */
	k=n1;
	for (ig=0; ig<ng; ig++){
		for (i1=0; i1<np[ig]; i1++){
			yy[k]=sxy[ig][i1][0];
			k++;
		}
	}

	sf_conjgrad(NULL, sf_fault_lop, pwsmooth_lop, p, xx, yy, niter);
	sf_conjgrad_close();

	sf_floatwrite(xx, n12, out);

	exit(0);
}
