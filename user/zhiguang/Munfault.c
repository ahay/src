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

int i1, i2, n1, n2, n12, ig, ng, i, j, off, ii, jj, k;
int *np, ***fxy, ***sxy; 
float u1, u2, ***u12, **dx, **dz, **dxx, **dzz;

void lin_int2(float **xx, float **mm)
/*< linear interpolation >*/
{
	int k1, k2, pos[ng+2];
	float sum, **tmp;

	tmp=sf_floatalloc2(n1, n2);
	memset(tmp[0], 0, n12*sizeof(float));

	// 1. smoothing fault value
	for(ig=0; ig<ng; ig++){
		for(i1=0; i1<np[ig]; i1++){
			i=fxy[ig][i1][0];
			j=fxy[ig][i1][1];
			for(k=j-off; k<=j+off; k++){
				sum=0.;
				for(k1=-5; k1<=5; k1++){
					k2=i1+k1;
					if(k2<0) k2=0;
					if(k2>np[ig]-1) k2=np[ig]-1;
					ii=fxy[ig][k2][0];
					jj=fxy[ig][k2][1]+k-j;
					sum += xx[jj][ii];
				}
				tmp[k][i]=sum/11.;
			}
		}
	} 

	// 2. linear interpolation
	for(i1=0; i1<n1; i1++){
		k=1;
		pos[0]=0;
		for(i2=0; i2<n2; i2++)
			if(mm[i2][i1]!=0){
				pos[k]=i2;
				k++;
			}
		pos[k]=n2-1;

		for(k1=0; k1<k; k1++){
			ii=pos[k1]+off;
			jj=pos[k1+1]-off;
			for(k2=ii+1; k2<jj; k2++){
				sum=1.*(k2-ii)/(jj-ii);
				tmp[k2][i1]=(1.-sum)*tmp[ii][i1]+sum*tmp[jj][i1];
			}
		}
	}

	// 3. copy to xx
	for(i2=0; i2<n2; i2++) 
		for(i1=0; i1<n1; i1++) 
			xx[i2][i1]=tmp[i2][i1];

	free(*tmp); free(tmp);
}

void sf_grad_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
	float *tmp;

	sf_adjnull(adj, add, nx, ny, x, y);
	if(adj){
		tmp=sf_floatalloc(n12);
		for(i=0; i<n12; i++){
			tmp[i]=x[i];
			x[i]=y[i];
			y[i]=0.;
		}
	}

	// gradient
	for(i2=0; i2<n2-1; i2++)
		for(i1=0; i1<n1; i1++)
			dx[i2][i1] = x[(i2+1)*n1+i1] - x[i2*n1+i1];
	for(i1=0; i1<n1; i1++)
		dx[n2-1][i1] = x[(n2-1)*n1+i1] - x[(n2-2)*n1+i1];

	for(i1=0; i1<n1-1; i1++)
		for(i2=0; i2<n2; i2++)
			dz[i2][i1] = x[i2*n1+i1+1] - x[i2*n1+i1];
	for(i2=0; i2<n2; i2++)
		dz[i2][n1-1] = x[i2*n1+n1-1] - x[i2*n1+n1-2];

	// weighting
	for(ig=0; ig<ng; ig++){
		for(i1=0; i1<np[ig]; i1++){
			i =fxy[ig][i1][0];
			j =fxy[ig][i1][1];
				dx[j][i]=0.;
				dz[j][i]=0.;
		}
	}

	// structure tensor
	for(i2=0; i2<n2; i2++){
		for(i1=0; i1<n1; i1++){
			dzz[i2][i1] = u12[i2][i1][0]*dz[i2][i1]+u12[i2][i1][1]*dx[i2][i1];
			dxx[i2][i1] = u12[i2][i1][1]*dz[i2][i1]+u12[i2][i1][2]*dx[i2][i1];
		}
	}

	// second gradient
	for(i2=1; i2<n2-2; i2++)
		for(i1=0; i1<n1; i1++)
			y[i2*n1+i1] = -dxx[i2][i1] + dxx[i2-1][i1];
	for(i1=0; i1<n1; i1++){
		y[i1] += -dxx[0][i1];
		y[(n2-2)*n1+i1] += -dxx[n2-1][i1] - dxx[n2-2][i1] + dxx[n2-3][i1];
		y[(n2-1)*n1+i1] += dxx[n2-1][i1] + dxx[n2-2][i1];
	}

	for(i1=1; i1<n1-2; i1++)
		for(i2=0; i2<n2; i2++)
			y[i2*n1+i1] += -dzz[i2][i1] + dzz[i2][i1-1];
	for(i2=0; i2<n2; i2++){
		y[i2*n1] += -dzz[i2][0];
		y[i2*n1+n1-2] += -dzz[i2][n1-1] - dzz[i2][n1-2] + dzz[i2][n1-3];
		y[i2*n1+n1-1] += dzz[i2][n1-1] + dzz[i2][n1-2];
	}
	
	if(adj){
		for(i=0; i<n12; i++) x[i]=tmp[i]+y[i];
		free(tmp);
	}
}

void sf_shape_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
	float *tmp;

	sf_adjnull(adj, add, nx, ny, x, y);

	tmp=sf_floatalloc(n12);
	if(adj){
		for(i=0; i<n12; i++) tmp[i]=x[i];
		sf_solver_reg(sf_copy_lop, sf_cgstep, sf_grad_lop, n12, n12, n12, x, y, 30, 1., "end");
		for(i=0; i<n12; i++) x[i] += tmp[i];
	}else{
		for(i=0; i<n12; i++) tmp[i]=y[i];
		sf_solver_reg(sf_copy_lop, sf_cgstep, sf_grad_lop, n12, n12, n12, y, x, 30, 1., "end");
		for(i=0; i<n12; i++) y[i] += tmp[i];
	}
	free(tmp);
}

void sf_fault_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
	sf_adjnull(adj, add, nx, ny, x, y);

	if(adj){
		for(i1=0; i1<n1; i1++)
			x[i1] += y[i1];

		k=n1;
		for(ig=0; ig<ng; ig++){
			for(i1=0; i1<np[ig]; i1++){
				i =fxy[ig][i1][0];
				j =fxy[ig][i1][1]-off;
				ii=i+sxy[ig][i1][0];
				jj=j+sxy[ig][i1][1]+2*off;

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
				j =fxy[ig][i1][1]-off;
				ii=i+sxy[ig][i1][0];
				jj=j+sxy[ig][i1][1]+2*off;

				y[k] += x[j*n1+i]-x[jj*n1+ii];
				k++;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	bool inv, mode;
	int niter, nknown;
	float **dd, **ss, **mm, **bb, **ff;
	float **xx, **zz, *p, *y, **pp, lam, norm;
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
	if(!sf_getint("off", &off)) off=2;
	/* offset to fault */
	if(!sf_getfloat("lam", &lam)) lam=1;
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
		for(i2=0; i2<n2; i2++){ // first line
			if(ss[i2][0]>-100.){
				if(ss[i2-1][0]>-100) mm[i2][0]=mm[i2-1][0];
				else {
					mm[i2][0]=ng+1;
					ng++;
				}
			}
		}
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
					bb[ig][i]=dd[i2-off][i1];
					bb[ig+ng][i]=dd[i2+off][i1];
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

	/* prepare initial value */
	xx=sf_floatalloc2(n1,n2);
	zz=sf_floatalloc2(n1,n2);
	p=sf_floatalloc(n12);
	memset(xx[0], 0, n12*sizeof(float));
	memset(zz[0], 0, n12*sizeof(float));
	for(ig=0; ig<ng; ig++){
		for(i1=0; i1<np[ig]; i1++){
			i=fxy[ig][i1][0];
			j=fxy[ig][i1][1]-off;
			ii=i+sxy[ig][i1][0];
			jj=j+sxy[ig][i1][1]+2*off;
			for(k=j; k<=j+off; k++){
				xx[j][i]=sxy[ig][i1][1]/2.;
				zz[j][i]=sxy[ig][i1][0]/2.;
			}
			for(k=jj; k>jj-off; k--){
				xx[jj][ii]=-sxy[ig][i1][1]/2.;
				zz[jj][ii]=-sxy[ig][i1][0]/2.;
			}
		}
	} 
	lin_int2(xx, mm);
	lin_int2(zz, mm);

	/* calculate tangent vector */
	pp=sf_floatalloc2(n1, n2);
	sf_floatread(pp[0], n12, dip);
	u12=sf_floatalloc3(3, n1, n2);
	dx=sf_floatalloc2(n1, n2);
	dz=sf_floatalloc2(n1, n2);
	dxx=sf_floatalloc2(n1, n2);
	dzz=sf_floatalloc2(n1, n2);
	for(i2=0; i2<n2; i2++){
		for(i1=0; i1<n1; i1++){
			norm=sqrtf(1.+pp[i2][i1]*pp[i2][i1]);
			u1=pp[i2][i1]/norm;
			u2=1./norm;
			u12[i2][i1][0]=u1*u1+0.1*u2*u2;
			u12[i2][i1][1]=u1*u2-0.1*u1*u2;
			u12[i2][i1][2]=u2*u2+0.1*u1*u1;
		}
	}

	/* structure oriented smoothing */
	nknown=n1;
	for(ig=0; ig<ng; ig++) nknown+=np[ig];
	y=sf_floatalloc(nknown);
	memset(y, nknown, sizeof(float));
	k=n1;
	for (ig=0; ig<ng; ig++){
		for (i1=0; i1<np[ig]; i1++){
			y[k]=sxy[ig][i1][1];
			k++;
		}
	}

	sf_conjgrad_init(n12, n12, nknown, nknown, lam, 10*FLT_EPSILON, true, true);
	for(i2=0; i2<n2; i2++) for(i1=0; i1<n1; i1++) p[i2*n1+i1]=xx[i2][i1];
	sf_conjgrad(NULL, sf_fault_lop, sf_shape_lop, p, xx[0], y, niter);
	//sf_solver_reg(sf_fault_lop, sf_cgstep, sf_grad_lop, n12, n12, nknown, p, y, niter, lam, "x0", xx[0], "end");
	
	sf_putint(out, "n3", 2);
	sf_floatwrite(xx[0], n12, out);

	k=n1;
	for (ig=0; ig<ng; ig++){
		for (i1=0; i1<np[ig]; i1++){
			y[k]=sxy[ig][i1][0];
			k++;
		}
	}
	
	//sf_solver_reg(sf_fault_lop, sf_cgstep, sf_grad_lop, n12, n12, nknown, p, y, niter, lam, "x0", zz[0], "end");
	for(i2=0; i2<n2; i2++) for(i1=0; i1<n1; i1++) p[i2*n1+i1]=zz[i2][i1];
	sf_conjgrad(NULL, sf_fault_lop, sf_shape_lop, p, zz[0], y, niter);
	sf_conjgrad_close();

	sf_floatwrite(zz[0], n12, out);

	exit(0);
}
