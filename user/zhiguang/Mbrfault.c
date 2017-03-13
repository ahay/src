/* Bridge fault zones with smooth transition */
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

int main(int argc, char* argv[])
{
	bool replace, mode;
	int i1, i2, n1, n2, ig, ng, i, j;
	int *np, width, ***sxy, order;
	float **dd, **ss, **mm, ***dg, **bb, **ff, *trace1, *trace2;
	sf_file in, out, slip, mask=NULL, shift=NULL;

	sf_init(argc, argv);

	if(!sf_getbool("replace", &replace)) replace=false;
	if(!sf_getbool("mode", &mode)) mode=true;

	in=sf_input("in");
	slip=sf_input("slip");
	out=sf_output("out");
	if(replace) shift=sf_input("shift");
	else mask=sf_output("mask");

	if(!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if(!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if(!sf_getint("width", &width)) width=9;
	if(!sf_getint("order", &order)) order=1;

	dd=sf_floatalloc2(n1, n2);
	ss=sf_floatalloc2(n1, n2);
	mm=sf_floatalloc2(n1, n2);

	sf_floatread(dd[0], n1*n2, in);
	sf_floatread(ss[0], n1*n2, slip);
	memset(mm[0], 0, n1*n2*sizeof(float));

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

	/* extract fault boundary, position, slip value, number of sample */
	np=sf_intalloc(ng);
	dg=sf_floatalloc3(n1, width, ng);
	bb=sf_floatalloc2(n1, 2*ng);
	sxy=sf_intalloc3(2, n1, ng);
	memset(np, 0, ng*sizeof(int));
	memset(dg[0][0], 0, ng*width*n1*sizeof(float));
	//memset(bb[0], 0, 2*ng*n1*sizeof(float));

	for (ig=0; ig<ng; ig++){
		for (i1=0; i1<n1; i1++){
			for (i2=0; i2<n2; i2++){
				if(mm[i2][i1]==ig+1){
					bb[ig][i1]=dg[ig][0][i1]=dd[i2-width/2][i1];
					bb[ig+ng][i1]=dg[ig][width-1][i1]=dd[i2+width/2][i1];
					sxy[ig][np[ig]][0]=i1;
					sxy[ig][np[ig]][1]=i2;
					np[ig]++;
				} // one fault
			} // end of i2
		} // end of i1

		for (i1=0; i1<sxy[ig][0][0]; i1++){
			j=sxy[ig][0][1];
			bb[ig][i1]=bb[ig+ng][i1]=dd[j][i1];
		}
		if (sxy[ig][np[ig]-1][0]<n1-1){
			for (i1=sxy[ig][np[ig]-1][0]+1; i1<n1; i1++){
				j=sxy[ig][np[ig]-1][1];
				bb[ig][i1]=bb[ig+ng][i1]=dd[j][i1];
			}
		}
	} // end of ig

	if(!replace){
		sf_floatwrite(mm[0], n1*n2, mask);
		sf_putint(out, "n2", 2*ng);
		sf_floatwrite(bb[0], n1*2*ng, out);
		exit(0);
	}

	/* calculate dip */
	ff=sf_floatalloc2(n1, 2*ng);
	sf_floatread(ff[0], n1*2*ng, shift);
	for (ig=0; ig<2*ng; ig++){
		for (i1=0; i1<n1; i1++){
			ff[ig][i1]/=(width-1);
		}
	}

	/* predictive painting */
	trace1=sf_floatalloc(n1);
	trace2=sf_floatalloc(n1);
	predict_init(n1, width, 0.0001, order, 1, false);
	for (ig=0; ig<ng; ig++){
		for (i1=0; i1<n1; i1++){
			dg[ig][0][i1]/=2.;
			dg[ig][width-1][i1]/=2.;
			trace1[i1]=dg[ig][0][i1];
			trace2[i1]=dg[ig][width-1][i1];
		}
		for (i2=1; i2<width; i2++){
			predict_step(false, true, trace1, ff[ig]); // right direction
			for (i1=0; i1<n1; i1++) dg[ig][i2][i1]+=trace1[i1];
			predict_step(false, true, ff[ig], ff[ig]);

			predict_step(false, false, trace2, ff[ig+ng]); // left direction
			for (i1=0; i1<n1; i1++) dg[ig][width-i2-1][i1]+=trace2[i1];
			predict_step(false, false, ff[ig+ng], ff[ig+ng]);
		}
	}

	/* replace fault */
	for (ig=0; ig<ng; ig++){
		for (i2=0; i2<width; i2++){
			for (i1=0; i1<np[ig]; i1++){
				i=sxy[ig][i1][0];
				j=sxy[ig][i1][1];
				dd[j+i2-width/2][i]=dg[ig][i2][i];
			}
		}
	}

	sf_floatwrite(dd[0], n1*n2, out);

	exit(0);
}
