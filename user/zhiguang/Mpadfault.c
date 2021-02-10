/* Horizontally pad fault  */
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
	sf_map4 mo;
	bool extend, mode;
	int i1, i2, n1, n2, nn2, n22, ig, ng, i, j, flag, conf, order;
	int *np, ***tips, *width, **sxy, *x2, *nl, label[20], id[20];
	float d1, o1, **dd1, **ss, **mm, **dd2, **bb, **ff, *trace, **pp;
	sf_file in, out, slip, mask=NULL, bound=NULL, dip=NULL, shift=NULL, newdip=NULL, ppbig=NULL;

	sf_init(argc, argv);
	
	if(!sf_getbool("extend", &extend)) extend=true;
	if(!sf_getbool("mode", &mode)) mode=true;
	if(!sf_getint("conf", &conf)) conf=1;
	if(!sf_getint("order", &order)) order=2;

	in=sf_input("in");
	slip=sf_input("slip");
	out=sf_output("out");
	if(extend){ 
		mask=sf_output("mask");
		bound=sf_output("bound");
	}else{
		dip=sf_input("dip");
		shift=sf_input("shift");
		newdip=sf_output("newdip");
		ppbig=sf_output("ppbig");
	}

	if(!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if(!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if(!sf_histfloat(in, "d1", &d1)) sf_error("No d1= in input");
	if(!sf_histfloat(in, "o1", &o1)) sf_error("No o1= in input");

	dd1=sf_floatalloc2(n1, n2);
	ss=sf_floatalloc2(n1, n2);
	mm=sf_floatalloc2(n1, n2);

	sf_floatread(dd1[0], n1*n2, in);
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
						label[ng]=i2;
						ng++;
					}
				} // slip!=0 
			} // end if i1
		} // end of i2

		// sorting
		for(i=0; i<ng; i++) id[i]=i+1;
		for(i=0; i<ng; i++){
			for(j=i+1; j<ng; j++){
				if(label[j]<label[i]){
					nn2=label[j];
					label[j]=label[i];
					label[i]=nn2;
					n22=id[j];
					id[j]=id[i];
					id[i]=n22;
				}
			}
		} // end of ig
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

		// sorting
		for(i=0; i<ng; i++) id[i]=i+1;
	}

	np=sf_intalloc(ng);
	tips=sf_intalloc3(2,2,ng);
	width=sf_intalloc(ng);
	sxy=sf_intalloc2(n1, ng);

	/* find fault position */
	for (ig=0; ig<ng; ig++){
		flag=0;
		np[ig]=0;
		for (i1=0; i1<n1; i1++){
			for (i2=0; i2<n2; i2++){
				if(mm[i2][i1]==id[ig]){
					if(flag==0){ 
						tips[ig][0][0]=i1; 
						tips[ig][0][1]=i2; 
						flag=1;
					}
					sxy[ig][i1]=i2;
					np[ig]++;
				} // if (one fault)
			} // end of i2
		} // end of i1
	} // end of ig

	/* end point, fault width, and pad fault position */
	nn2=n2;
	for (ig=0; ig<ng; ig++){
		i=tips[ig][0][0]+np[ig]-1;
		tips[ig][1][0]=i;
		tips[ig][1][1]=sxy[ig][i];

		width[ig]=abs(tips[ig][1][1]-tips[ig][0][1])+30;
		nn2+=width[ig];

		for(i1=0; i1<tips[ig][0][0]; i1++) sxy[ig][i1]=tips[ig][0][1];
		for(i1=tips[ig][1][0]; i1<n1; i1++) sxy[ig][i1]=tips[ig][1][1];
	}

	/* re-organization */
	dd2=sf_floatalloc2(n1, nn2);
	nl=sf_intalloc(ng);
	x2=sf_intalloc(n1);
	memset(x2, 0, n1*sizeof(int));
	n22=0;
	bb=sf_floatalloc2(n1, 2*ng);

	for (ig=0; ig<ng; ig++){

		if(tips[ig][0][1]<tips[ig][1][1]) nl[ig]=n22+tips[ig][1][1]+15;
		else nl[ig]=n22+tips[ig][0][1]+15;

		for(i1=0; i1<n1; i1++){
			for (i2=x2[i1]; i2<n22+sxy[ig][i1]-conf; i2++) dd2[i2][i1]=dd1[i2-n22][i1];
			i=n22+sxy[ig][i1]-conf;
			j=n22+sxy[ig][i1]+width[ig]+conf;
			for (i2=i; i2<=nl[ig]; i2++) dd2[i2][i1]=dd1[sxy[ig][i1]-1-conf][i1];
			for (i2=nl[ig]+1; i2<=j; i2++) dd2[i2][i1]=dd1[sxy[ig][i1]+1+conf][i1];
		}
		n22+=width[ig];
		for(i1=0; i1<n1; i1++) x2[i1]=sxy[ig][i1]+n22+1+conf;

		// fault boundary
		for(i1=0; i1<n1; i1++){
			bb[ig][i1]=dd2[nl[ig]][i1];
			bb[ng+ig][i1]=dd2[nl[ig]+1][i1];
		}
	}

	for(i1=0; i1<n1; i1++){
		for (i2=x2[i1]; i2<nn2; i2++) dd2[i2][i1]=dd1[i2-n22][i1];
	}

	if(extend){ // extend the model to remove the fault
		sf_floatwrite(mm[0], n1*n2, mask);
		sf_putint(bound, "n2", 2*ng);
		sf_floatwrite(bb[0], n1*2*ng, bound);
		sf_putint(out, "n2", nn2);
		sf_floatwrite(dd2[0], n1*nn2, out);
		exit(0);
	}

	/* 'not extend' starts here ... */

	/* initialize dip of fault zones */
	sf_floatread(dd2[0], n1*nn2, dip);
	for (ig=0; ig<ng; ig++){
		for (i2=nl[ig]-10; i2<=nl[ig]+10; i2++)
			for (i1=0; i1<n1; i1++)
				dd2[i2][i1]=0.;
	}

	sf_putint(newdip, "n2", nn2);
	sf_floatwrite(dd2[0], n1*nn2, newdip);

	/* find vertical shifts */
	ff=sf_floatalloc2(n1, ng);
	sf_floatread(ff[0], n1*ng, shift);
	trace=sf_floatalloc(n1);

	for (i1=0; i1<n1; i1++) trace[i1]=i1*d1+o1;
	for (ig=0; ig<ng; ig++){
		for (i1=0; i1<n1; i1++) ff[ig][i1]+=trace[i1];
	}

	/* predictive painting */
	n22=0;
	pp=sf_floatalloc2(n1, nn2);
	predict_init(n1, nn2, 0.0001, order, 1, false);
	mo=sf_stretch4_init(n1, o1, d1, n1, 0.1);

	   // first column 
	for (i1=0; i1<n1; i1++) pp[0][i1]=trace[i1];

	   // central part
	for (ig=0; ig<ng; ig++){

		for(i2=n22+1; i2<nl[ig]; i2++){
			predict_step(false, true, trace, dd2[i2-1]);
			for (i1=0; i1<n1; i1++) pp[i2][i1]=trace[i1];
		}
		
		// inverse interpolation
		n22=nl[ig];
		sf_stretch4_define(mo, ff[ig], false);
		sf_stretch4_invert(false, mo, pp[n22], trace);
		
		//for (i1=n1-5; i1<n1; i1++) pp[n22][i1]=pp[n22][n1-6]; // stabilization

		for (i1=0; i1<n1; i1++) trace[i1]=pp[n22][i1];
	}

	    // last part
	for (i2=n22+1; i2<nn2; i2++){
		predict_step(false, true, trace, dd2[i2-1]);
		for (i1=0; i1<n1; i1++) pp[i2][i1]=trace[i1];
	}

	/* output expanded PP */
	sf_putint(ppbig, "n2", nn2);
	sf_floatwrite(pp[0], n1*nn2, ppbig);

	/* reduce model size */
	n22=0;
	memset(x2, 0, n1*sizeof(int));
	for (ig=0; ig<ng; ig++){
		for(i1=0; i1<n1; i1++)
			for (i2=x2[i1]; i2<n22+sxy[ig][i1]; i2++) 
				dd1[i2-n22][i1]=pp[i2][i1];
		n22+=width[ig];
		for(i1=0; i1<n1; i1++) x2[i1]=sxy[ig][i1]+n22;
	}

	for(i1=0; i1<n1; i1++){
		for (i2=x2[i1]; i2<nn2; i2++) dd1[i2-n22][i1]=pp[i2][i1];
	}

	sf_floatwrite(dd1[0], n1*n2, out);

	exit(0);
}
