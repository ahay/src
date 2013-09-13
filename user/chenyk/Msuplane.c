/* Create common offset data file with up to 3 planes 
Command line should be: sfsuplane in=NULL >file.rsf 
(There should be at least one parameter in command line because of sf_init(argc,argv))
*/
/*
  Copyright (C) 2007 Colorado School of Mines
  Copyright (C) 2013 University of Texas at Austin
   
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

#define NT	64
#define NTR	32
#define DT	0.004
#define OFF	400
#define NPL	3

int main(int argc, char *argv[])
{
	float dip1;		/* time-dip of plane 1 (ms/trace)	*/
	float dip2;		/* time-dip of plane 2 (ms/trace)	*/
	float dip3;		/* time-dip of plane 3 (ms/trace)	*/
	float dt;		/* time sample interval in seconds 	*/
	float eps;		/* fit - itless (for linear interpol)	*/
	float fit;		/* fractional time sample		*/
	float offset;		/* constant offset		*/
	float *data=NULL;		/* outdata */
	int ct1,cx1;		/* center of plane 1 (sample and trace)	*/
	int ct2,cx2;		/* center of plane 2 (sample and trace)	*/
	int ct3,cx3;		/* center of plane 3 (sample and trace)	*/
	int itless;		/* sample above plane intersection	*/
	int itmore;		/* sample below plane intersection 	*/
	int itr;		/* trace counter			*/
	int len1;		/* HORIZ extent of plane 1		*/
	int len2;		/* HORIZ extent of plane 2		*/
	int len3;		/* HORIZ extent of plane 3		*/
	int liner;		/* flag for special output section	*/
	float msps;		/* milliseconds per sample 		*/
	int npl;		/* number of planes 			*/
	int nt;			/* time samples in outdata		*/
	int ntr;		/* traces in outdata			*/
	int tfe1;		/* traces-from-end of plane 1 (for taper) */
	int tfe2;		/* traces-from-end of plane 2 (for taper) */
	int tfe3;		/* traces-from-end of plane 3 (for taper) */
	int taper;		/* flag to taper plane ends to zero	*/

	sf_init(argc,argv);

	sf_file out;

	/* set parameters  */
	if(!sf_getint("nt",&nt)) 			nt = NT;
	if(!sf_getint("ntr", &ntr)) 		ntr = NTR;
	if(!sf_getfloat("dt", &dt)) 		dt=DT;	
	if(!sf_getfloat("offset", &offset)) offset=OFF;	
	if(!sf_getint("npl", &npl)) 		npl = NPL;
	if(!sf_getint("taper", &taper)) 	taper = 0;

	/* set defaults and/or get parameters for plane 1 */
	if(!sf_getfloat("dip1", &dip1)) dip1=0;
	if(!sf_getint("len1", &len1))   len1 = 3*ntr/4;	
	if(!sf_getint("ct1", &ct1))		ct1 = nt/2-1;	
	if(!sf_getint("cx1", &cx1))		cx1 = ntr/2-1;

	/* set defaults and/or get parameters for plane 2 */
	if(!sf_getfloat("dip2", &dip2)) dip2 = 4;	
	if(!sf_getint("len2", &len2)) 	len2 = 3*ntr/4;	
	if(!sf_getint("ct2", &ct2))		ct2 = nt/2-1;
	if(!sf_getint("cx2", &cx2))		cx2 = ntr/2-1;		

	/* set defaults and/or get parameters for plane 3 */
	if(!sf_getfloat("dip3", &dip3)) dip3 = 8;	
	if(!sf_getint("len3", &len3)) 	len3 = 3*ntr/4;	
	if(!sf_getint("ct3", &ct3)) 	ct3 = nt/2-1;		
	if(!sf_getint("cx3", &cx3))		cx3 = ntr/2-1;		

	/* check if user wants the special output specified */
        /* by liner=1; if so, set parameters accordingly    */
	if(!sf_getint("liner", &liner))	liner = 0;	
	if (liner == 1) {
		nt = 64;
		ntr = 64;
		npl = 3;	

		dip1 = 0;
		len1 = ntr/4;	
		ct1 = nt/2;	ct1 -= 1;		
		cx1 = 3*ntr/4;	cx1 -= 1;

		dip2 = 4;
		len2 = ntr/4;	
		ct2 = nt/2;	ct2 -= 1;
		cx2 = ntr/2;	cx2 -= 1;	

		dip3 = 8;
		len3 = ntr/4;	
		ct3 = nt/2;	ct3 -= 1;
		cx3 = ntr/4;	cx3 -= 1;	
	}

	/* calculate milliseconds per sample */
	msps = dt*1000.0;	

	tfe1 = 0; tfe2 = 0; tfe3 = 0;

	/* allocate memory for output data */
	data=sf_floatalloc(nt);

    out = sf_output("out");
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",0);
	sf_putint(out,"n2",ntr);
	sf_putfloat(out,"d2",0.01);
	sf_putfloat(out,"o2",0);

	for (itr = 0; itr < ntr; itr++) {
		memset( (void *) data, 0, nt * sizeof(float));

		/* plane 1 */
		if (itr >= cx1-len1/2 && itr <= cx1+len1/2) {
		    ++tfe1;

		    /* fit is fractional sample of plane intersection */
		    fit = ct1 - ( cx1 - itr ) * dip1 / msps; 
		    if (fit >= 0 && fit <= (float) nt) {

			/* linear interpolation */
			itless = fit;
			eps = fit - itless;
			itmore = fit + 1;
			data[itless] += 1.0 - eps;	 
			data[itmore] += eps;	 

			/* taper option */
			if (taper == 1) {
			  /* last point */
			  if (tfe1 == 1 || tfe1 == len1 + 1) {
				data[itless] /= 6.0;	 
				data[itmore] /= 6.0;	 
			  } 
			  /*  next-to-last point */
			  if (tfe1 == 2 || tfe1 == len1) {
				data[itless] /= 3.0;	 
				data[itmore] /= 3.0;	 
			  } 
		    }
		  }
		}

		/*  plane 2  */
		if (npl > 1) {
		  if (itr >= cx2-len2/2 && itr <= cx2+len2/2) {
		    ++tfe2;

		    /* fit is fractional sample of plane intersection */
		    fit = ct2 - ( cx2 - itr ) * dip2 / msps; 
		    if (fit >= 0 && fit <= (float) nt) {

			/* linear interpolation */
			itless = fit;
			eps = fit - itless;
			itmore = fit + 1;
			data[itless] += 1.0 - eps;	 
			data[itmore] += eps;	 

			/* taper option */
			if (taper == 1) {
			  /* last point */
			  if (tfe2 == 1 || tfe2 == len2 + 1) {
				data[itless] /= 6.0;	 
				data[itmore] /= 6.0;	 
			  } 
			  /*  next-to-last point */
			  if (tfe2 == 2 || tfe2 == len2) {
				data[itless] /= 3.0;	 
				data[itmore] /= 3.0;	 
			  } 
		        }
		    }
		  }
		}

		/* plane 3  */
		if (npl > 2) {
		  if (itr >= cx3-len3/2 && itr <= cx3+len3/2) {
		    ++tfe3;

		    /* fit is fractional sample of plane intersection */
		    fit = ct3 - ( cx3 - itr ) * dip3 / msps; 
		    if (fit >= 0 && fit <= (float) nt) {

			/* linear interpolation */
			itless = fit;
			eps = fit - itless;
			itmore = fit + 1;
			data[itless] += 1.0 - eps;	 
			data[itmore] += eps;	 

			/* taper option */
			if (taper == 1) {
			  /* last point */
			  if (tfe3 == 1 || tfe3 == len3 + 1) {
				data[itless] /= 6.0;	 
				data[itmore] /= 6.0;	 
			  } 
			  /* next-to-last point */
			  if (tfe3 == 2 || tfe3 == len3) {
				data[itless] /= 3.0;	 
				data[itmore] /= 3.0;	 
			  } 
		        }
		    }
		  }
		}

	sf_floatwrite(data,nt,out);
	}
	/* free allocated memory */	
	free(data);

	exit(0);
}
