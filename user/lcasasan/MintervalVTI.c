/* Interval/Effective VTI parameters from Effective/Interval profiles */
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

#include <math.h>
#include <rsf.h>

int main (int argc, char* argv[])
{
    bool interval;
    int it, nt;
    float dt, t0, vnmo4;
    float v4th[2];
    float *velN=NULL, *e=NULL, *velN_out=NULL, *velH_out=NULL, *e_out=NULL ,*t=NULL, *v4th_out=NULL;
    sf_file vn=NULL, eta=NULL, vn_out=NULL, vh_out=NULL, eta_out=NULL;

    sf_init (argc,argv);
    vn = sf_input("in");
    
    vn_out = sf_output("out");

    
    if (NULL != sf_getstring("vH_out")) vh_out = sf_output("vH_out"); /*output HOR vel*/
    else vh_out = NULL ;
    
    if (NULL != sf_getstring("eta_out")) eta_out = sf_output("eta_out"); /*output eta*/
    else eta_out = NULL ;
	
    if (NULL != sf_getstring("eta")) eta = sf_input("eta"); /*input eta*/
    else eta = NULL ;

    if (!sf_getbool("interval",&interval)) interval=true;
    /* output are interval [y] or effective [n] profiles */

    if (SF_FLOAT != sf_gettype(vn)) sf_error("Need float input");
    if (!sf_histint(vn,"n1",&nt))   sf_error("No n1= in input");
    if (!sf_histfloat(vn,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(vn,"o1",&t0)) sf_error("No o1= in input");


    velN = sf_floatalloc(nt);
    velN_out = sf_floatalloc(nt);
    t = sf_floatalloc(nt);      

    sf_floatread(velN,nt,vn);
    velN_out[0]=velN[0]; /*initialization*/

    if (eta!=NULL) {
        e = sf_floatalloc(nt);
        e_out=sf_floatalloc(nt);
        velH_out = sf_floatalloc(nt);
        v4th_out = sf_floatalloc(nt);

    	sf_floatread(e,nt,eta); /*read eta from file*/

    	e_out[0]=e[0]; /*initialization*/
        v4th_out[0] = (1+8*e_out[0])*(velN_out[0]*velN_out[0]*velN_out[0]*velN_out[0]);
	}
    for (it=1;it<nt;it++){

    	t[it]=t0+it*dt;

    	/* Vnmo */
    	if (interval)
    		velN_out[it] = sqrt(fabs( ( (velN[it]*velN[it]*t[it]) -(velN[it-1]*velN[it-1]*t[it-1]))/ (t[it]-t[it-1]) ) );
    	else
    		velN_out[it] = sqrt(fabs( (velN_out[it-1]*velN_out[it-1]*t[it-1] + velN[it]*velN[it]*(t[it]-t[it-1])  ) / t[it] ) );

    	/* eta and Vhor */
		if (eta!=NULL) {
	        vnmo4=velN_out[it]*velN_out[it]*velN_out[it]*velN_out[it];
	        v4th[0]=(1+8*e[it])*velN[it]*velN[it]*velN[it]*velN[it];

	        if (interval) {
		        v4th[1]=(1+8*e[it-1])*velN[it-1]*velN[it-1]*velN[it-1]*velN[it-1];
		        v4th_out[it] = (v4th[0]*t[it]-v4th[1]*t[it-1])/(t[it]-t[it-1]);
	        } else
		        v4th_out[it] = (v4th_out[it-1]*t[it-1] +  v4th[0] *(t[it]-t[it-1]) )/t[it];

	        e_out[it] = 0.125 * ( ( v4th_out[it] / vnmo4 ) - 1 );
			velH_out[it] = velN_out[it]*sqrt(fabs(1+2*e_out[it]) );
		}
    }
    
    sf_floatwrite (velN_out,nt,vn_out);
    if ((eta_out!=NULL) && (eta!=NULL))
		sf_floatwrite (e_out,nt,eta_out);
    if ((vh_out!=NULL) && (eta!=NULL))
    	sf_floatwrite (velH_out,nt,vh_out);

    exit (0);
}
