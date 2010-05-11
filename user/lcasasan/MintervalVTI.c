/* Interval VTI parameters from effective profiles */
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
    
    int it, nt;
    float dt, t0, vNint4;
    float *velN=NULL, *e=NULL, *velN_int=NULL, *velH_int=NULL, *e_int=NULL ,*t=NULL, *f=NULL;
    sf_file vn=NULL, eta=NULL, vn_int=NULL, vh_int=NULL, eta_int=NULL; 

    sf_init (argc,argv);
    vn = sf_input("in");
    
    vn_int = sf_output("out"); 

    
	if (NULL != sf_getstring("vH_int")) vh_int = sf_output("vH_int"); /*interval HOR vel*/
	else vh_int = NULL ;
    
	if (NULL != sf_getstring("eta_int")) eta_int = sf_output("eta_int"); /*interval eta*/
	else eta_int = NULL ; 
	
    if (NULL != sf_getstring("eta")) eta = sf_input("eta"); /*effective eta*/ 
	else eta = NULL ;


    if (SF_FLOAT != sf_gettype(vn)) sf_error("Need float input");
    if (!sf_histint(vn,"n1",&nt))   sf_error("No n1= in input");
    if (!sf_histfloat(vn,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(vn,"o1",&t0)) sf_error("No o1= in input");

    //sf_warning("\nnt %d",nt);

    velN = sf_floatalloc(nt);
    e = sf_floatalloc(nt);
    velN_int = sf_floatalloc(nt);
    velH_int = sf_floatalloc(nt);
    e_int=sf_floatalloc(nt);  
    t = sf_floatalloc(nt);      
    f = sf_floatalloc(nt);    

    sf_floatread(velN,nt,vn);
	/*velN */
    sf_floatread(e,nt,eta); 
    //sf_warning("velN %f eta %f ",velN[1],e[1]);
    
    for (it=0;it<nt;it++){
    t[it]=t0+it*dt;
    }
    
    velN_int[0]=velN[0];
    e_int[0]=e[0];
    velH_int[0]=velN[0]*sqrt(fabs(1+2*e[0]));
    for (it=1;it<nt;it++){
    velN_int[it]=sqrt(fabs( ( (velN[it]*velN[it]*t[it]) -(velN[it-1]*velN[it-1]*t[it-1]))/ (t[it]-t[it-1]) ) );   
    

    /*f(i) = Vnmo_eff(i)^2 * (4*Vh_eff(i)^2 - 3*Vnmo_eff(i)^2);*/
    f[it] = (velN[it]*velN[it])* ( 4 * (velN[it]*velN[it])*(1+2*e[it]) - 3 * (velN[it]*velN[it]) );

    vNint4 = velN_int[it]*velN_int[it]*velN_int[it]*velN_int[it];
    //velH_int[it] = velN_int[it] * sqrt ( (1/(4*vNint4)) * (f[it]*t[it]-f[it-1]*t[it-1])/(t[it]-t[it-1]) + 3/4  );    

    /* eta_int(i) = 1/(8*Vnmo_int(i)^4) * ( (f(i)*to(i)-f(i-1)*to(i-1))/(to(i)-to(i-1)) - Vnmo_int(i)^4  );*/
    e_int[it] = 1/(8*vNint4) * ( (f[it]*t[it]-f[it-1]*t[it-1])/(t[it]-t[it-1]) - vNint4  );
    velH_int[it] = velN_int[it]*sqrt(fabs(1+2*e_int[it]) );
    }
    
    sf_floatwrite (velN_int,nt,vn_int);
    sf_floatwrite (velH_int,nt,vh_int);
    sf_floatwrite (e_int,nt,eta_int);
    
    
    sf_close();
    exit (0);
}
