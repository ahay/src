/* Fast marching eikonal solver (3-D) in spherical coordinates. 

Also see sfeikonal (cartesian version)
http://ahay.org/blog/2014/06/11/program-of-the-month-sfeikonal/
*/
/*
  Copyright (C) 2022 University of Texas at Austin
  
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

#include<rsf.h>

#include "fastmarchrtp.h"

int main (int argc,char* argv[]) 
{
    int b1, b2, b3, n1, n2, n3, i, nshot, ndim, is,order,n123, *p;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3, slow;
    float **s, *t, *v;
    char *sfile;
    bool isvel, sweep, plane[3];
    sf_file vel, time, shots;

    sf_init (argc, argv);
    vel = sf_input("in");
    time = sf_output("out");

    if (SF_FLOAT != sf_gettype(vel)) 
	sf_error("Need float input");
    if(!sf_histint(vel,"n1",&n1)) sf_error("No n1= in input");
    if(!sf_histint(vel,"n2",&n2)) sf_error("No n2= in input");
    if(!sf_histint(vel,"n3",&n3)) n3=1;

    if(!sf_histfloat(vel,"d1",&d1)) sf_error("No d1= in input");
    if(!sf_histfloat(vel,"d2",&d2)) sf_error("No d2= in input");
    if(!sf_histfloat(vel,"d3",&d3)) d3=d2;

    if(!sf_histfloat(vel,"o1",&o1)) o1=0.;
    if(!sf_histfloat(vel,"o2",&o2)) o2=0.;
    if(!sf_histfloat(vel,"o3",&o3)) o3=0.;

    if(!sf_getbool("vel",&isvel)) isvel=true;
    /* if y, the input is velocity; n, slowness squared */

    if(!sf_getint("order",&order)) order=2;
    /* [1,2] Accuracy order */

    if (!sf_getbool("sweep",&sweep)) sweep=false;
    /* if y, use fast sweeping instead of fast marching */

    if(!sf_getfloat("br1",&br1)) br1=d1;    
    if(!sf_getfloat("br2",&br2)) br2=d2; 
    if(!sf_getfloat("br3",&br3)) br3=d3;
    /* Constant-velocity box around the source (in physical dimensions) */
 
    if(!sf_getbool("plane1",&plane[2])) plane[2]=false;
    if(!sf_getbool("plane2",&plane[1])) plane[1]=false;
    if(!sf_getbool("plane3",&plane[0])) plane[0]=false;
    /* plane-wave source */

    if(!sf_getint("b1",&b1)) b1= plane[2]? n1: (int) (br1/d1+0.5); 
    if(!sf_getint("b2",&b2)) b2= plane[1]? n2: (int) (br2/d2+0.5); 
    if(!sf_getint("b3",&b3)) b3= plane[0]? n3: (int) (br3/d3+0.5); 
    /* Constant-velocity box around the source (in samples) */

    if( b1<1 ) b1=1;  
    if( b2<1 ) b2=1;  
    if( b3<1 ) b3=1;

    sfile = sf_getstring("shotfile");
    /* File with shot locations (n2=number of shots, n1=3) */

    if(NULL != sfile) {
	shots = sf_input("shotfile");

	if (SF_FLOAT != sf_gettype(shots)) 
	    sf_error("Need float shotfile");
	if(!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
	if(!sf_histint(shots,"n1",&ndim) || ndim != 3) 
	    sf_error("Need n1=3 in shotfile");
  
	s = sf_floatalloc2 (ndim,nshot);
	sf_floatread(s[0],nshot*ndim,shots);
	sf_fileclose(shots);
    
	sf_putint (time,"n4",nshot);
	free (sfile);
    } else {
	nshot = 1;
	ndim = 3;
    
	s = sf_floatalloc2 (ndim,nshot);     

	if(!sf_getfloat("zshot",&s[0][0])  ) s[0][0]=0.; 
	/* Shot location (used if no shotfile) */
	if(!sf_getfloat("yshot",&s[0][1])) s[0][1]=o2 + 0.5*(n2-1)*d2;
	if(!sf_getfloat("xshot",&s[0][2])) s[0][2]=o3 + 0.5*(n3-1)*d3;

	sf_warning("Shooting from zshot=%g yshot=%g xshot=%g",
		   s[0][0],s[0][1],s[0][2]);
    }

    n123 = n1*n2*n3;

    t  = sf_floatalloc (n123);
    v  = sf_floatalloc (n123);
    p  = sf_intalloc   (n123);

    sf_floatread(v,n123,vel);
    
	/*extra variables for spherical case*/
    int nr,nt,np;
    float dr,dt,dp;
    float or,ot,op;
	float tmp, rabs, rat, delr; 
	float tabs, theta, delt; 
	float pabs, ph, delp;
    float *r, *the, *phi, *z, *x, *y;
    float *trtp, *vrtp;
    float xs[3];
    int ir,it,ip, i1, i2, i3;
    
    nr=(n1-1)*1.8+1;/*2.0 can be polished*/
    nt=n2;
    if(n3==1) np=n3;
    else np=n3;
    dr=d1;
    dt=SF_PI*2/(nt-1);
    if(n3==1) 
    {dp=1;op=0;}
    else 
    {dp=SF_PI*2/(np-1);op=-SF_PI;}
    or=0;
    ot=-SF_PI;
    
    b1=1;
    b2=nt;
	if(n3==1)
		b3=1;
	else
		b3=np;
	
    sf_file time2;
    time2 = sf_output("time2");
    sf_putint(time2,"n1",nr);
    sf_putint(time2,"n2",nt);
    sf_putint(time2,"n3",np);
    sf_putfloat(time2,"d1",dr);
    sf_putfloat(time2,"d2",dt);
    sf_putfloat(time2,"d3",dp);
    sf_putfloat(time2,"o1",or);
    sf_putfloat(time2,"o2",ot);
    sf_putfloat(time2,"o3",op);
    
	/*extra arrays for spherical case*/
    trtp  = sf_floatalloc (nr*nt*np);
    vrtp  = sf_floatalloc (nr*nt*np);
    
    r=sf_floatalloc(nr);
    the=sf_floatalloc(nt);
    phi=sf_floatalloc(np);
    z=sf_floatalloc(n1);
    x=sf_floatalloc(n2);
    y=sf_floatalloc(n3);
    for(ir=0;ir<nr;ir++) r[ir]=or+ir*dr;
    for(it=0;it<nt;it++) the[it]=ot+it*dt;
    for(ip=0;ip<np;ip++) phi[ip]=op+ip*dp;
    for(i1=0;i1<n1;i1++) z[i1]=o1+i1*d1; 
    for(i2=0;i2<n2;i2++) x[i2]=o2+i2*d2; 
    for(i3=0;i3<n3;i3++) y[i3]=o3+i3*d3;  
    
    if (isvel) {
	/* transform velocity to slowness squared */
	for(i = 0; i < n123; i++) {
	    slow = v[i];
	    v[i] = 1./(slow*slow);
	}
    } 
    
    if (!sweep) fastmarchrtp_init (np,nt,nr);
 
    /* loop over shots */
    for( is = 0; is < nshot; is++) {
    xs[0]=s[is][0];xs[1]=s[is][1];xs[2]=s[is][2];
    
    sf_warning("xs=[%g,%g,%g]",xs[0],xs[1],xs[2]);
    sf_warning("nr=%d,nt=%d,np=%d",nr,nt,np);
    sf_warning("n1=%d,n2=%d,n3=%d",n1,n2,n3);
    sf_warning("dr=%g,dt=%g,dp=%g",dr,dt,dp);
    sf_warning("d1=%g,d2=%g,d3=%g",d1,d2,d3);
    sf_warning("or=%g,ot=%g,op=%g",or,ot,op);
    sf_warning("o1=%g,o2=%g,o3=%g",o1,o2,o3);
    sf_warning("b1=%d,b2=%d,b3=%d",b1,b2,b3);
    sf_warning("xs0=%g,xs1=%g,xs2=%g",xs[0],xs[1],xs[2]);
    sf_warning("pi=%g,sin(pi)=%g,cos(pi)=%g",SF_PI,sinf(SF_PI),cosf(SF_PI));  
    
	/*velocity transformation from Cartesian to Spherical*/
    for(ir=0;ir<nr;ir++)
    	for(it=0;it<nt;it++)
    		for(ip=0;ip<np;ip++)
    		{
    			i1=floorf(r[ir]*cosf(the[it])+s[is][0])/d1;
    			i2=floorf(r[ir]*sinf(the[it])*cosf(phi[ip])+s[is][1])/d2;
    			i3=floorf(r[ir]*sinf(the[it])*cosf(phi[ip])+s[is][2])/d3;
    			
    			if(i1>n1-1) i1=n1-1;
    			if(i1<0) 	i1=0;
    			if(i2>n2-1) i2=n2-1;
    			if(i2<0)	i2=0;
    			if(i3>n3-1)	i3=n3-1;
    			if(i3<0)	i3=0;
    			vrtp[ir+it*nr+ip*nr*nt]=v[i1+i2*n1+i3*n1*n2];
// 				vrtp[ir+it*nr+ip*nr*nt]=1/1.5/1.5;
    		}
//     for(ir=0;ir<nr*nt*np;ir++)
//     if(vrtp[ir]!=1/1.5/1.5) sf_warning("vrtp is %g",vrtp[ir]);
    
	sf_warning("shot %d of %d;",is+1,nshot);
	if (sweep) {
	    continue;
	} else {
	    fastmarchrtp(trtp,vrtp,p, plane,
		      np,nt,nr,
		      op,ot,or,
		      dp,dt,dr,
		      0,0,0, 
		      b3,b2,b1,
		      order);
	}	

	/*traveltime transformation from Spherical to Cartesian*/
	
	for(i1=0;i1<n1;i1++)
		for(i2=0;i2<n2;i2++)
			for(i3=0;i3<n3;i3++)
				{
				tmp=sqrtf((z[i1]-xs[0])*(z[i1]-xs[0])+(x[i2]-xs[1])*(x[i2]-xs[1])+(y[i3]-xs[2])*(y[i3]-xs[2]) );
        		rabs=tmp;   /*absolute r*/
        		rat=(tmp-r[0])/dr;
				ir=floorf(rat);
				delr=rat-ir*dr;
// 				sf_warning("r is finished = %d", ir);
        		tmp=atanf(sqrtf((x[i2]-xs[1])*(x[i2]-xs[1])+(y[i3]-xs[2])*(y[i3]-xs[2]))/(z[i1]-xs[0]+0.0000000000001));
				if(z[i1]-xs[0]<0)
                	tmp=-SF_PI-tmp;
        		tabs=tmp;   	/*absolute t*/
        		theta=(tmp-the[0])/dt;
        		it=floorf(theta);
        		delt=theta-it*dt;
// 				sf_warning("t is finished = %d", it);
				
				
				tmp=atanf((y[i3]-xs[2])/(x[i2]-xs[1]+0.0000000000001));
				if(x[i2]-xs[1] <0)
					if (y[i2]-xs[2]>0)
						tmp=-SF_PI-tmp;
        		pabs=tmp;   /*absolute p*/
        		ph=(tmp-phi[0])/dt;
        		ip=floorf(ph);
        		delp=ph-ip*dp;
// 				sf_warning("p is finished = %d", ip);
				
				t[i1+i2*n1+i3*n1*n2] = trtp[ir+it*nr+ip*nr*nt];
				}
				
	sf_floatwrite (trtp,nr*nt*np,time2);
	sf_floatwrite (t,n123,time);
    }
    sf_warning(".");
	
	fastmarchrtp_close();
	
    exit (0);
}

