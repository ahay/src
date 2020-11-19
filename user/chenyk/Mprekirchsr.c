/* 2-D Prestack Kirchhoff time migration with antialiasing in source-receiver domain. 
The axes in the input are {time,receiver,source}
The axes in the output are {time,receiver}
I(t0,x0)=\int (W(t0,y0,s,r) * u(s,r,T(s,t0,x0)+T(r,t0,x0)) dsdr
*/
/*
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

void aal_doubint(int nt, float *trace)
/* causal and anticausal integration */
{
    int it;
    float tt;

    tt = trace[0];
    for (it=1; it < nt; it++) {
		tt += trace[it];
		trace[it] = tt;
    }
    tt = trace[nt-1];
    for (it=nt-2; it >=0; it--) {
		tt += trace[it];
		trace[it] = tt;
    }
}
    
float aal_pick(float ti, float deltat, 
	       const float *trace,
	       int nt, float dt, float t0)
/* pick a traveltime sample from a trace */
{
    int it, itm, itp;
    float ft, tm, tp, ftm, ftp, imp;

    ft = (ti-t0)/dt; it = floorf(ft); ft -= it; 
    if ( it < 0 || it >= nt-1) return 0.0;
 
    tm = ti-deltat-dt;
    ftm = (tm-t0)/dt; itm = floorf(ftm); ftm -= itm; 
    if (itm < 0) return 0.0;
                 
    tp = ti+deltat+dt;
    ftp = (tp-t0)/dt; itp = floorf(ftp); ftp -= itp; 
    if (itp >= nt-1) return 0.0;

    imp = dt/(dt+tp-tm);
    imp *= imp;
		
    return imp*(
	2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
	(1.-ftm)*trace[itm] - ftm*trace[itm+1]    - 
	(1.-ftp)*trace[itp] - ftp*trace[itp+1]);    
}

int main(int argc, char* argv[])
{
    int nt, nx, ns, nz, ix, is, iz, ir, i, it;
    float *trace, **mod=NULL, **dat=NULL, **v;
    float z, dx, ti, x0, t0, dt, z0, dz, vi, aal, ft, f,g;
	float s, s0, ds, t1, t2;
	bool adj;
    sf_file inp,outp, vel;

    sf_init (argc,argv);
    vel = sf_input("vel");
	inp = sf_input("in");
    outp = sf_output("out");
	

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2=");
    if (!sf_histint(inp,"n3",&ns)) sf_error("No n3=");

    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"o2",&x0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1=");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(inp,"o3",&s0)) sf_error("No o3=");
    if (!sf_histfloat(inp,"d3",&ds)) sf_error("No d3=");


    if (!sf_getint("nz",&nz)) 		nz=nt;
    if (!sf_getfloat("dz",&dz)) 	dz=dt;
    if (!sf_getfloat("z0",&z0)) 	z0=t0;
	if (!sf_getbool("adj",&adj))   	adj=false;

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */

    v = sf_floatalloc2(nz,nx);
    sf_floatread(v[0],nz*nx,vel);

    trace = sf_floatalloc(nt);

	if(adj)
	{
    	mod = sf_floatalloc2(nz,nx);
	}else{
    	dat = sf_floatalloc2(nz,nx);
	}



	for(is=0;is<ns;is++){
		sf_warning("Shot %d of %d", is+1, ns);

		s=s0+is*ds;
		

    	for (i=0; i < nz*nx; i++) {
			if(adj)
			{
				mod[0][i] = 0.;
			}else{
				dat[0][i] = 0.;
			}			
    	}	

    	/* loop over input traces (receiver position) */ 
    	for (ir=0; ir < nx; ir++) {

				sf_floatread (trace,nt,inp);

			// aal_doubint(nt,trace);
        
			/* loop over output traces (receiver position)*/
			for (ix=0; ix < nx; ix++) { 

	    		/* loop over output time */
	    		for (iz=0; iz < nz; iz++) {
					z = z0 + iz*dz;  
					vi = v[ix][iz];

					/* hypot(a,b) = sqrt(a*a+b*b) */
					t1 = hypotf(0.5*z,((ix-ir)*dx)/vi);		
					t2 = hypotf(0.5*z,(x0+ix*dx-s)/vi);		
					ti=t1+t2;
					
					ft=(ti-t0)/dt;
					it=floorf(ft);
					
					if(adj)
					{
						f=it+1-ft;
						g=ft-it;
						if (f<0) f=0; 
						if(f>1) f=1; 
						if (g<0) g=0; 
						if (g>1) g=1;
						//sf_warning("f=%g,g=%g",f,g);
						mod[ix][iz] += trace[it]*f + trace[it+1]*g;
					}
					else
					{
						f=it+1-ft;
						g=ft-it;
						if (f<0) f=0; 
						if(f>1) f=1; 
						if (g<0) g=0; 
						if (g>1) g=1;
						dat[ix][it] += trace[iz]*f;		
						dat[ix][it+1] += trace[iz]*g;
					}


	    		} 
			} 
    	}   
		if(adj)
		{
	    	sf_floatwrite(mod[0],nz*nx,outp); 
		}else{
	    	sf_floatwrite(dat[0],nz*nx,outp); 
		}

	}/* is */       

    exit(0);
}


