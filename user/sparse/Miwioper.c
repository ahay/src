/* Image-domain waveform tomography (linear operator). */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "iwioper.h"

int main(int argc, char* argv[])
{
    bool adj, load, mass;
    int n1, n2, npml; 
    int nh, ns, nw;
    float d1, d2, **vel, dw, ow;
    float *di, *dm, ***wght, **prec;
    char *datapath;
    sf_file in, out, model, us, ur;
    sf_file weight, precon;
    int uts, mts;
    char *order;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");    

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (!sf_getint("nh",&nh)) nh=0;
    /* horizontal space-lag */

    if (!sf_getbool("load",&load)) load=false;
    /* load LU */

    if (!sf_getint("uts",&uts)) uts=0;
    /* number of OMP threads */

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    uts = (uts < 1)? mts: uts;

    if (!sf_getint("npml",&npml)) npml=10;
    /* PML width */

    if (NULL == (order = sf_getstring("order"))) order="j";
    /* discretization scheme (default optimal 9-point) */

    if (!sf_getbool("mass",&mass)) mass=false;
    /* if y, use discretization-based mass term */

    /* read model */
    if (NULL == sf_getstring("model"))
	sf_error("Need model=");
    model = sf_input("model");

    if (!sf_histint(model,"n1",&n1)) sf_error("No n1= in model.");
    if (!sf_histint(model,"n2",&n2)) sf_error("No n2= in model.");

    if (!sf_histfloat(model,"d1",&d1)) sf_error("No d1= in model.");
    if (!sf_histfloat(model,"d2",&d2)) sf_error("No d2= in model.");

    vel = sf_floatalloc2(n1,n2);
    sf_floatread(vel[0],n1*n2,model);
    
    if (load)
	datapath = sf_histstring(model,"in");	
    else
	datapath = NULL;

    /* allocate memory */
    dm = sf_floatalloc(n1*n2);
    di = sf_floatalloc(n1*n2*(2*nh+1));

    /* read input */
    if (adj)
	sf_floatread(di,n1*n2*(2*nh+1),in);
    else
	sf_floatread(dm,n1*n2,in);

    /* read source wavefield */
    if (NULL == sf_getstring("us"))
	sf_error("Need source wavefield us=");
    us = sf_input("us");

    if (!sf_histint(us,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histint(us,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(us,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(us,"o4",&ow)) sf_error("No ow=.");    
    
    /* read receiver wavefield */	
    if (NULL == sf_getstring("ur"))
	sf_error("Need receiver wavefield ur=");
    ur = sf_input("ur");

    /* read image weight */
    if (NULL == sf_getstring("weight")) {
	weight = NULL;
	wght = NULL;
    } else {
	weight = sf_input("weight");
	wght = sf_floatalloc3(n1,n2,2*nh+1);
	sf_floatread(wght[0][0],n1*n2*(2*nh+1),weight);
	sf_fileclose(weight);
    }

    /* read model preconditioner */
    if (NULL == sf_getstring("precon")) {
	precon = NULL;
	prec = NULL;
    } else {
	precon = sf_input("precon");
	prec = sf_floatalloc2(n1,n2);
	sf_floatread(prec[0],n1*n2,precon);
	sf_fileclose(precon);
    }

    /* write output header */
    if (adj)
	sf_putint(out,"n3",1);
    else
	sf_putint(out,"n3",2*nh+1);    
    
    /* initialize */
    iwi_init(npml, n1,n2,d1,d2, nh,ns,ow,dw,nw,
	     us,ur, load,datapath, uts, order,mass);

    /* set velocity and weight */
    iwi_set(vel,wght,prec);

    /* linear operator */
    iwi_oper(adj,false,n1*n2,n1*n2*(2*nh+1),dm,di);

    /* write output */
    if (adj)
	sf_floatwrite(dm,n1*n2,out);
    else
	sf_floatwrite(di,n1*n2*(2*nh+1),out);

    exit(0);
}
