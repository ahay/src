/* Reconstruct right-hand side from wavefield. */
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

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    int i, n1, j, n2; 
    int is, ns, iw, nw;
    float d1, d2, dw, ow;
    double omega;
    sf_complex ***wave, ***rhs;
    float **vel;
    double complex neib, cent;
    sf_file in, out, modl;
    int uts, mts;
    char *order;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("uts",&uts)) uts=0;
    /* number of OMP threads */

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    uts = (uts < 1)? mts: uts;

    if (NULL == (order = sf_getstring("order"))) order="j";
    /* discretization scheme (default optimal 9-point) */

    /* read model */
    if (NULL == sf_getstring("model"))
	sf_error("Need model=");
    modl = sf_input("model");

    if (!sf_histint(modl,"n1",&n1)) sf_error("No n1=.");
    if (!sf_histint(modl,"n2",&n2)) sf_error("No n2=.");

    if (!sf_histfloat(modl,"d1",&d1)) sf_error("No d1=.");
    if (!sf_histfloat(modl,"d2",&d2)) sf_error("No d2=.");

    vel = sf_floatalloc2(n1,n2);
    sf_floatread(vel[0],n1*n2,modl);

    /* read wavefield */
    if (!sf_histint(in,"n3",&ns)) sf_error("No ns= in input.");
    if (!sf_histint(in,"n4",&nw)) sf_error("No nw= in input.");
    if (!sf_histfloat(in,"d4",&dw)) sf_error("No dw= in input.");
    if (!sf_histfloat(in,"o4",&ow)) sf_error("No ow= in input.");

    wave = sf_complexalloc3(n1,n2,ns);
    
    /* allocate memory */
    rhs = sf_complexalloc3(n1,n2,ns);

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);

	sf_complexread(wave[0][0],n1*n2*ns,in);

	/* loop over shots */
#ifdef _OPENMP
#pragma omp parallel for num_threads(uts) private(j,i,neib,cent)
#endif
	for (is=0; is < ns; is++) {
	    switch (order[0]) {
		case '5':
		    for (j=0; j < n2; j++) {
			for (i=0; i < n1; i++) {
			    neib = 0.+I*0.;
			    cent = 0.+I*0.;
			    
			    /* left */
			    if (i > 0)    neib += wave[is][j][i-1]/(d1*d1);
			    cent += -1./(d1*d1);
			    /* right */
			    if (i < n1-1) neib += wave[is][j][i+1]/(d1*d1);
			    cent += -1./(d1*d1);
			    /* down */
			    if (j > 0)    neib += wave[is][j-1][i]/(d2*d2);
			    cent += -1./(d2*d2);
			    /* up */
			    if (j < n2-1) neib += wave[is][j+1][i]/(d2*d2);
			    cent += -1./(d2*d2);
			    /* center */
			    cent += pow(omega/vel[j][i],2.);
			    cent *= wave[is][j][i];
			    
			    rhs[is][j][i] = neib+cent;
			}
		    }
		    break;
		    
		case '9':
		    for (j=0; j < n2; j++) {
			for (i=0; i < n1; i++) {
			    neib = 0.+I*0.;
			    cent = 0.+I*0.;
			    
			    /* left left */
			    if (i > 1)    neib += (-1./12)*wave[is][j][i-2]/(d1*d1);
			    cent += -(-1./12)/(d1*d1);
			    /* left */
			    if (i > 0)    neib += (4./3)*wave[is][j][i-1]/(d1*d1);
			    cent += -(4./3)/(d1*d1);
			    /* right */
			    if (i < n1-1) neib += (4./3)*wave[is][j][i+1]/(d1*d1);
			    cent += -(4./3)/(d1*d1);
			    /* right right */
			    if (i < n1-2) neib += (-1./12)*wave[is][j][i+2]/(d1*d1);
			    cent += -(-1./12)/(d1*d1);
			    /* down down */
			    if (j > 1)    neib += (-1./12)*wave[is][j-2][i]/(d2*d2);
			    cent += -(-1./12)/(d2*d2);
			    /* down */
			    if (j > 0)    neib += (4./3)*wave[is][j-1][i]/(d2*d2);
			    cent += -(4./3)/(d2*d2);
			    /* up */
			    if (j < n2-1) neib += (4./3)*wave[is][j+1][i]/(d2*d2);
			    cent += -(4./3)/(d2*d2);
			    /* up up */
			    if (j < n2-2) neib += (-1./12)*wave[is][j+2][i]/(d2*d2);
			    cent += -(-1./12)/(d2*d2);
			    /* center */
			    cent += pow(omega/vel[j][i],2.);
			    cent *= wave[is][j][i];
			    
			    rhs[is][j][i] = neib+cent;
			}
		    }
		    break;
		    
		case 'j':
		    for (j=0; j < n2; j++) {
			for (i=0; i < n1; i++) {
			    neib = 0.+I*0.;
			    cent = 0.+I*0.;
			    
			    /* left */
			    if (i > 0)
				neib += (0.7926/(d1*d1)-
					 0.1037/(d2*d2)-
					 0.1037/(d2*d2)+
					 0.0942*pow(omega/vel[j][i-1],2.))*
				    wave[is][j][i-1];
			    cent += -0.7926/(d1*d1);
			    /* right */
			    if (i < n1-1)
				neib += (0.7926/(d1*d1)-
					 0.1037/(d2*d2)-
					 0.1037/(d2*d2)+
					 0.0942*pow(omega/vel[j][i+1],2.))*
				    wave[is][j][i+1];
			    cent += -0.7926/(d1*d1);
			    /* down */
			    if (j > 0)
				neib += (0.7926/(d2*d2)-
					 0.1037/(d1*d1)-
					 0.1037/(d1*d1)+
					 0.0942*pow(omega/vel[j-1][i],2.))*
				    wave[is][j-1][i];
			    cent += -0.7926/(d2*d2);
			    /* up */
			    if (j < n2-1)
				neib += (0.7926/(d2*d2)-
					 0.1037/(d1*d1)-
					 0.1037/(d1*d1)+
					 0.0942*pow(omega/vel[j+1][i],2.))*
				    wave[is][j+1][i];
			    cent += -0.7926/(d2*d2);
			    /* left down */
			    if (i > 0 && j > 0)
				neib += (0.1037/(d1*d1)+
					 0.1037/(d2*d2)-
					 0.0016*pow(omega/vel[j-1][i-1],2.))*
				    wave[is][j-1][i-1];
			    /* right up */
			    if (i < n1-1 && j < n2-1)
				neib += (0.1037/(d1*d1)+
					 0.1037/(d2*d2)-
					 0.0016*pow(omega/vel[j+1][i+1],2.))*
				    wave[is][j+1][i+1];
			    /* left up */
			    if (i > 0 && j < n2-1)
				neib += (0.1037/(d1*d1)+
					 0.1037/(d2*d2)-
					 0.0016*pow(omega/vel[j+1][i-1],2.))*
				    wave[is][j+1][i-1];
			    /* right down */
			    if (i < n1-1 && j > 0)
				neib += (0.1037/(d1*d1)+
					 0.1037/(d2*d2)-
					 0.0016*pow(omega/vel[j-1][i+1],2.))*
				    wave[is][j-1][i+1];
			    /* center */
			    cent += 0.6296*pow(omega/vel[j][i],2.);
			    cent *= wave[is][j][i];
			    
			    rhs[is][j][i] = neib+cent;
			}
		    }
		    break;
		    
		case 'c':
		    for (j=0; j < n2; j++) {
			for (i=0; i < n1; i++) {
			    neib = 0.+I*0.;
			    cent = 0.+I*0.;
			    
			    /* left up */
			    if (i > 0 && j > 0)
				neib += (0.06181325/(d1*d1)+
					 0.06181325/(d2*d2)+
					 0.0424801*pow(omega/vel[j-1][i-1],2.))*
				    wave[is][j-1][i-1];
			    /* up */
			    if (i > 0)
				neib += (0.2880195/(d1*d1)-
					 0.06181325/(d2*d2)-
					 0.06181325/(d2*d2)-
					 0.1389664/(4.*d2*d2)-
					 0.1389664/(4.*d2*d2)+
					 0.108598*pow(omega/vel[j][i-1],2.))*
				    wave[is][j][i-1];
			    cent += -0.2880195/(d1*d1);
			    /* right up */
			    if (i > 0 && j < n2-1)
				neib += (0.06181325/(d1*d1)+
					 0.06181325/(d2*d2)+
					 0.0424801*pow(omega/vel[j+1][i-1],2.))*
				    wave[is][j+1][i-1];
			    /* right */
			    if (j < n2-1)
				neib += (0.2880195/(d2*d2)-
					 0.06181325/(d1*d1)-
					 0.06181325/(d1*d1)-
					 0.1389664/(4.*d1*d1)-
					 0.1389664/(4.*d1*d1)+
					 0.108598*pow(omega/vel[j+1][i],2.))*
				    wave[is][j+1][i];
			    cent += -0.2880195/(d2*d2);
			    /* right down */
			    if (i < n1-1 && j < n2-1)
				neib += (0.06181325/(d1*d1)+
					 0.06181325/(d2*d2)+
					 0.0424801*pow(omega/vel[j+1][i+1],2.))*
				    wave[is][j+1][i+1];
			    /* down */
			    if (i < n1-1)
				neib += (0.2880195/(d1*d1)-
					 0.06181325/(d2*d2)-
					 0.06181325/(d2*d2)-
					 0.1389664/(4.*d2*d2)-
					 0.1389664/(4.*d2*d2)+
					 0.108598*pow(omega/vel[j][i+1],2.))*
				    wave[is][j][i+1];
			    cent += -0.2880195/(d1*d1);
			    /* left down */
			    if (i < n1-1 && j > 0)
				neib += (0.06181325/(d1*d1)+
					 0.06181325/(d2*d2)+
					 0.0424801*pow(omega/vel[j-1][i+1],2.))*
				    wave[is][j-1][i+1];
			    /* left */
			    if (j > 0)
				neib += (0.2880195/(d2*d2)-
					 0.06181325/(d1*d1)-
					 0.06181325/(d1*d1)-
					 0.1389664/(4.*d1*d1)-
					 0.1389664/(4.*d1*d1)+
					 0.108598*pow(omega/vel[j-1][i],2.))*
				    wave[is][j-1][i];
			    cent += -0.2880195/(d2*d2);
			    /* left left up up */
			    if (i > 1 && j > 1)
				neib += (0.007436025/(4.*d1*d1)+
					 0.007436025/(4.*d2*d2)+
					 0.000206312*pow(omega/vel[j-2][i-2],2.))*
				    wave[is][j-2][i-2];
			    /* left up up */
			    if (i > 1 && j > 0)
				neib += (0.1389664/(4.*d1*d1)+
					 0.00188342*pow(omega/vel[j-1][i-2],2.))*
				    wave[is][j-1][i-2];
			    /* up up */
			    if (i > 1)
				neib += (0.29554905/(4.*d1*d1)-
					 0.007436025/(4.*d2*d2)-
					 0.007436025/(4.*d2*d2)+
					 0.0041487*pow(omega/vel[j][i-2],2.))*
				    wave[is][j][i-2];
			    cent += -0.29554905/(4.*d1*d1);
			    /* right up up */
			    if (i > 1 && j < n2-1)
				neib += (0.1389664/(4.*d1*d1)+
					 0.00187765*pow(omega/vel[j+1][i-2],2.))*
				    wave[is][j+1][i-2];
			    /* right right up up */
			    if (i > 1 && j < n2-2)
				neib += (0.007436025/(4.*d1*d1)+
					 0.007436025/(4.*d2*d2)+
					 0.000206312*pow(omega/vel[j+2][i-2],2.))*
				    wave[is][j+2][i-2];
			    /* right right up */
			    if (i > 0 && j < n2-2)
				neib += (0.1389664/(4.*d2*d2)+
					 0.00188342*pow(omega/vel[j+2][i-1],2.))*
				    wave[is][j+2][i-1];
			    /* right right */
			    if (j < n2-2)
				neib += (0.29554905/(4.*d2*d2)-
					 0.007436025/(4.*d1*d1)-
					 0.007436025/(4.*d1*d1)+
					 0.0041487*pow(omega/vel[j+2][i],2.))*
				    wave[is][j+2][i];
			    cent += -0.29554905/(4.*d2*d2);
			    /* right right down */
			    if (i < n1-1 && j < n2-2)
				neib += (0.1389664/(4.*d2*d2)+
					 0.00187765*pow(omega/vel[j+2][i+1],2.))*
				    wave[is][j+2][i+1];
			    /* right right down down */
			    if (i < n1-2 && j < n2-2)
				neib += (0.007436025/(4.*d1*d1)+
					 0.007436025/(4.*d2*d2)+
					 0.000206312*pow(omega/vel[j+2][i+2],2.))*
				    wave[is][j+2][i+2];
			    /* right down down */
			    if (i < n1-2 && j < n2-1)
				neib += (0.1389664/(4.*d1*d1)+
					 0.00188342*pow(omega/vel[j+1][i+2],2.))*
				    wave[is][j+1][i+2];
			    /* down down */
			    if (i < n1-2)
				neib += (0.29554905/(4.*d1*d1)-
					 0.007436025/(4.*d2*d2)-
					 0.007436025/(4.*d2*d2)+
					 0.0041487*pow(omega/vel[j][i+2],2.))*
				    wave[is][j][i+2];
			    cent += -0.29554905/(4.*d1*d1);
			    /* left down down */
			    if (i < n1-2 && j > 0)
				neib += (0.1389664/(4.*d1*d1)+
					 0.00187765*pow(omega/vel[j-1][i+2],2.))*
				    wave[is][j-1][i+2];
			    /* left left down down */
			    if (i < n1-2 && j > 1)
				neib += (0.007436025/(4.*d1*d1)+
					 0.007436025/(4.*d2*d2)+
					 0.000206312*pow(omega/vel[j-2][i+2],2.))*
				    wave[is][j-2][i+2];
			    /* left left down */
			    if (i < n1-1 && j > 1)
				neib += (0.1389664/(4.*d2*d2)+
					 0.00188342*pow(omega/vel[j-2][i+1],2.))*
				    wave[is][j-2][i+1];
			    /* left left */
			    if (j > 1)
				neib += (0.29554905/(4.*d2*d2)-
					 0.007436025/(4.*d1*d1)-
					 0.007436025/(4.*d1*d1)+
					 0.0041487*pow(omega/vel[j-2][i],2.))*
				    wave[is][j-2][i];
			    cent += -0.29554905/(4.*d2*d2);
			    /* left left up */
			    if (i > 0 && j > 1)
				neib += (0.1389664/(4.*d2*d2)+
					 0.00187765*pow(omega/vel[j-2][i-1],2.))*
				    wave[is][j-2][i-1];
			    /* center */
			    cent += 0.363276*pow(omega/vel[j][i],2.);
			    cent *= wave[is][j][i];
			    
			    rhs[is][j][i] = neib+cent;
			}
		    }
		    break;
		    
		default:
		    sf_error("Fail to load discretization scheme.");
	    }
	}

	sf_complexwrite(rhs[0][0],n1*n2*ns,out);
    }
    
    exit(0);
}
