/* 2-D Kirchhoff zero-offset modeling/migration antialiased by parameterization */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
  Foundation, Inc., 59 Temple Place, Suite 330, Bston, MA  02111-1307  USA
*/

#include <math.h>

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void mig2_lop (bool adj /* adjoint flag */,
	       bool add /* addition flag - for x0! */,
               int nt, int nx /* data size */,
               float dt, float dx,
	       float ot, float ox,
               float *trace /* zero-offset */,
               float *out /* image */,
               float *v /* velocity */,
               float rho /* phase factor */,
               bool hd /* half derivative */,
               char antialias /* antialiasing method */,
	       bool doomp /* perform OpenMP optimization */,
	       int apt /* distance aperture in samples */,
	       float angle /* angle aperture in degrees */,
               bool ps /* spherical divergence */,
               bool dd /* if perform differentiation in the data domain */)
/*< Apply >*/
{

    int nn, ix, iz;
    int nothreads, ch=0, chomp=0;
    float *pp, *trace_temp;
    float x_in, ftm, ftp, imp, xtemp, veltemp, vel, amp;
    float sq, ti, tx, tm, tp, x, z, t;
    float rx, ft;
    int it, itp, itm, itrx, index, checkomp=0;

    /* add flag is enabled when starting model is given */
    /* refer to api/c/bigsolver.c 1392: "line oper (false, true, nx, ny, x, rr);"*/
    sf_adjnull(adj,add,nt*nx,nt*nx,out,trace);

/* _____________________________________________________________________________________________________________________________________________________ */
//sf_warning("checking sample (%d,%d,%d)",izc,ixc,iyc);
//if (!ch2 && (fabsf(trace[indexc])>0.) ){ sf_warning("adjnull: adj=%d add=%d trace[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,trace[indexc]);}
//if (!ch3 && (fabsf(out[indexc])>0.) ){ sf_warning("adjnull: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]);}
/* _____________________________________________________________________________________________________________________________________________________ */

    
    trace_temp = sf_floatalloc(nt);

    for (iz=0; iz < nt; iz++) {

	trace_temp[iz] = 0.;

    }

    angle = fabsf(tanf(angle*SF_PI/180.0));

    /* half differentiation variables */
    nn = 2*kiss_fft_next_fast_size((nt+1)/2);
    //pp = sf_floatalloc(nt);
    pp = sf_floatalloc(nn);

    /*^^^^^*/
    //sf_halfint_init (true, nt, rho);
    sf_halfint_init (true, nn, rho);

    if (!adj) {

	for (itrx=0; itrx < nx; itrx++) {

		/* half differentiation in the model domain */
		if (hd && !dd) {
			for (iz=0; iz < nt; iz++) {
	
				index = itrx*nt + iz;
	
				pp[iz] = out[index];

			}

			for (iz=nt; iz < nn; iz++) {
				pp[iz] = 0.;
	        	}

			sf_halfint (false, pp);

			for (iz=0; iz < nt; iz++) {

				index = itrx*nt + iz;
		
				out[index] = pp[iz];

			}

		}/* hd */

	}

    }

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch3 && (fabsf(out[indexc])>0.) && adj ){ sf_warning("out halfint start: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]);}
/* ________________________________________________________________________________________________________________________________________________________ */

    for (itrx=0; itrx < nx; itrx++) {

	x_in = ox + itrx*dx;

	/* if adj - we actually copy traces */
	/* if fwd - we just copy zeroes */
	/* however when add=true (when starting model is given) */
	/* it will be properly taken into account  */
	/* refer to api/c/bigsolver.c 1392: "line oper (false, true, nx, ny, x, rr);" */
	for(iz=0; iz<nt; iz++){

		if (adj){		
			trace_temp[iz] = trace[itrx*nt + iz];
		}
		if (!adj){
			trace_temp[iz] = 0.;
		}

	}/* iz */

	if (adj){

		/* half differentiation in the data domain */
		if (hd && dd) {

			for (iz=0; iz < nt; iz++) {
	
				pp[iz] = trace_temp[iz];

			}

			for (iz=nt; iz < nn; iz++) {
				pp[iz] = 0.;
		        }

			sf_halfint (true, pp);

			for (iz=0; iz < nt; iz++) {
		
				trace_temp[iz] = pp[iz];

			}

		}/* hd */

		if ('t'==antialias) { 
        	    	
			sf_doubint(true,nt,trace_temp);

		}

	}/* if */

	sf_warning("inline=%d out of %d;",itrx,nx);

	/*if(doomp) {do OMP flag // working progress to account for OMP defficiency */

		#ifdef _OPENMP
		checkomp = 1;
		#endif

		if(checkomp !=1) sf_error("OMP is not supported");

		#pragma omp parallel
		{
#ifdef _OPENMP
    			nothreads = omp_get_num_threads();
#else
			nothreads = 1;
#endif			

			#pragma omp single
			{
    			if(!chomp) {
					sf_warning("Using %d threads",nothreads); 
					chomp=1;
				}
			}

			float *trace_thread;
			int ithr;

			trace_thread = sf_floatalloc(nt);

			for (ithr = 0; ithr < nt; ithr++){

				if (adj){

					trace_thread[ithr] = trace_temp[ithr];

				} else {/* we do have a reduction clause */
					/* if we set "trace_thread[ithr] = trace_temp[ithr];" */
					/* for adj=false - starting model data will be summed twice */

					trace_thread[ithr] = 0.;

				}

			}/* for */
	
			#pragma omp for private(ix,x,rx,iz,z,t,sq,ti,tx,ft,it,tm,ftm,itm,tp,ftp,itp,imp,index,xtemp,ch,veltemp,vel,amp)
			for (ix=0; ix < nx; ix++) {
				
				x = ox + ix*dx - x_in;

				if (SF_ABS(ix-itrx) > apt) continue;

				for (iz=0; iz < nt; iz++) {
		    			z = ot + iz*dt;              

					veltemp = v[ix*nt + iz]/2.0;
					vel = 2./v[ix*nt + iz]; 
					vel *= vel;

					rx = fabsf(x*dx) * vel;
					/* x *= x * vel; // moved to t */

					xtemp = (ix - itrx)*dx;
	
					if ((fabsf(xtemp) > angle*veltemp*z)) continue;		    				
				 
					t = z*z + x*x*vel; 
	    				if (t <= 0.) continue;
              
	    				sq = t*t; 
	    				if (sq <= 0.) continue;

	    				sq = sqrtf(sq);
	    				ti = sqrtf(0.5*(t + sq));

	    				tx = rx / ti * 0.5 * (1 + (t)/sq);

		    			if ('f'==antialias && tx > dt) continue;
		    
		    			ft = (ti-ot)/dt; it = floorf(ft); ft -= it; 

		    			if ( it < 0) continue;
		    			if ( it >= nt-1) break; 

					/* spherical divergence 2D */
					amp = ps ? (z/(ti+dt))*sqrtf(nt*dt / (ti+dt)) : 1.;              

		    			if ('t'==antialias) {

						tm = ti-tx-dt;
						ftm = (tm-ot)/dt; itm = floorf(ftm); ftm = ftm - itm; 
						if (itm < 0) continue;
                 
						tp = ti+tx+dt;
						ftp = (tp-ot)/dt; itp = floorf(ftp); ftp = ftp - itp; 
						if (itp >= nt-1) continue;

						imp = dt/(dt+tp-tm);
						imp *= imp;

						index = ix*nt + iz;

						if (adj) {

							out[index] += imp*amp*(
			    				2.*(1.-ft)*trace_thread[it] + 2.*ft*trace_thread[it+1] -
			    				(1.-ftm)*trace_thread[itm] - ftm*trace_thread[itm+1]    - 
			    				(1.-ftp)*trace_thread[itp] - ftp*trace_thread[itp+1]);

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch4 && (fabsf(out[indexc])>0.) ){ sf_warning("out+=tracethr: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]); ch4=1; }
/* ________________________________________________________________________________________________________________________________________________________ */

						} else {

							trace_thread[it]    += amp*imp*2.*(1.-ft)*out[index];
							trace_thread[it+1]  += amp*imp*2.*ft*out[index];
							trace_thread[itm]   -= amp*imp*(1.-ftm)*out[index];
							trace_thread[itm+1] -= amp*imp*ftm*out[index];
							trace_thread[itp]   -= amp*imp*(1.-ftp)*out[index];
							trace_thread[itp+1] -= amp*imp*ftp*out[index];

/* __________________________________________________________________________________________________________________________________________________________ */	
//if (!ch4 && !adj && (fabsf(trace_thread[225])>0.) && (index==indexc) ){ sf_warning("tracethr: adj=%d add=%d trace[225]=%f\n",adj,add,trace[225]); ch4=1; }
/* __________________________________________________________________________________________________________________________________________________________ */
							
						}

		    			} else {

						index = ix*nt + iz;

						if (adj){

							out[index] +=
			  				amp*(1.-ft)*trace_thread[it] + amp*ft*trace_thread[it+1];

						} else {

							trace_thread[it] += amp*(1.-ft)*out[index];
							trace_thread[it+1] += amp*ft*out[index];

						}

		    			} /* antialiasing flag */
				} /* iz */
	    		} /* ix */

			if (adj==false) {
				#pragma omp critical
				{

				int jthr;

				for(jthr = 0; jthr < nt; jthr++) {
    						 
					trace_temp[jthr] += trace_thread[jthr];

					if( ((fabsf(trace_thread[jthr]) > 3.0e+38)  || (isinf(trace_thread[jthr])) || (isnan(trace_thread[jthr]))) && !ch )
						{ 
							sf_warning("mig2 adj=false: problems in zero-offset OMP reduction inline=%d",itrx); ch=1;
							sf_warning("mig2 adj=false: trace_thread=%f",trace_thread[jthr]); ch=1;
						}
				
					if( ((fabsf(trace_temp[jthr]) > 3.0e+38)  || (isinf(trace_temp[jthr])) || (isnan(trace_temp[jthr]))) && !ch)
						{
							sf_warning("mig2 adj=false: problems in zero-offset OMP reduction inline=%d",itrx); ch=1;
							sf_warning("mig2 adj=false: trace_temp=%f",trace_temp[jthr]); ch=1;
						}

  				}/* jthr */
				}/* critical */

			}/* adj */

			free(trace_thread);

		}/* omp region */

	/*} doomp */

    	if (adj == false){

		if ('t'==antialias) { 
        	    				
			sf_doubint(true,nt,trace_temp);

		}/* antialias */

		/* half differentiation in the data domain */
		if (hd && dd) {
			for (iz=0; iz < nt; iz++) {
	
				pp[iz] = trace_temp[iz];

			}

			for (iz=nt; iz < nn; iz++) {
				pp[iz] = 0.;
	        	}

			sf_halfint (false, pp);

			for (iz=0; iz < nt; iz++) {
		
				trace_temp[iz] = pp[iz];

			}

		}/* hd */

		/* with add flag we keep data modeled from the initial guess */
		/* here newly modeled data is added to it */
 		for(iz=0; iz<nt; iz++){

			trace[itrx*nt + iz] += trace_temp[iz];

		}/* iz */

	}/* if */


    }/* itrx */

    if (adj) {

	for (itrx=0; itrx < nx; itrx++) {

		/* half differentiation in the model domain */
		if (hd && !dd) {

			for (iz=0; iz < nt; iz++) {
	
				index = itrx*nt + iz;
	
				pp[iz] = out[index];

			}

			for (iz=nt; iz < nn; iz++) {
				pp[iz] = 0.;
	        	}

			sf_halfint (adj, pp);

			for (iz=0; iz < nt; iz++) {

				index = itrx*nt + iz;
		
				out[index] = pp[iz];

			}

		}/* hd */

	}

    }

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch2 && !adj && (fabsf(trace[indexc])>0.) ){ sf_warning("trace+=tracetemp: adj=%d add=%d trace[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,trace[indexc]);}
/* ________________________________________________________________________________________________________________________________________________________ */

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch3 && (fabsf(out[indexc])>0.) ){ sf_warning("out halfint final: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]);}
/* ________________________________________________________________________________________________________________________________________________________ */

    free(trace_temp);

}
