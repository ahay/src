/* 3-D Kirchhoff zero-offset modeling/migration antialiased by parameterization */
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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void mig3_lop (bool adj /* adjoint flag */,
	       bool add /* addition flag - for x0! */,
               int nt, int nx, int ny /* data size */,
               float dt, float dx, float dy,
	       float ot, float ox, float oy,
               float *trace /* zero-offset */,
               float *out /* image */,
               float vel /* velocity */,
               float rho /* phase factor */,
               char antialias /* antialiasing method */,
	       bool doomp /* perform OpenMP optimization */,
	       int apt /* distance aperture in samples */,
	       float angle /* angle aperture in degrees */)
/*< Apply >*/
{

    int nn, iy, ix, iz, iyc, ixc, izc, indexc;
    int mythread, nothreads, ch=0, chomp=0, deltat, ch2=0, ch3=0, ch4=0;
    float *pp, *qq, *trace_temp;
    float x_in, y_in, ftm, ftp, imp, xtemp, ytemp, veltemp;
    float sq, ti, tx, tm, tp, x,y, z, t;
    float rx, ry, ft;
    int it, itp, itm, itrx, itry, index, checkomp=0;

    /* add flag is enabled when starting model is given */
    /* refer to api/c/bigsolver.c 1392: "line oper (false, true, nx, ny, x, rr);"*/
    sf_adjnull(adj,add,nt*nx*ny,nt*nx*ny,out,trace);

    /* looking for artifacts */
    izc = 225;
    ixc = 175;
    iyc = 0;
    indexc = iyc*nx*nt + ixc*nt + izc;

/* _____________________________________________________________________________________________________________________________________________________ */
//sf_warning("checking sample (%d,%d,%d)",izc,ixc,iyc);
//if (!ch2 && (fabsf(trace[indexc])>0.) ){ sf_warning("adjnull: adj=%d add=%d trace[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,trace[indexc]);}
//if (!ch3 && (fabsf(out[indexc])>0.) ){ sf_warning("adjnull: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]);}
/* _____________________________________________________________________________________________________________________________________________________ */

    
    trace_temp = sf_floatalloc(nt);

    for (iz=0; iz < nt; iz++) {

	trace_temp[iz] = 0.;

    }
    
    veltemp = vel/2.0;

    vel = 2./vel; 
    vel *= vel;

    angle = fabsf(tanf(angle*SF_PI/180.0));

    /* full differentiation variables */
    //nn = 2*kiss_fft_next_fast_size((nt+1)/2);
    //pp = sf_floatalloc(nn);

    /* igrad */
    pp = sf_floatalloc(nt);
    qq = sf_floatalloc(nt);
    /* sf_halfint_init (true, nn, rho); */
    /* is repeated each time just in case */ 

    for (iz=0; iz < nt /*nn*/; iz++) {
	pp[iz] = 0.;
	qq[iz] = 0.;
    }

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch3 && (fabsf(out[indexc])>0.) && adj ){ sf_warning("out halfint start: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]);}
/* ________________________________________________________________________________________________________________________________________________________ */

    for (itry=0; itry < ny; itry++) {
    for (itrx=0; itrx < nx; itrx++) {

	x_in = ox + itrx*dx;
	y_in = oy + itry*dy;

	/* if adj - we actually copy traces */
	/* if fwd - we just copy zeroes */
	/* however when add=true (when starting model is given) */
	/* it will be properly taken into account  */
	/* refer to api/c/bigsolver.c 1392: "line oper (false, true, nx, ny, x, rr);" */
	for(iz=0; iz<nt; iz++){

		if (adj){		
			trace_temp[iz] = trace[itry*nx*nt + itrx*nt + iz];
		}
		if (!adj){
			trace_temp[iz] = 0.;
		}

	}/* iz */

	if (adj){

		/* full-differentiation in the data domain */
		for (iz=0; iz < nt; iz++) {
	
		    		pp[iz] = trace_temp[iz];

		}

		sf_igrad1_lop (true,false,nt,nt,qq,pp);

		for (iz=0; iz < nt; iz++) {
		
			trace_temp[iz] = qq[iz];

		}

		if ('t'==antialias) { 
        	    	
			sf_doubint(true,nt,trace_temp);

		}

	}/* if */

	sf_warning("inline=%d out of %d, crossline=%d out of %d;",itrx,nx,itry,ny);

	if(doomp) {//do OMP flag

		#ifdef _OPENMP
		checkomp = 1;
		#endif

		if(checkomp !=1) sf_error("OMP is not supported");

		#pragma omp parallel
		{

			mythread = omp_get_thread_num();
    			nothreads = omp_get_num_threads();

			#pragma omp single
			{
    			if(!chomp) sf_warning("Using %d threads",nothreads); chomp=1;
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
	
			#pragma omp for private(iy,y,ry,ix,x,rx,iz,z,t,sq,ti,tx,ft,it,tm,ftm,itm,tp,ftp,itp,imp,index,xtemp,ytemp,ch)
			for (iy=0; iy < ny; iy++) {
	    			y = oy + iy*dy - y_in;
	    			ry = fabsf(y*dy) * vel;
	    			y *= y * vel;

				if (SF_ABS(iy-itry) > apt) continue;
	    
            			for (ix=0; ix < nx; ix++) {
					x = ox + ix*dx - x_in;
					rx = fabsf(x*dx) * vel;
					x *= x * vel;
					x += y;

					if (SF_ABS(ix-itrx) > apt) continue;

					for (iz=0; iz < nt; iz++) {
		    				z = ot + iz*dt;              

						xtemp = (ix - itrx)*dx;
						ytemp = (iy - itry)*dy;
	
						if ((fabsf(xtemp) > angle*veltemp*z) || (fabsf(ytemp) > angle*veltemp*z)) continue;		    				
				
						t = z*z + x; 
		    				if (t <= 0.) continue;
              
		    				sq = t*t; 
		    				if (sq <= 0.) continue;

		    				sq = sqrtf(sq);
		    				ti = sqrtf(0.5*(t + sq));

		    				tx = SF_MAX(
						rx / ti * 0.5 * (1 + (t)/sq), 
						ry * ti / sq);

		    				if ('f'==antialias && tx > dt) continue;
		    
		    				ft = (ti-ot)/dt; it = floorf(ft); ft -= it; 

		    				if ( it < 0) continue;
		    				if ( it >= nt-1) break; 
              
		    				if ('t'==antialias) {

							tm = ti-tx-dt;
							ftm = (tm-ot)/dt; itm = floorf(ftm); ftm = ftm - itm; 
							if (itm < 0) continue;
                 
							tp = ti+tx+dt;
							ftp = (tp-ot)/dt; itp = floorf(ftp); ftp = ftp - itp; 
							if (itp >= nt-1) continue;

							imp = dt/(dt+tp-tm);
							imp *= imp;

							index = iy*nx*nt + ix*nt + iz;

							if (adj) {

								out[index] += imp*(
			    					2.*(1.-ft)*trace_thread[it] + 2.*ft*trace_thread[it+1] -
			    					(1.-ftm)*trace_thread[itm] - ftm*trace_thread[itm+1]    - 
			    					(1.-ftp)*trace_thread[itp] - ftp*trace_thread[itp+1]);

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch4 && (fabsf(out[indexc])>0.) ){ sf_warning("out+=tracethr: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]); ch4=1; }
/* ________________________________________________________________________________________________________________________________________________________ */

							} else {

								trace_thread[it]    += imp*2.*(1.-ft)*out[index];
								trace_thread[it+1]  += imp*2.*ft*out[index];
								trace_thread[itm]   -= imp*(1.-ftm)*out[index];
								trace_thread[itm+1] -= imp*ftm*out[index];
								trace_thread[itp]   -= imp*(1.-ftp)*out[index];
								trace_thread[itp+1] -= imp*ftp*out[index];

/* __________________________________________________________________________________________________________________________________________________________ */	
//if (!ch4 && !adj && (fabsf(trace_thread[225])>0.) && (index==indexc) ){ sf_warning("tracethr: adj=%d add=%d trace[225]=%f\n",adj,add,trace[225]); ch4=1; }
/* __________________________________________________________________________________________________________________________________________________________ */
							
							}

		    				} else {

							index = iy*nx*nt + ix*nt + iz;

							if (adj){

								out[index] +=
			  				  	(1.-ft)*trace_thread[it] + ft*trace_thread[it+1];

							} else {

								trace_thread[it] += (1.-ft)*out[index];
								trace_thread[it+1] += ft*out[index];

							}

		    				} /* antialiasing flag */
					} /* iz */
	    			} /* ix */
			} /* iy */

			if (adj==false) {
				#pragma omp critical
				{

				int jthr;

				for(jthr = 0; jthr < nt; jthr++) {
    						 
					trace_temp[jthr] += trace_thread[jthr];

					if( ((fabsf(trace_thread[jthr]) > 3.0e+38)  || (isinf(trace_thread[jthr])) || (isnan(trace_thread[jthr]))) && !ch )
						{ 
							sf_warning("mig3 adj=false: problems in zero-offset OMP reduction"); ch=1;
						}
				
					if( ((fabsf(trace_temp[jthr]) > 3.0e+38)  || (isinf(trace_temp[jthr])) || (isnan(trace_temp[jthr]))) && !ch)
						{
							sf_warning("mig3 adj=false: problems in zero-offset OMP reduction"); ch=1;
						}

  				}/* jthr */
				}/* critical */

			}/* adj */

			free(trace_thread);

		}/* omp region */

	}/* doomp */ else {

			for (iy=0; iy < ny; iy++) {
	    			y = oy + iy*dy - y_in;
	    			ry = fabsf(y*dy) * vel;
	    			y *= y * vel;
	    
            			for (ix=0; ix < nx; ix++) {
					x = ox + ix*dx - x_in;
					rx = fabsf(x*dx) * vel;
					x *= x * vel;
					x += y;

					for (iz=0; iz < nt; iz++) {
		    				z = ot + iz*dt;              

		    				t = z*z + x; 
		    				if (t <= 0.) continue;
              
		    				sq = t*t; 
		    				if (sq <= 0.) continue;

		    				sq = sqrtf(sq);
		    				ti = sqrtf(0.5*(t + sq));

		    				tx = SF_MAX(
						rx / ti * 0.5 * (1 + (t)/sq), 
						ry * ti / sq);

		    				if ('f'==antialias && tx > dt) continue;
		    
		    				ft = (ti-ot)/dt; it = floorf(ft); ft -= it; 

		    				if ( it < 0) continue;
		    				if ( it >= nt-1) break; 
              
		    				if ('t'==antialias) {
							tm = ti-tx-dt;
							ftm = (tm-ot)/dt; itm = floorf(ftm); ftm = ftm - itm; 
							if (itm < 0) continue;
                 
							tp = ti+tx+dt;
							ftp = (tp-ot)/dt; itp = floorf(ftp); ftp = ftp - itp; 
							if (itp >= nt-1) continue;

							imp = dt/(dt+tp-tm);
							imp *= imp;

							index = iy*nx*nt + ix*nt + iz;

							if (adj) {
			
								out[index] += imp*(
			    					2.*(1.-ft)*trace_temp[it] + 2.*ft*trace_temp[it+1] -
			    					(1.-ftm)*trace_temp[itm] - ftm*trace_temp[itm+1]    - 
			    					(1.-ftp)*trace_temp[itp] - ftp*trace_temp[itp+1]);

							} else {

								trace_temp[it]    += imp*2.*(1.-ft)*out[index];
								trace_temp[it+1]  += imp*2.*ft*out[index];
								trace_temp[itm]   -= imp*(1.-ftm)*out[index];
								trace_temp[itm+1] -= imp*ftm*out[index];
								trace_temp[itp]   -= imp*(1.-ftp)*out[index];
								trace_temp[itp+1] -= imp*ftp*out[index];
	
							}

		    				} else {

							index = iy*nx*nt + ix*nt + iz;

							if (adj){

								out[index] +=
			  				  	(1.-ft)*trace_temp[it] + ft*trace_temp[it+1];

							} else {

								trace_temp[it] += (1.-ft)*out[index];
								trace_temp[it+1] += ft*out[index];

							}

		    				} /* antialiasing flag */
					} /* iz */
	    			} /* ix */
			} /* iy */	

	}/* no omp */

    	if (adj == false){

		if ('t'==antialias) { 
        	    				
			sf_doubint(true,nt,trace_temp);

		}/* antialias */

		/* full-differentiation in the data domain */
		for (iz=0; iz < nt; iz++) {
	
		    	pp[iz] = trace_temp[iz];

		}

		sf_igrad1_lop (false,false,nt,nt,pp,qq);

		for (iz=0; iz < nt; iz++) {
		
			trace_temp[iz] = qq[iz];

		}

		/* with add flag we keep data modeled from the initial guess */
		/* here newly modeled data is added to it */
 		for(iz=0; iz<nt; iz++){

			trace[itry*nx*nt + itrx*nt + iz] += trace_temp[iz];

		}/* iz */

	}/* if */

    } }/* itrx and itry */

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch2 && !adj && (fabsf(trace[indexc])>0.) ){ sf_warning("trace+=tracetemp: adj=%d add=%d trace[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,trace[indexc]);}
/* ________________________________________________________________________________________________________________________________________________________ */

/* ________________________________________________________________________________________________________________________________________________________ */
//if (!ch3 && (fabsf(out[indexc])>0.) ){ sf_warning("out halfint final: adj=%d add=%d out[%d,%d,%d]=%f\n",adj,add,izc,ixc,iyc,out[indexc]);}
/* ________________________________________________________________________________________________________________________________________________________ */

    free(trace_temp);

}
