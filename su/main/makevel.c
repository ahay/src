/* Make a velocity function v(x,y,z) */
/*
  Copyright (c) 2007, Colorado School of Mines,
  All rights reserved.
  
  
  Redistribution and use in source and binary forms, with or 
  without modification, are permitted provided that the following 
  conditions are met:
  
  *  Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
  *  Redistributions in binary form must reproduce the above 
  copyright notice, this list of conditions and the following 
  disclaimer in the documentation and/or other materials provided 
  with the distribution.
  *  Neither the name of the Colorado School of Mines nor the names of
  its contributors may be used to endorse or promote products 
  derived from this software without specific prior written permission.
  
  Warranty Disclaimer:
  THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
  COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
  POSSIBILITY OF SUCH DAMAGE.
  
  
  Export Restriction Disclaimer:
  We believe that CWP/SU: Seismic Un*x is a low technology product that does
  not appear on the Department of Commerce CCL list of restricted exports.
  Accordingly, we believe that our product meets the qualifications of
  an ECCN (export control classification number) of EAR99 and we believe
  it fits the qualifications of NRR (no restrictions required), and
  is thus not subject to export restrictions of any variety.
  
  Approved Reference Format:
  In publications, please refer to SU as per the following example:
  Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
  Release No. __: an open source software  package for seismic 
  research and processing, 
  Center for Wave Phenomena, Colorado School of Mines.
  
  Articles about SU in peer-reviewed journals:
  Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
  Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
  Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
  Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
  
  Acknowledgements:
  SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
  School of Mines, partially based on Stanford Exploration Project (SEP) 
  software.
*/
/* Modified by Tariq Alkhalifah for inclusion with Madagascar. */

#include <math.h>

#include <rsf.h>
#include <rsf_su.h>

/*
 * Extracted from CWP Author: Dave Hale
 *
 */

int
main (int argc, char **argv)
{
    int nx,ny,nz,ix,iy,iz;
    float dx,dy,dz,fx,fy,fz,x,y,z,v000,dvdx,dvdy,dvdz,
	xlens,ylens,zlens,dlens,tlens,vlens,xn,ynn,zn,
	/* abot,bbot,atop,btop, */
	vran,vzran, 
	vzc,z1c,z2c,l1c,l2c,exc,ac,bc,vtemp,
	*v,*vz;
    /* float zbot,ztop;*/
    char *vzfile;
    sf_file out, vzf;

    /* hook up getpar to handle the parameters */
    sf_init (argc, argv);

    out = sf_output("out");

    if (NULL != (vzfile = sf_getstring("vzfile"))) {
	vzf = sf_input(vzfile);
	free(vzfile);
    } else {
	vzf = NULL;
    }
	
    /* get required parameters */
    if (!sf_getint("n2",&nx)) sf_error("must specify n2!\n");
    /* number of x samples (2nd dimension), must be provided!*/
    if (!sf_getint("n1",&nz)) sf_error("must specify n1!\n");
    /* number of z samples (1st dimension)), must be provided!*/
	
    /* get optional parameters */
    if (!sf_getint("n3",&ny)) ny = 1;
    /* number of y samples (3rd dimension)*/
    if (!sf_getfloat("d2",&dx)) dx = 1.0;
    /* 2nd dimension sampling interval*/
    if (!sf_getfloat("d3",&dy)) dy = 1.0;
    /* 3rd dimension sampling interval*/
    if (!sf_getfloat("d1",&dz)) dz = 1.0;
    /* 1st dimension sampling interval*/
    if (!sf_getfloat("o2",&fx)) fx = 0.0;
    /* Origin 2nd dimension*/
    if (!sf_getfloat("o3",&fy)) fy = 0.0;
    /* Origin 3rd dimension*/
    if (!sf_getfloat("o1",&fz)) fz = 0.0;
    /* Origin 1st dimension*/
    if (!sf_getfloat("v000",&v000)) v000 = 2.0;
    /* velocity at (x=0,y=0,z=0)*/
    if (!sf_getfloat("dvdx2",&dvdx)) dvdx = 0.0;
    /* velocity gradient with respect to 2nd dimension*/
    if (!sf_getfloat("dvdx3",&dvdy)) dvdy = 0.0;
    /* velocity gradient with respect to 3rd dimension*/
    if (!sf_getfloat("dvdx1",&dvdz)) dvdz = 0.0;
    /* velocity gradient with respect to 1st dimension*/
    if (!sf_getfloat("x2lens",&xlens)) xlens = fx;
    /* 2nd dimension coordinate of center of parabolic lens*/
    if (!sf_getfloat("x3lens",&ylens)) ylens = fy;
    /* 3rd dimension coordinate of center of parabolic lens*/
    if (!sf_getfloat("x1lens",&zlens)) zlens = fz;
    /* 1st dimension coordinate of center of parabolic lens*/
    if (!sf_getfloat("vlens",&vlens)) vlens = 0.0;
    /* velocity perturbation in parabolic lens*/
    if (!sf_getfloat("dlens",&dlens)) dlens = 1.0;
    /* diameter of parabolic lens*/
    if (!sf_getfloat("tlens",&tlens)) tlens = 1.0;
    /* thickness of parabolic lens*/
    if (!sf_getfloat("vran",&vran)) vran = 0.0;
    /* standard deviation of random perturbation*/
    vzfile = sf_getstring("vx1file");
    /* file containing v(z) 1st dimension profile*/
    if (!sf_getfloat("vx1ran",&vzran)) vzran = 0.0;
    /* standard deviation of random perturbation to 1st dimension*/
    if (!sf_getfloat("vx1c",&vzc)) vzc = 0.0;
    /* 1st dimension v(z) chirp amplitude*/
    if (!sf_getfloat("x11c",&z1c)) z1c = fz;
    /* 1st dimension at which to begin chirp*/
    if (!sf_getfloat("x12c",&z2c)) z2c = fz+(nz-1)*dz;
    /* 1st dimension at which to end chirp*/
    if (!sf_getfloat("l1c",&l1c)) l1c = dz;
    /* wavelength at beginning of chirp*/
    if (!sf_getfloat("l2c",&l2c)) l2c = dz;
    /*wavelength at end of chirp*/
    if (!sf_getfloat("exc",&exc)) exc = 1.0;
    /* exponent of chirp*/

    sf_setformat(out,"native_float");
    sf_putint(out,"n1",nz);
    sf_putint(out,"n2",nx);
    sf_putint(out,"n3",ny);
    sf_putfloat(out,"d1",dz);
    sf_putfloat(out,"d2",dx);
    sf_putfloat(out,"d3",dy);
    sf_putfloat(out,"o1",fz);
    sf_putfloat(out,"o2",fx);
    sf_putfloat(out,"o3",fy);
	
    /* compute chirp constants */
    bc = SF_PI/(z2c-z1c)*(1.0/l2c-1.0/l1c);
    ac = 2.0*SF_PI/l1c - 2.0*bc*z1c;
	
    /* allocate space */
    v = sf_floatalloc(nz);
    vz = sf_floatalloc(nz);
	
    if (NULL != vzf) {
	sf_floatread(vz,nz,vzf);
	sf_fileclose(vzf);
    } else {
	for (iz=0; iz<nz; ++iz)
	    vz[iz] = 0.0;
    }

    /* random v(z) perturbation */
    /*for (iz=0; iz<nz; ++iz)
      vz[iz] += vzran*frannor();*/


    /* loop over x */
    for (ix=0,x=fx; ix<nx; ++ix,x+=dx) {
	
	/* loop over y */
	for (iy=0,y=fy; iy<ny; ++iy,y+=dy) {
		
	    /* compute top and bottom of lens */
/*
  ztop = atop+btop*(pow(x-xlens,2)+pow(y-ylens,2));
  zbot = abot+bbot*(pow(x-xlens,2)+pow(y-ylens,2));
*/
			
	    /* loop over z */
	    for (iz=0,z=fz; iz<nz; ++iz,z+=dz) {
				
		/* v(z) profile */
		v[iz] = vz[iz];
				
		/* constant + constant gradient */
		v[iz] += v000+x*dvdx+y*dvdy+z*dvdz;
				
		/* lens */
		xn = 2.0*(x-xlens)/dlens;
		ynn = 2.0*(y-ylens)/dlens;
		zn = 2.0*(z-zlens)/tlens;
		v[iz] += vlens*exp(-xn*xn-ynn*ynn-zn*zn);
		/*
		  if (z>zbot && z<ztop) v[iz] += vlens;
		*/

		/* chirp */
		if (z>z1c && z<z2c) {
		    vtemp = sinf((ac+bc*z)*z);
		    if (vtemp<0.0)
			v[iz] -= vzc*powf(-vtemp,exc);
		    else
			v[iz] += vzc*powf(vtemp,exc);
		}

		/* random perturbation */
		/*v[iz] += vran*frannor();*/
		/*warn("v00=%f dvdx=%f dvdz=%f vz=%f",v000,dvdx,dvdz,v[iz]);*/
	    }
			
	    /* write velocity function */
	    sf_floatwrite (v,nz,out);
	}
    }
	
    exit (0);
}
