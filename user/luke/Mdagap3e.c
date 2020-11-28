/* Reflection event apex protector/removal for dip-angle gathers.

   May be used for migration aperture optimization or for reflected energy
   supression. For the last multiply the output on -1.

   Input:
   dagFile.rsf - input dip-angle gathers;
   dipFile.rsf - dips esitimated in the image domain. The dips are in degree (!)

   Output:
   taperFile.rsf - mask for input dip-angle gathers
*/

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

int main (int argc, char* argv[]) 
{
// dip-angle gathers dimensions
    int zNum_;   float zStart_;   float zStep_;
    int dipNum1_; float dipStart1_; float dipStep1_;
    int dipNum2_; float dipStart2_; float dipStep2_;
    int xNum_;	 float xStart_;   float xStep_;
    int yNum_;   float yStart_;   float yStep_;
    sf_file dagFile, dipFile=NULL, taperFile;
    bool isDepthDep = true;
    float pwidth1 = 10.f;   float pwidth2 = 10.f;
    float greyarea1 = 0.f;   float greyarea2 = 0.f;
    float greyarea1a = 0.f;   float greyarea2a = 0.f;
    float dz = 0.f;
    int panelSize, taperSize, ix, iz, iy, ind, idI, idX, indt, cube;
    float dipI, dipX, curZ, curDip1, curDip2, dipRight1, dipRight2;
    float R1, R2, L1, L2, one, secondLeft1, secondRight1, secondLeft2, secondRight2, taper;
    float fixd, X, Y, radius, radiusC, slope;
    float y0, y1, x0, x1, e0, e1, esqr1, esqr0, ey0, ey1, t, t0, t1, r0, r1, f, d0, d1, distance, e0y0, x0de0, x0de0sqr, denom0;
    float y0a, y1a, inner_distance, outer_distance;
    const int imax = 2 * FLT_MAX_EXP;
    int i;

// Initialize RSF 
    sf_init (argc,argv);
// Input files
    dagFile = sf_input ("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dagFile) ) sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

    if ( NULL != sf_getstring("dips") ) {
	/* dips esitimated in the image domain (in degree) */ 
	dipFile  = sf_input ("dips");
    }

// Output file
    taperFile = sf_output("out");

//(Depth,InlineDip,CrosslineDip,Xaxis,Yaxis)
// Depth/time axis 
    if ( !sf_histint   (dagFile, "n1", &zNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dagFile, "d1", &zStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dagFile, "o1", &zStart_) ) sf_error ("Need o1= in input");
// Dip angle axis Inline
    if ( !sf_histint   (dagFile, "n2", &dipNum1_) )     sf_error ("Need n2= in input");
    if ( !sf_histfloat (dagFile, "d2", &dipStep1_) )    sf_error ("Need d2= in input");
    if ( !sf_histfloat (dagFile, "o2", &dipStart1_) )   sf_error ("Need o2= in input");
// Dip angle axis Crossline
    if ( !sf_histint   (dagFile, "n3", &dipNum2_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (dagFile, "d3", &dipStep2_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (dagFile, "o3", &dipStart2_) )   sf_error ("Need o3= in input");
// x axis 
    if ( !sf_histint   (dagFile, "n4", &xNum_) )   sf_error ("Need n4= in input");
    if ( !sf_histfloat (dagFile, "d4", &xStep_) )  sf_error ("Need d4= in input");
    if ( !sf_histfloat (dagFile, "o4", &xStart_) ) sf_error ("Need o4= in input");
// y axis 
    if ( !sf_histint   (dagFile, "n5", &yNum_) )   sf_error ("Need n5= in input");
    if ( !sf_histfloat (dagFile, "d5", &yStep_) )  sf_error ("Need d5= in input");
    if ( !sf_histfloat (dagFile, "o5", &yStart_) ) sf_error ("Need o5= in input");
// tapering paremeters
    if (!sf_getbool ("ddep", &isDepthDep)) isDepthDep = true;
    /* if y, taper depends on depth; if n, no */

    if ( !sf_getfloat ("pwidth1", &pwidth1) ) pwidth1 = 10.f;
    if ( !sf_getfloat ("pwidth2", &pwidth2) ) pwidth2 = 10.f;
    /* protected width (in degree) */

    if ( !sf_getfloat ("greyarea1", &greyarea1) ) greyarea1 = 10.f;
    if ( !sf_getfloat ("greyarea2", &greyarea2) ) greyarea2 = 10.f;
    /* width of event tail taper (in degree) */
    if ( greyarea1 < 0 ) {sf_warning ("inline greyarea value is changed to 10"); greyarea1 = 10.f;}
    if ( greyarea2 < 0 ) {sf_warning ("crossline greyarea value is changed to 10"); greyarea2 = 10.f;}
    if ( !sf_getfloat ("dz", &dz) ) dz = 20.f;
    /* half of a migrated wave length */
    if (dz < 0) {sf_warning ("dz value is changed to 20"); dz = 20.f;}

    if (greyarea1 > pwidth1) {
        greyarea1 = pwidth1;
        sf_warning ("greyarea1 changed to %g", pwidth1);
       }
    if (greyarea2 > pwidth2) {
        greyarea2 = pwidth2;
        sf_warning ("greyarea1 changed to %g", pwidth2);
       }

    // input dips
    panelSize = xNum_ * yNum_ * zNum_ * 2;
    float* dipPanel   = sf_floatalloc (panelSize);
    sf_seek (dipFile, 0, SEEK_SET);		
    sf_floatread (dipPanel, panelSize, dipFile);
    // output taper
    taperSize = xNum_ * yNum_* dipNum1_ * dipNum2_ * zNum_;
    float* taperPanel = sf_floatalloc (taperSize);
    memset (taperPanel, 0, taperSize * sizeof (float));

//define cubesize for shift
    cube = xNum_ * yNum_ * zNum_;

for (iy = 0; iy < yNum_; ++iy){
    for (ix = 0; ix < xNum_; ++ix){

	for (iz = 0; iz < zNum_; ++iz) {
	    ind = iy * xNum_ * zNum_ + ix * zNum_ + iz; 
	    dipI = dipPanel [ind]; // local inline slope in degree
            dipX = dipPanel [ind + cube];     // local crossline slope in degrees

	    curZ = zStart_ + iz * zStep_;
	    if (! (curZ - dz > 0) ) continue; // out from data



            for (idX = 0; idX < dipNum2_; ++idX){
                curDip2 = dipStart2_ + idX * dipStep2_;
	     for (idI = 0; idI < dipNum1_; ++idI) {

		curDip1 = dipStart1_ + idI * dipStep1_;
                //define index
		indt = iy *xNum_ * dipNum2_ * dipNum1_ * zNum_  + ix * dipNum2_ * dipNum1_ * zNum_ + idX * dipNum1_ * zNum_ + idI * zNum_ + iz;
//		int indDip = (idX * dipNum1_ + idI) * zNum_ + iz;


		dipRight1 = 0.f;
		dipRight2 = 0.f;

		if (isDepthDep) { // depth-dependent taper 
                    one = iz*1.;
                    fixd = (one/zNum_+0.5)*(one/zNum_+0.5);
                    dipRight1=pwidth1+pwidth1/fixd;
                    dipRight2=pwidth2+pwidth2/fixd;

                    greyarea1a=greyarea1+greyarea1/fixd;
                    greyarea2a=greyarea2+greyarea2/fixd;
		} else { // taper width is constant along depths

               
                    dipRight1=pwidth1;
                    dipRight2=pwidth2;
                    greyarea1a=greyarea1;
                    greyarea2a=greyarea2;
                       }

		    secondLeft1  = dipRight1  - greyarea1a;
		    secondRight1 = dipRight1 + greyarea1a;


		    secondLeft2  = dipRight2  - greyarea2a;
		    secondRight2 = dipRight2 + greyarea2a;
		taper = 1.;

                R1 = (curDip1 - dipI)/secondRight1;
                R2 = (curDip2 - dipX)/secondRight2;
               if ( R1*R1 + R2*R2 > 1) {
                  taper = 0.;
                    } //Out of Elipse
               else {
                     L1 = (curDip1 - dipI)/secondLeft1;
                     L2 = (curDip2 - dipX)/secondLeft2;
                           

                     if ( L1*L1 + L2*L2 > 1) {

				//Position Coordinate
                           X = curDip1 - dipI;
                           Y = curDip2 - dipX;


				//Has the User Defined a circle?
				if ((greyarea1 == greyarea2) && (pwidth1 == pwidth2)) {
                           		radius = sqrt(X*X + Y*Y);
                           		radiusC = radius - secondLeft1;
                           
                           		slope = 1./(secondRight1 - secondLeft1);
                           		taper = 1. - radiusC*slope;//Linear Fit between gray area elipses
				}else{
				//This solves for the distance to an inner elipse.  Only works in the first quadrant.
				//Should not be an issue because only concerned about distance, not the point on the elipse

				//initialize variables
					y0 = 0;
					y1 = 0;
					e0 = 0;
					d1 = 0;
				
				if ( X < 0.){ y0a = -X;
				}else{ y0a = X;}
				if ( Y < 0.){ y1a = -Y;
				}else{ y1a = Y;}
				//This alogrythm assumes y0 is the coordinate with the major axis
				//And y1 is the coordinate with the minor axis of the elipse
				if (secondLeft1 > secondLeft2){
					y0 = y0a;
					y1 = y1a;
					e0 = secondLeft1;
					e1 = secondLeft2;
				}else{
					y1 = y0a;
					y0 = y1a;
					e1 = secondLeft1;
					e0 = secondLeft2;
				}

				if (y1 > 0){
		
					if (y0 > 0){
					//Bisect to compute root for parameter t
						esqr0 = e0*e0;
                			        esqr1 = e1*e1;
  						ey0 = e0*y0;
  						ey1 = e1*y1;
                        
                   			        t0 = -esqr1 + ey1;
                     			        t1 = -esqr1 + sqrt(ey0*ey0+ey1*ey1);
                        			t = t0;
						for (i=0; i < imax; i++){
							t = (0.5)*(t0+t1);
                        				if (t == t0 || t == t1){
								break;
							}			
							r0 = ey0/(t + esqr0);
							r1 = ey1/(t + esqr1);

							f = r0*r0 + r1*r1 -1.;
							if (f > 0.){
								t0=t;
							}else if (f < 0.){
								t1=t;
							}else{
								break;
							}
						}
						x0 = esqr0*y0/(t+esqr0);
						x1 = esqr1*y1/(t+esqr1);
						d0 = x0 - y0;
						d1 = x1 - y1;

						distance = sqrt(d0*d0 + d1*d1);
					} else { //y0 == 0
						x0=0.;
						x1=e1;
						distance = fabs(y1-e1);
						}
				}else{ //y1==0
					denom0 = e0*e0 - e1*e1;
					e0y0 = e0*y0;
					if (e0y0 < denom0){
					// y0 is inside the subinterval
						x0de0 = e0y0/denom0;
						x0de0sqr = x0de0*x0de0;
						x0 = e0*x0de0;
						x1 = e1*sqrt(fabs(1. - x0de0sqr));
						d0 = x0 - y0;
						distance = sqrt(d0*d0 + x1*x1);
					}else{
					// y0 is outside the subinterval.  The closest ellipse point has
					//x1 ==0 and is on the boundary interval
						x0 = e0;
						x1 = 0.;
						distance = fabs(y0 - e0);
					}
				}
				
				inner_distance=distance;
					//repeat for outer elipse distance
				if (secondRight1 > secondRight2){
					y0 = y0a;
					y1 = y1a;
					e0 = secondRight1;
					e1 = secondRight2;
				}else{
					y1 = y0a;
					y0 = y1a;
					e1 = secondRight1;
					e0 = secondRight2;
				}

				if (y1 > 0){
		
					if (y0 > 0){
					//Bisect to compute root for parameter t
						esqr0 = e0*e0;
                			        esqr1 = e1*e1;
  						ey0 = e0*y0;
  						ey1 = e1*y1;
                        
                   			        t0 = -esqr1 + ey1;
                     			        t1 = -esqr1 + sqrt(ey0*ey0+ey1*ey1);
                        			t = t0;
						for (i=0; i < imax; i++){
							t = (0.5)*(t0+t1);
                        				if (t == t0 || t == t1){
								break;
							}			
							r0 = ey0/(t + esqr0);
							r1 = ey1/(t + esqr1);

							f = r0*r0 + r1*r1 -1.;
							if (f > 0.){
								t0=t;
							}else if (f < 0.){
								t1=t;
							}else{
								break;
							}
						}
						x0 = esqr0*y0/(t+esqr0);
						x1 = esqr1*y1/(t+esqr1);
						d0 = x0 - y0;
						d1 = x1 - y1;

						distance = sqrt(d0*d0 + d1*d1);
					} else { //y0 == 0
						x0=0.;
						x1=e1;
						distance = fabs(y1-e1);
					}
				}else{ //y1==0
					denom0 = e0*e0 - e1*e1;
					e0y0 = e0*y0;
					if (e0y0 < denom0){
					// y0 is inside the subinterval
						x0de0 = e0y0/denom0;
						x0de0sqr = x0de0*x0de0;
						x0 = e0*x0de0;
						x1 = e1*sqrt(fabs(1. - x0de0sqr));
						d0 = x0 - y0;
						distance = sqrt(d0*d0 + x1*x1);
					}else{
					// y0 is outside the subinterval.  The closest ellipse point has
					//x1 ==0 and is on the boundary interval
						x0 = e0;
						x1 = 0.;
						distance = fabs(y0 - e0);
					}
				}
			
				outer_distance=distance;
					
				slope = 1./(inner_distance + outer_distance);
                           	taper = 1. - inner_distance*slope;//Linear Fit between gray area elipses	

				}//end elipse			
                          }  //end grey area 
                    else {
                             taper = 1.;
                             } //In Preserved, Non Gray Area

                    } 		
		taperPanel[indt] = taper;			

	    } //closes idI
          }//closes idX
	}//closes idz
    }//closes dx
}//closes dy
    sf_floatwrite (taperPanel, taperSize, taperFile);
	

    sf_warning (".");

    sf_fileclose (dipFile);
    sf_fileclose (dagFile);
    sf_fileclose (taperFile);

    free (dipPanel);
    free (taperPanel);

    return 0;
}
