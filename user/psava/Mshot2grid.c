/* Synthesize shot/receiver wavefields for 3-D SR migration */
/*
  Copyright (C) 2007 Colorado School of Mines
  
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

int main (int argc, char *argv[])
{
    sf_axis aw,at;
    sf_axis arx,asx,ax;
    sf_axis ary,asy,ay;
    sf_axis ae;

    sf_file Fs;  /* source   file sou[nw]         */
    sf_file Fr;  /* receiver file rec[nw][nrx][nry][nsx][nsy] */
    sf_file Fsou; /* source wavefield [nw][ nx][ ny][ne] */
    sf_file Frec; /* source wavefield [nw][ nx][ ny][ne] */

    sf_complex ***c_irec;
    sf_complex   *c_isou;
    sf_complex ***c_osou;
    sf_complex ***c_orec;

    float      ***r_irec;
    float        *r_isou;
    float      ***r_osou;
    float      ***r_orec;

    int   isx,isy,irx,iry,ix,iy,iw,it;
    int   nx,ny;
    float ox,dx,oy,dy;
    float  sx, sy, rx, ry,xx,yy;

    bool cflag;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    Fr  = sf_input ( "in");
    Fs  = sf_input ("wav");
    Fsou= sf_output("swf"); 
    Frec= sf_output("out");

    cflag= (SF_COMPLEX == sf_gettype(Fs))?true:false;
	
    if(cflag) {
	sf_settype(Fsou,SF_COMPLEX);
	sf_settype(Frec,SF_COMPLEX);
    }
    
    arx = sf_iaxa(Fr,2); sf_setlabel(arx,"rx");
    ary = sf_iaxa(Fr,3); sf_setlabel(ary,"ry");
    asx = sf_iaxa(Fr,4); sf_setlabel(asx,"sx");
    asy = sf_iaxa(Fr,5); sf_setlabel(asy,"sy");

    ae = sf_maxa(sf_n(asx)*sf_n(asy),0,1); sf_setlabel(ae,"e");  /* experiments axis */

    if(cflag) {
	aw  = sf_iaxa(Fr,1); sf_setlabel(aw,"w");
	sf_oaxa(Fsou,aw,1); 
	sf_oaxa(Frec,aw,1);
	at = NULL;
    } else {
	at  = sf_iaxa(Fr,1); sf_setlabel(at,"t");
	sf_oaxa(Fsou,at,1); 
	sf_oaxa(Frec,at,1);
	aw = NULL;
    }

    if (!sf_getint  ("nx",&nx)) sf_error ("Need nx="); /* x samples */
    if (!sf_getfloat("dx",&dx)) sf_error ("Need dx="); /* x sampling */
    if (!sf_getfloat("ox",&ox)) sf_error ("Need ox="); /* x origin */
    ax = sf_maxa(nx,ox,dx); sf_setlabel(ax,"x");

    if (!sf_getint  ("ny",&ny)) ny=1; /* y samples */
    if (!sf_getfloat("dy",&dy)) dy=1; /* y sampling */
    if (!sf_getfloat("oy",&oy)) oy=0; /* y origin */
    ay = sf_maxa(ny,oy,dy); sf_setlabel(ay,"y");

    sf_oaxa(Fsou, ax,2); sf_oaxa(Frec, ax,2);
    sf_oaxa(Fsou, ay,3); sf_oaxa(Frec, ay,3);
    sf_oaxa(Fsou, ae,4); sf_oaxa(Frec, ae,4);

    /*------------------------------------------------------------*/
    if(cflag) {
	c_isou = sf_complexalloc (sf_n(aw));
	c_irec = sf_complexalloc3(sf_n(aw),sf_n(arx),sf_n(ary));
	c_osou = sf_complexalloc3(sf_n(aw),sf_n(ax),sf_n(ay));
	c_orec = sf_complexalloc3(sf_n(aw),sf_n(ax),sf_n(ay));
	r_isou = NULL;
	r_irec = NULL;
	r_osou = NULL;
	r_orec = NULL;
    } else {
	r_isou = sf_floatalloc (sf_n(at));
	r_irec = sf_floatalloc3(sf_n(at),sf_n(arx),sf_n(ary));
	r_osou = sf_floatalloc3(sf_n(at),sf_n(ax),sf_n(ay));
	r_orec = sf_floatalloc3(sf_n(at),sf_n(ax),sf_n(ay));
    }

    /*------------------------------------------------------------*/
    /* SOURCE wavefield */
    /*------------------------------------------------------------*/
    if(cflag) {
	sf_complexread(c_isou,sf_n(aw),Fs);

	for(    isy=0;isy<sf_n(asy);isy++) { sy = sf_o(asy) + isy * sf_d(asy);
	    for(isx=0;isx<sf_n(asx);isx++) { sx = sf_o(asx) + isx * sf_d(asx);
		
		for        (iy=0;iy<sf_n(ay);iy++) {
		    for    (ix=0;ix<sf_n(ax);ix++) {
			for(iw=0;iw<sf_n(aw);iw++) {
			    c_osou[iy][ix][iw] = sf_cmplx(0.,0.);
			}
		    }
		}
		
		yy = sy; iy = 0.5+(yy-oy)/dy;
		xx = sx; ix = 0.5+(xx-ox)/dx;
		if(iy>=0 && iy<sf_n(ay) && ix>=0 && ix<sf_n(ax)) {
		    
		    for( iw=0;iw<sf_n(aw);iw++) {
			c_osou[iy][ix][iw] = c_isou[iw];
		    }
		}
		
		sf_complexwrite(c_osou[0][0],sf_n(ax)*sf_n(ay)*sf_n(aw),Fsou);	
	    } /* isx */
	}     /* isy */	

    } else {
	sf_floatread  (r_isou,sf_n(at),Fs);

	for(    isy=0;isy<sf_n(asy);isy++) { sy = sf_o(asy) + isy * sf_d(asy);
	    for(isx=0;isx<sf_n(asx);isx++) { sx = sf_o(asx) + isx * sf_d(asx);
		
		for        (iy=0;iy<sf_n(ay);iy++) {
		    for    (ix=0;ix<sf_n(ax);ix++) {
			for(it=0;it<sf_n(at);it++) {
			    r_osou[iy][ix][it] = 0.0;
			}
		    }
		}
		
		yy = sy; iy = 0.5+(yy-oy)/dy;
		xx = sx; ix = 0.5+(xx-ox)/dx;
		if(iy>=0 && iy<sf_n(ay) && ix>=0 && ix<sf_n(ax)) {
		    
		    for( it=0;it<sf_n(at);it++) {
			r_osou[iy][ix][it] = r_isou[it];
		    }   
		}
		
		sf_floatwrite  (r_osou[0][0],sf_n(ax)*sf_n(ay)*sf_n(at),Fsou);	
	    } /* isx */
	}     /* isy */
    }
    
    /*------------------------------------------------------------*/
    /* RECEIVER wavefield */
    /*------------------------------------------------------------*/

    if(cflag) {
	
	for(    isy=0;isy<sf_n(asy);isy++) { sy = sf_o(asy) + isy * sf_d(asy);
	    for(isx=0;isx<sf_n(asx);isx++) { sx = sf_o(asx) + isx * sf_d(asx);
				
		for        (iy=0;iy<sf_n(ay);iy++) {
		    for    (ix=0;ix<sf_n(ax);ix++) {
			for(iw=0;iw<sf_n(aw);iw++) {
			    c_orec[iy][ix][iw] = sf_cmplx(0.,0.);
			}
		    }   
		}

		sf_complexread(c_irec[0][0],sf_n(arx)*sf_n(ary)*sf_n(aw),Fr);

		for(    iry=0;iry<sf_n(ary);iry++) { ry = sf_o(ary) + iry * sf_d(ary);
		    for(irx=0;irx<sf_n(arx);irx++) { rx = sf_o(arx) + irx * sf_d(arx);
			
			yy = sy + ry; iy = 0.5+(yy-oy)/dy;
			xx = sx + rx; ix = 0.5+(xx-ox)/dx;
			
			if(iy>=0 && iy<sf_n(ay) && ix>=0 && ix<sf_n(ax)) {
			    for( iw=0;iw<sf_n(aw);iw++) {
				c_orec[iy][ix][iw] = c_irec[iry][irx][iw];
			    }
			}
		    } /* irx */
		}     /* iry */
		
		sf_complexwrite(c_orec[0][0],sf_n(ax)*sf_n(ay)*sf_n(aw),Frec);		
	    } /* isx */
	}     /* isy */

    } else {

	for(    isy=0;isy<sf_n(asy);isy++) { sy = sf_o(asy) + isy * sf_d(asy);
	    for(isx=0;isx<sf_n(asx);isx++) { sx = sf_o(asx) + isx * sf_d(asx);
		
		for        (iy=0;iy<sf_n(ay);iy++) {
		    for    (ix=0;ix<sf_n(ax);ix++) {
			for(it=0;it<sf_n(at);it++) {
			    r_orec[iy][ix][it] = 0.0;
			}
		    }
		}
		
		sf_floatread  (r_irec[0][0],sf_n(arx)*sf_n(ary)*sf_n(at),Fr);
		
		for(    iry=0;iry<sf_n(ary);iry++) { ry = sf_o(ary) + iry * sf_d(ary);
		    for(irx=0;irx<sf_n(arx);irx++) { rx = sf_o(arx) + irx * sf_d(arx);
			
			yy = sy + ry; iy = 0.5+(yy-oy)/dy;
			xx = sx + rx; ix = 0.5+(xx-ox)/dx;
			
			if(iy>=0 && iy<sf_n(ay) && ix>=0 && ix<sf_n(ax)) {  
			    for( it=0;it<sf_n(at);it++) {
				r_orec[iy][ix][it] = r_irec[iry][irx][it];
			    }
			}
		    } /* irx */
		}     /* iry */
		
		sf_floatwrite  (r_orec[0][0],sf_n(ax)*sf_n(ay)*sf_n(at),Frec);
	    } /* isx */
	}     /* isy */
    }
    
    /*------------------------------------------------------------*/    
    if(cflag) {
	free(c_isou);
	free(**c_irec); free(*c_irec); free(c_irec);
	free(**c_osou); free(*c_osou); free(c_osou);
	free(**c_orec); free(*c_orec); free(c_orec);
    } else {
	free(r_isou);
	free(**r_irec); free(*r_irec); free(r_irec);
	free(**r_osou); free(*r_osou); free(r_osou);
	free(**r_orec); free(*r_orec); free(r_orec);
    }
    /*------------------------------------------------------------*/
    

    exit (0);
}

