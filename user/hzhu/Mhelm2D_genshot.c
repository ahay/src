/* Generate shot file for Helmholtz solver. */
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

int main(int argc, char* argv[])
{
	sf_file out, fmag;
	int n1,n2,ns,nw,i,j,k,iw;
	int srcz, srcx0, srcdx;
	int isource, nsource, dsource;
	float d1,d2,ds,dw,ow;
	float mag, *mags;
	float ***f;

	sf_init(argc, argv);
	out = sf_output("out");

	if (!sf_getint("n1",&n1)) n1=1;
	if (!sf_getint("n2",&n2)) n2=1;
	if (!sf_getint("ns",&ns)) ns=1;
	if (!sf_getfloat("d1",&d1)) d1=0.1;
    if (!sf_getfloat("d2",&d2)) d2=0.1;
        
    if (!sf_getint("nw",&nw)) nw=1;
    if (!sf_getfloat("dw",&dw)) dw=1.0;
    if (!sf_getfloat("ow",&ow)) ow=1.0;
    
    if (!sf_getint("nsource",&nsource)) nsource=1;
    if (!sf_getint("dsource", &dsource)) dsource=0;
    if (!sf_getint("srcz",&srcz)) sf_error("No srcz=.");
    if (!sf_getint("srcx0",&srcx0)) sf_error("No srcx0=.");
    if (!sf_getint("srcdx",&srcdx)) sf_error("No srcdx=.");

	if(NULL != sf_getstring("fmag")){
		fmag=sf_input("fmag");
		mags=sf_floatalloc(nw);
		sf_floatread(mags, nw, fmag);
	}else{
		fmag=NULL;
		mags=NULL;
		if (!sf_getfloat("mag",&mag)) mag=1.0;
	}

    f=sf_floatalloc3(n1,n2,ns);
	ds=d2*srcdx;

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putint(out,"n3",ns);
    sf_putint(out,"n4",nw);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"d3",ds);
    sf_putfloat(out,"d4",dw);
	sf_putfloat(out,"o4",ow);

	for (iw=0; iw< nw; iw++ ) {
		if(fmag!=NULL) mag=mags[iw];
		for (k=0; k< ns; k++) { 
			for (j=0; j<n2; j++) { 
				for (i=0; i<n1; i++) {
					f[k][j][i]=0.0;
				} /* for i */
			} /* for j */

			for(j=0; j<n2; j++){
				for(i=0; i<n1; i++){
					if (i == srcz ) {
						if ( j == srcx0+k*srcdx) {
							for(isource=0; isource<nsource; isource++){
								f[k][j+isource*dsource][i]=mag;
							}
						}
					}
				} /* for i */
			} /* for j */
		} /* for k */

		sf_floatwrite(f[0][0],n1*n2*ns,out);

	} /* for iw */
	exit(0);
}
