/* Double square-root eikonal solver (2D + explicit) */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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
    bool velocity;
    int dim, i, ss[3], n[SF_MAX_DIM], nm, nt, j, k;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *s, *t, upd[2], temp[2];
    char key[6];
    sf_file in, out;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    /* read input dimension */
    dim = sf_filedims(in,n);

    nm = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input.",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nm *= n[i];
    }
    if (dim < 3) {
	/* extend the third dimension for output (copy second dimension) */
	n[2] = n[1]; d[2] = d[1]; o[2] = o[1];
    }

    /* read input */
    s = sf_floatalloc(nm);
    sf_floatread(s,nm,in);

    /* write output dimension */
    sf_putint(out,"n3",n[2]);
    sf_putfloat(out,"d3",d[2]);
    sf_putfloat(out,"o3",o[2]);

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */

    /* convert to slowness squared */
    if (velocity) {
	for (i=0; i < nm; i++)
	    s[i] = 1./s[i]*1./s[i];
    }

    nt = nm*n[2];
    ss[0] = 1; ss[1] = n[0]; ss[2] = n[0]*n[1];

    /* allocate memory for output */
    t = sf_floatalloc(nt);

    /* initialize */
    for (k=0; k < n[0]; k++) {
	for (i=0; i < n[1]; i++) {
	    t[i*ss[2]+i*ss[1]+k] = 0.;
	}
    }

    for (j=1; j < n[1]; j++) {
	for (i=0; i < n[1]-j; i++) {
	    t[i*ss[2]+(i+j)*ss[1]+n[0]-1] = 
		SF_MIN(t[i*ss[2]+(i+j-1)*ss[1]+n[0]-1]+sqrtf(s[(i+j)*ss[1]+n[0]-1])*d[1],
		       t[(i+1)*ss[2]+(i+j)*ss[1]+n[0]-1]+sqrtf(s[i*ss[1]+n[0]-1])*d[2]);
	}
    }

    /* sweep in -z direction */
    for (k=n[0]-1; k > 0; k--) {
	for (j=1; j < n[1]; j++) {
	    for (i=0; i < n[1]-j; i++) {
		/* singular = true; */

		upd[0] = SF_MAX(0.,(t[i*ss[2]+(i+j)*ss[1]+k]-t[i*ss[2]+(i+j-1)*ss[1]+k])/d[1]);
		if (s[(i+j)*ss[1]+k]-powf(upd[0],2.) >= 0.) {
		    temp[0] = sqrtf(s[(i+j)*ss[1]+k]-powf(upd[0],2.));
		    /* singular = false; */
		} else {
		    temp[0] = sqrtf(s[(i+j)*ss[1]+k]);
		}

		upd[1] = SF_MAX(0.,(t[i*ss[2]+(i+j)*ss[1]+k]-t[(i+1)*ss[2]+(i+j)*ss[1]+k])/d[2]);
		if (s[i*ss[1]+k]-powf(upd[1],2.) >= 0.) {
		    temp[1] = sqrtf(s[i*ss[1]+k]-powf(upd[1],2.));
		    /* singular = false; */
		} else {
		    temp[1] = sqrtf(s[i*ss[1]+k]);
		}
		
		t[i*ss[2]+(i+j)*ss[1]+k-1] = 
		    SF_MIN(SF_MIN(t[i*ss[2]+(i+j-1)*ss[1]+k-1]+sqrtf(s[(i+j)*ss[1]+k-1])*d[1],
				  t[(i+1)*ss[2]+(i+j)*ss[1]+k-1]+sqrtf(s[i*ss[1]+k-1])*d[2]),
			   t[i*ss[2]+(i+j)*ss[1]+k]+d[0]*(temp[0]+temp[1]));
	    }
	}
    }

    /* mirror */
    for (k=1; k < n[1]; k++) {
	for (j=0; j < k; j++) {
	    for (i=0; i < n[0]; i++) {
		t[k*ss[2]+j*ss[1]+i] = t[j*ss[2]+k*ss[1]+i];
	    }
	}
    }

    /* write output */
    sf_floatwrite(t,nt,out);

    exit(0);
}
