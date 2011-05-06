/* Make receiver mask for Mfatomo */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
    int n[SF_MAX_DIM], nt, dim, i, j, is, nshot, offset, *topo, **mask, left, right;
    float **source, d[SF_MAX_DIM], o[SF_MAX_DIM], temp;
    char key[4];
    sf_file in, shot, out;
    
    sf_init(argc,argv);
    in   = sf_input("in");
    out  = sf_output("out");
    
    /* read input dimension */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i];
    }

    topo = sf_intalloc(nt);
    sf_intread(topo,nt,in);

    /* read in shot file */
    if (NULL == sf_getstring("shot"))
	sf_error("Need source shot=");
    shot = sf_input("shot");

    if (!sf_histint(shot,"n2",&nshot)) nshot=1;

    source = sf_floatalloc2(3,nshot);
    sf_floatread(source[0],3*nshot,shot);
    sf_fileclose(shot);

    if (!sf_getint("offset",&offset)) offset=1;
    /* receiver offset (on each side) */

    /* allocate memory for output */
    mask = sf_intalloc2(nt,nshot);

    sf_putint(out,"n4",nshot);

    /* modify for 3D */
    for (is=0; is < nshot; is++) {
	temp = (source[is][1]-o[1])/d[1];

	left  = (int)temp-offset;
	right = (int)temp+offset;
	
	for (j=0; j < n[1]; j++) {
	    
	    if (j < left || j > right) {
		
		for (i=0; i < n[0]; i++) {
		    mask[is][j*n[0]+i] = 0;
		}

	    } else {

		for (i=0; i < n[0]; i++) {
		    mask[is][j*n[0]+i] = topo[j*n[0]+i];
		}
		
	    }
	}
        
    }
    
    sf_intwrite(mask[0],nt*nshot,out);
    
    exit(0);
}
