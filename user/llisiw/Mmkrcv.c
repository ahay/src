/* Make receiver list for first-arrival traveltime tomography */
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
    int n[SF_MAX_DIM], nt, dim, i, j, k, is, nshot, nrecv, *mask, **recv;
    int offset1, offset2, left1, right1, left2, right2, count;
    float **source, d[SF_MAX_DIM], o[SF_MAX_DIM], temp, **t0, **reco;
    char key[4];
    sf_file in, topo, table, out, record;
    
    sf_init(argc,argv);
    in     = sf_input("in");
    out    = sf_output("out");
    record = sf_output("record");
    
    /* read input shot file */
    if (!sf_histint(in,"n2",&nshot)) nshot=1;

    source = sf_floatalloc2(3,nshot);
    sf_floatread(source[0],3*nshot,in);

    /* read model dimension from topography file */
    if (NULL == sf_getstring("topo"))
	sf_error("Need topography file topo=");
    topo = sf_input("topo");
    
    dim = sf_filedims(topo,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(topo,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(topo,key,o+i)) o[i]=0.;
	nt *= n[i];
    }
    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1];
    }

    mask = sf_intalloc(nt);
    sf_intread(mask,nt,topo);
    sf_fileclose(topo);

    /* read complete time table */
    if (NULL == sf_getstring("table"))
	sf_error("Need complete time file table=");
    table = sf_input("table");

    t0 = sf_floatalloc2(nt,nshot);
    sf_floatread(t0[0],nt*nshot,table);
    sf_fileclose(table);

    if (!sf_getint("offset1",&offset1)) offset1=0;
    /* receiver offset inline (on each side) */

    if (!sf_getint("offset2",&offset2)) offset2=0;
    /* receiver offset crossline (on each side) */

    nrecv = ((2*offset1+1)*(2*offset2+1)>n[1]*n[2])?n[1]*n[2]:(2*offset1+1)*(2*offset2+1);

    /* allocate memory for output */
    recv = sf_intalloc2(nrecv,nshot);
    reco = sf_floatalloc2(nrecv,nshot);

    /* write output header */
    sf_putint(out,"n1",nrecv);
    sf_putint(record,"n1",nrecv);
    sf_settype(out,SF_INT);

    for (is=0; is < nshot; is++) {
	temp   = (source[is][1]-o[1])/d[1];
	left1  = (int)temp-offset1;
	right1 = (int)temp+offset1;

	temp   = (source[is][2]-o[2])/d[2];
	left2  = (int)temp-offset2;
	right2 = (int)temp+offset2;

	count = 0;
	for (k=0; k < n[2]; k++) {
		for (j=0; j < n[1]; j++) {
		    if (left1 <= j && j <= right1 && left2 <= k && k <= right2) {
			for (i=0; i < n[0]; i++) {
			    if (mask[k*n[1]*n[0]+j*n[0]+i] == 1) {
				recv[is][count] = k*n[1]*n[0]+j*n[0]+i;
				reco[is][count] = t0[is][k*n[1]*n[0]+j*n[0]+i];
				count++;
			    }
			}
		    }
		}	    
	}
	
	/* negative flag for void space */
	if (count < nrecv) {
	    for (i=count; i < nrecv; i++) {
		recv[is][i] = -1;
		reco[is][i] = -1.0;
	    }
	}
    }
    
    sf_intwrite(recv[0],nrecv*nshot,out);
    sf_floatwrite(reco[0],nrecv*nshot,record);
    
    exit(0);
}
