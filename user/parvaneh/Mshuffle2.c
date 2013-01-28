/*Shuffling an array */
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


#include <rsf.h>
//#include <time.h>

static void bubble(int *in, int *idx, int n)
{
    int i1, i2, t1;
    for(i1=0; i1<n; i1++)
	idx[i1] = i1;
	
    for(i2=n; i2>0; i2--)
	for(i1=0; i1<i2; i1++)
	{
	    if(in[idx[i1]] > in[idx[i1+1]])
	    {
		t1 = idx[i1+1];
		idx[i1+1] = idx[i1];
		idx[i1] = t1;
	    }
	}
}


int main (int argc, char *argv[]) 
{    
    int m, k, l, seed;
    sf_file  bshuffle,ashuffle;
    sf_axis ax,at,ap,av;
    int nx,nt,np,nv,iteration,*a1, *a2;
    float ***bsh, ***bsh2;
    sf_init(argc,argv);

    bshuffle=sf_input("in");
    ashuffle=sf_output("out");

    if (!sf_getint("iteration",&iteration)) iteration=1;
    if (!sf_getint("seed",&seed)) seed=2012;
    at=sf_iaxa( bshuffle,1); nt=sf_n(at);
    ax=sf_iaxa( bshuffle,3); nx=sf_n(ax);
    ap=sf_iaxa( bshuffle,2); np=sf_n(ap);
    av=sf_iaxa( bshuffle,4); nv=sf_n(av);
    bsh=sf_floatalloc3(nt,np,nx);
    bsh2=sf_floatalloc3(nt,np,nx);
    a1=sf_intalloc(np);
    a2=sf_intalloc(np);
    sf_warning("ntpx=%d",nt*np*nx);

    srand(seed);


    for (m=0; m<np; m++) {
	a1[m]=rand();
    }

    bubble(a1, a2, np);


    for (k=0; k<nv; k++) {
	sf_floatread(bsh[0][0],nt*np*nx,bshuffle);
	for(l=0; l<nx; l++) {
	    for (m=0; m<np; m++) {
		memcpy(bsh2[l][m], bsh[l][a2[m]], nt*sizeof(float));
	    }
	}
	sf_floatwrite(bsh2[0][0],nt*np*nx,ashuffle);
    }

    free(bsh[0][0]);
    free(bsh[0]);
    free(bsh);
    free(bsh2[0][0]);
    free(bsh2[0]);
    free(bsh2);
    free(a1);
    free(a2);
    exit (0);
}
