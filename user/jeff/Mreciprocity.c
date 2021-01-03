/* Create a dipping layer model for HTI testing purposes.  Has fixed velocity structure, but can change dip of layer and degree of anisotropy.*/

/*
  Copyrite (C) 2012 University of Western Australia

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public Licences as published by
  the Free Software Foundation; either version 2 of the License, or 
  (at you option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License 
  along with this program; if not, write to the Free Software Foundation
  Inc., 59 Temle Place, Suite 330, Boston, MA 02111-1307, USA.
*/

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{

    sf_file Fin=NULL; /* Input model size */
    sf_file Fout=NULL; /* Stiffness model */

    float ***in=NULL;
    float ***out=NULL;

    int nt,nr,ns,nrout;
    float dr,ds,or,maxr;

    float rr;
    int rloc, /* sloc, */ isr,irr;

    sf_axis at,ar,as,aout; /* Cube axes */
    int is, ir, it;
  
    /*----------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
	
    Fin  = sf_input ("in") ; /* Input field */
    Fout = sf_output("out"); /* Output field */

    at=sf_iaxa(Fin,1); sf_setlabel(at,"at"); sf_raxa(at);
    ar=sf_iaxa(Fin,2); sf_setlabel(ar,"ar"); sf_raxa(ar);
    as=sf_iaxa(Fin,3); sf_setlabel(as,"as"); sf_raxa(as);

    aout=ar;

    nt = sf_n(at); 
    nr = sf_n(ar); dr = sf_d(ar); or = sf_o(ar);
    ns = sf_n(as); ds = sf_d(as); /* os = sf_o(as); */

    maxr=-((float)(nr-1)*dr+or);
    nrout = 2*(int)(-maxr/dr)+1;

    /* Allocate input/output */
    in  = sf_floatalloc3(nt,nr   ,ns);
    out = sf_floatalloc3(nt,nrout,ns);
   
    sf_setn(aout,nrout);
    sf_seto(aout,maxr);
    sf_setd(aout,dr);	
    sf_setlabel(aout,"aout"); sf_raxa(aout);

    sf_oaxa(Fout,at  ,1);
    sf_oaxa(Fout,aout,2);
    sf_oaxa(Fout,as  ,3);

    /* Read in Data */
    sf_floatread( in[0][0],ns*nr*nt,Fin );
	
    for (is=0; is < ns; is++) {
	/* sloc = (int)((os+is*ds)/dr); */
		
	for (ir=0; ir < nr; ir++) {

	    /* receiver location in split-spread */
	    rr = ((float)(ir-1)*dr+or-maxr);
			
	    /* receiver integer location */
	    rloc = (int)(rr/dr)+1;
			
	    /* Compute location of true trace */
	    for (it=0; it < nt; it++) {
		out[is][rloc][it]=in[is][ir][it];
	    }
			
	    /* Reciprocal source and receiver*/
	    isr=(int)(ir*dr/ds)+(int)(or/ds)+is;
	    irr=(nrout-1)/2-ir-(int)(or/dr);
			
	    if (isr > 0 && isr < ns){
		for (it=0; it < nt; it++) {
		    out[isr][irr][it]=in[is][ir][it];
		}			
	    }
			
	    /* Compute reciprocal trace */
			
			
	}
    }

    /* Output section */
    sf_floatwrite(out[0][0],nt*nrout*ns,Fout);

    /*----------------------------------------*/
    /* deallocate */
    free(**in); free(*in); free(in);
    free(**out);free(*out);free(out);

    exit(0);
}
