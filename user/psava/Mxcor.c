/* Cross-correlation of time series (zero-lag output) */
/*
  Copyright (C) 2006 Colorado School of Mines
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;

    /* I/O files */
    sf_file Fi,Fs,Fr;

    /* cube axes */
    sf_axis ac,a1,a2,aa;
    int     nc,nn;
    int     in;
    int     nbuf,ibuf;

    /* arrays */
    float *ii=NULL, **us=NULL,**ur=NULL;

    int ompchunk; 
    int axis;

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    if(! sf_getint("axis",&axis)) axis=2;              /* cross-correlation axis */
    if(! sf_getint("nbuf",&nbuf)) nbuf=1;              /* buffer size */

    Fs = sf_input ("in" );
    Fr = sf_input ("uu" );
    Fi = sf_output("out");

    aa=sf_maxa(1,0,1); 

    /* read axes */
    if (axis==3) {
	a1=sf_iaxa(Fs,1); sf_setlabel(a1,"a1"); if(verb) sf_raxa(a1);
	a2=sf_iaxa(Fs,2); sf_setlabel(a2,"a2"); if(verb) sf_raxa(a2);
	ac=sf_iaxa(Fs,3); sf_setlabel(ac,"ac"); if(verb) sf_raxa(ac);
	sf_oaxa(Fi,aa,3);

	nn = sf_n(a1)*sf_n(a2);
    } else {
	a1=sf_iaxa(Fs,1); sf_setlabel(a1,"a1"); if(verb) sf_raxa(a1);
	ac=sf_iaxa(Fs,2); sf_setlabel(ac,"ac"); if(verb) sf_raxa(ac);
	sf_oaxa(Fi,aa,2);

	nn = sf_n(a1);
    }
    nc = sf_n(ac);

    nbuf = SF_MIN(nbuf,nc);

    /* allocate work arrays */
    us=sf_floatalloc2(nn,nbuf);
    ur=sf_floatalloc2(nn,nbuf);
    ii=sf_floatalloc (nn);

    /* init output */
    for (in=0; in<nn; in++) {
	ii[in]=0;
    }
    
    for (; nc > 0; nc -= nbuf) {
	if(verb) sf_warning("nsiz=%ld nbuf=%ld",nc,nbuf);
	if (nbuf > nc) nbuf=nc;

	sf_floatread(us[0],nn*nbuf,Fs);
	sf_floatread(ur[0],nn*nbuf,Fr);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(ibuf,in) shared(nn,ii,us,ur)
#endif
	for(ibuf=0; ibuf<nbuf; ibuf++) {
	    for(in=0; in<nn; in++) {
		ii[in] += us[ibuf][in] * ur[ibuf][in];
	    }
	}
    }

    sf_floatwrite(ii,nn,Fi);

    exit (0);
}
