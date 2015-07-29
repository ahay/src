/* 3-D modeling/migration with extended SSF */

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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "wex.h"
#include "wextap.h"
#include "wexssr.h"
#include "wexslo.h"
#include "wexutl.h"
/*^*/

int main (int argc, char *argv[])
{
    bool verb;            /* verbosity */
    bool causal;          /* causal/anti-causal flag */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    int wexsign;
    int ompnth=1;

    sf_axis amx,amy,az;
    sf_axis alx,aly;
    sf_axis aw,ae;

    /* I/O files */
    sf_file Fs=NULL;    /*  slowness file S(nlx,nly,nz   ) */
    sf_file Fd=NULL;    /*      data file D(nmx,nmy,   nw) */
    sf_file Fw=NULL; 

    /* I/O slices */
    sf_fslice slow=NULL;
    sf_fslice data=NULL;
    sf_fslice wfld=NULL;


    wexcub3d cub; /* wavefield hypercube */
    wextap3d tap; /* tapering */
    wexssr3d ssr; /* SSR operator */
    wexslo3d slo; /* slowness */

    wexop3d weop;

    float dsmax;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
#endif
    
    /* default mode is migration/modeling */
    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */
    if (!sf_getfloat(  "eps",&eps   ))   eps =  0.01; /* stability parameter */
    if (!sf_getbool("causal",&causal)) causal= false; /* causality flag */
    if (!sf_getint(  "nrmax",&nrmax )) nrmax =     1; /* maximum references */
    if (!sf_getfloat("dtmax",&dtmax )) dtmax = 0.004; /* max time error */

    if (!sf_getint(    "pmx",&pmx   ))   pmx =     0; /* padding on x */
    if (!sf_getint(    "pmy",&pmy   ))   pmy =     0; /* padding on y*/

    if (!sf_getint(    "tmx",&tmx   ))   tmx =     0; /* taper on x*/
    if (!sf_getint(    "tmy",&tmy   ))   tmy =     0; /* taper on y */

    /*------------------------------------------------------------*/
    /* data files */
    Fs = sf_input ("slo");
    Fd = sf_input ( "in");
    Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
    if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");

    alx = sf_iaxa(Fs,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs,2); sf_setlabel(aly,"ly");
    az = sf_iaxa(Fs,3); sf_setlabel(az,"z" );

    amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fw,amx,1);
    amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fw,amy,2);
    ;                                           sf_oaxa(Fw,az,3);
    aw  = sf_iaxa(Fd,3); sf_setlabel(aw ,"w" ); sf_oaxa(Fw,aw ,4);
    ae  = sf_maxa(1,0,1);
   
    /*------------------------------------------------------------*/
    /* init temp files */
    slow = sf_fslice_init(sf_n(alx)*sf_n(aly)*sf_n(az),1       ,sizeof(float));
    data = sf_fslice_init(sf_n(amx)*sf_n(amy),          sf_n(aw),sizeof(sf_complex));
    wfld = sf_fslice_init(sf_n(amx)*sf_n(amy)*sf_n(az),sf_n(aw),sizeof(sf_complex));
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    sf_fslice_load(Fs,slow,SF_FLOAT);
    sf_fslice_load(Fd,data,SF_COMPLEX);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    cub = wex_cube(verb,
		   amx,amy,az,
		   alx,aly,
		   aw,
		   ae,
		   eps,
		   ompnth);
    
    dsmax = dtmax/cub->az.d;
    
    /*------------------------------------------------------------*/
    /* init structures */
    tap = wextap_init(cub->amx.n,
		      cub->amy.n,
		      SF_MIN(tmx,cub->amx.n-1), /* tmx */
		      SF_MIN(tmy,cub->amy.n-1), /* tmy */
		      true,true);
    ssr = wexssr_init(cub,pmx,pmy,tmx,tmy,dsmax);
    sf_warning("nrmax=%d",nrmax);
    slo = wexslo_init(cub,slow,nrmax,dsmax);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    wexsign=causal?+1:-1;

    weop = wex_init(cub);
    wex(weop,cub,ssr,tap,slo,wexsign,data,wfld);
    wex_close(weop);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    sf_fslice_dump(Fw,wfld,SF_COMPLEX);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* close temp files */
    sf_fslice_close(data);
    sf_fslice_close(wfld);
    sf_fslice_close(slow);    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* close structures */
    wexslo_close(slo);
    wexssr_close(cub,ssr);
    wextap_close(tap);
    /*------------------------------------------------------------*/


    exit (0);
}
