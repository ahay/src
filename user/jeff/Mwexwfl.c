/* 3-D wavefield extrapolation with extended SSF */

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
#include "omputil.h"
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
    int  datum;           /* datuming flag */
    int    inv;           /* down/upward flag */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    int wexsign;

    sf_axis amx,amy,az;
    sf_axis alx,aly;
    sf_axis aw,ae;

    /* I/O files */
    sf_file Fs=NULL;    /*  slowness file S(nlx,nly,nz   ) */
    sf_file Fd=NULL;    /*      data file D(nmx,nmy,   nw) */
    sf_file Fw=NULL; 
    sf_file Ft=NULL;
    sf_file Fr=NULL;

    int ompnth=1;

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
    
   // ompnth=6;
    
    /* default mode is migration/modeling */
    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */
    if (!sf_getfloat(  "eps",&eps   ))   eps =  0.01; /* stability parameter */
    if (!sf_getint(  "datum",&datum )) datum =     0; /* datuming flag */
    if (!sf_getint(    "inv",&inv ))     inv =     1; /* down/upward flag */
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
    
   // printf("slowfile %d\n",Fs[1][1][1]);
    
    Fd = sf_input ( "in");
    Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
    if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");

    alx = sf_iaxa(Fs,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs,2); sf_setlabel(aly,"ly");
    az = sf_iaxa(Fs,3); sf_setlabel(az,"z" );

    amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx");
    amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); 
    aw  = sf_iaxa(Fd,3); sf_setlabel(aw ,"w" ); 
   
    if(datum)
      { ae = sf_iaxa(Fd,4); sf_setlabel(ae ,"e" ); }
    else 
      ae = sf_maxa(1,0,1);
   
    Ft = sf_tmpfile(NULL); sf_settype(Ft,SF_COMPLEX);
    Fr = sf_tmpfile(NULL); sf_settype(Fr,SF_COMPLEX);
    /*------------------------------------------------------------*/
    /* output wavefield files */
    sf_oaxa(Fw,amx,1);
    sf_oaxa(Fw,amy,2);
    if(datum){
      sf_oaxa(Fw,aw,3);
      sf_oaxa(Fw,ae,4);
    }
    else{
      sf_oaxa(Fw,az,3);
      sf_oaxa(Fw,aw,4);
    }
    
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
                      1,
		      SF_MIN(tmx,cub->amx.n-1), /* tmx */
		      SF_MIN(tmy,cub->amy.n-1), /* tmy */
                      0,
		      true,true,false);
    slo = wexslo_init(cub,Fs,nrmax,dsmax);
    ssr = wexssr_init(cub,slo,pmx,pmy,tmx,tmy,dsmax);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    wexsign=causal?+1:-1;

    weop = wex_init(cub);

    if(datum){
      sf_filecopy(Ft,Fd,SF_COMPLEX);
      sf_seek(Ft,(off_t)0,SEEK_SET);
      wexdatum(weop,cub,ssr,tap,slo,wexsign,Ft,Fr,inv);
    }
    else
      wex(weop,cub,ssr,tap,slo,wexsign,Fd,Fr,true);

    sf_filefresh(Fr);
    sf_filecopy(Fw,Fr,SF_COMPLEX);

    /*------------------------------------------------------------*/
    /* close structures */
    wexslo_close(slo);
    wex_close(weop);
    wexssr_close(cub,ssr);
    wextap2D_close(tap);
    /*------------------------------------------------------------*/

    if (Fd!=NULL) sf_fileclose(Fd);
    if (Fs!=NULL) sf_fileclose(Fs);
    if (Fw!=NULL) sf_fileclose(Fw);
    if (Ft!=NULL) sf_tmpfileclose(Ft);
    if (Fr!=NULL) sf_tmpfileclose(Fr);

    exit (0);
}
