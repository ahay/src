/*
reformat gridded maps to point clouds
Paul Sava
Copyright (C) 2022 Colorado School of Mines
*/
#include <rsf.h>

#define NCO 6

int main(int argc, char* argv[])
{
  bool verb,sphc;
  bool mask = false;

  /* I/O files */
  sf_file Fin = NULL;      /* input   map */
  sf_file Fou = NULL;      /* output  cloud */
  sf_file Fmk = NULL;      /* cloud mask */

  sf_axis a1,aa,ac;
  int     n1;

  float * map = NULL;
  float *   x = NULL;
  float *   z = NULL;
  float * msk = NULL; /* mask array */
  int    nout;        /* count output */

  float dou[NCO];

  float deg2rad = SF_PI / 180;
  float lat;

  float ax,az;
  float nx,nz, nn; /* normal */

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  /* default behavior */
  if (!sf_getbool("verb",   &verb))     verb=false; /* verbosity  */
  if (!sf_getbool("sphc",   &sphc))     sphc=false; /* spherical coordinates */

  /* setup i/o files */
  Fin = sf_input ( "in");
  Fou = sf_output("out");

  /* coordinate axes */
  a1 = sf_iaxa(Fin,1); n1 = sf_n(a1);

  /* use absent-value mask */
  if(NULL != sf_getstring("msk")){
    mask = true;
    if(verb) sf_warning("mask=%d",mask);

    Fmk = sf_input("msk");
    msk = sf_floatalloc( n1 );
    sf_floatread(msk,n1,Fmk);

    // count present values
    nout = 0;
    for(int i1 = 0; i1 < n1; i1++) {
      if(msk[i1] != 0) nout++;
    }

    sf_fileclose(Fmk); 
  } else {
    nout = n1;
  }

  aa = sf_maxa(  NCO, 0.0, 1.0);
  ac = sf_maxa( nout, 0.0, 1.0);

  sf_oaxa(Fou,aa,1);
  sf_oaxa(Fou,ac,2);

  /*------------------------------------------------------------*/
  map = sf_floatalloc( n1 );
  sf_floatread(map,n1,Fin);

  x   = sf_floatalloc( n1 );
  z   = sf_floatalloc( n1 );

  if(sphc) {

    for(int i1 = 0; i1 < n1; i1++) {
      lat = sf_o(a1) + i1 * sf_d(a1);
      lat *= deg2rad;

      x[i1] = map[i1] * cosf(lat);
      z[i1] = map[i1] * sinf(lat);
    }

  } else {

    for(int i1 = 0; i1 < n1; i1++) {
      x[i1] = sf_o(a1) + i1 * sf_d(a1);
      z[i1] = map[i1];
    }

  }

  /*------------------------------------------------------------*/
  for(int i1 = 0; i1 < n1; i1++) {

    /* compute vector a */
    if( i1 < n1 - 1) {
      ax = x[i1 + 1] - x[i1];
      az = z[i1 + 1] - z[i1];
    } else {
      ax = x[i1] - x[i1 - 1];
      az = z[i1] - z[i1 - 1];
    }

    /* compute the normal */
    nx = -az;
    nz = +ax;
    nn = sqrtf(nx*nx + nz*nz);

    /* output cloud point */
    dou[0] = x[i1];
    dou[1] = z[i1];

    dou[2] = - nx/nn;
    dou[3] = - nz/nn;

    dou[4] = 0.0;
    dou[5] = 0.0;

    if( mask == true ) {
      if( msk[i1] != 0  ) {
        sf_floatwrite(dou, NCO, Fou);
      }
    } else {
      sf_floatwrite(dou, NCO, Fou);
    }
  }

  /*------------------------------------------------------------*/
  free(map);
  free(  x);
  free(  z);

  exit (0);
}
