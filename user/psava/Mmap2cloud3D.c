/*
reformat gridded maps to point clouds
Paul Sava
Copyright (C) 2022 Colorado School of Mines
*/
#include <rsf.h>

#define NCO 9

int main(int argc, char* argv[])
{
  bool verb,sphc;
  bool mask = false;

  /* I/O files */
  sf_file Fin = NULL;      /* input   map */
  sf_file Fou = NULL;      /* output  cloud */
  sf_file Fmk = NULL;      /* cloud mask */

  sf_axis a1,a2,aa,ac;
  int     n1,n2;

  float ** map = NULL;
  float **   x = NULL;
  float **   y = NULL;
  float **   z = NULL;
  float ** msk = NULL; /* mask array */
  int     nout;        /* count output */

  float dou[NCO];

  float deg2rad = SF_PI / 180;
  float lat, lon;

  float ax,ay,az;
  float bx,by,bz;
  float nx,ny,nz, nn; /* normal */

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
  a2 = sf_iaxa(Fin,2); n2 = sf_n(a2);

  /* use absent-value mask */
  if(NULL != sf_getstring("msk")){
    mask = true;
    if(verb) sf_warning("mask=%d",mask);

    Fmk = sf_input("msk");
    msk = sf_floatalloc2( n1, n2 );
    sf_floatread(msk[0],n1*n2,Fmk);

    // count present values
    nout = 0;
    for(int i2 = 0; i2 < n2; i2++) {
      for(int i1 = 0; i1 < n1; i1++) {
        if(msk[i2][i1] != 0) nout++;
      }
    }

    sf_fileclose(Fmk); 
  } else {
    nout = n1 * n2;
  }

  aa = sf_maxa(  NCO, 0.0, 1.0);
  ac = sf_maxa( nout, 0.0, 1.0);

  sf_oaxa(Fou,aa,1);
  sf_oaxa(Fou,ac,2);

  /*------------------------------------------------------------*/
  map = sf_floatalloc2( n1, n2 );
  sf_floatread(map[0],n1*n2,Fin);

  x   = sf_floatalloc2( n1, n2 );
  y   = sf_floatalloc2( n1, n2 );
  z   = sf_floatalloc2( n1, n2 );

  if(sphc) {

    for(int i2 = 0; i2 < n2; i2++) {
      lon = sf_o(a2) + i2 * sf_d(a2);
      lon *= deg2rad;

      for(int i1 = 0; i1 < n1; i1++) {
        lat = sf_o(a1) + i1 * sf_d(a1);
        lat *= deg2rad;

        x[i2][i1] = map[i2][i1] * cosf(lat) * cosf(lon);
        y[i2][i1] = map[i2][i1] * cosf(lat) * sinf(lon);
        z[i2][i1] = map[i2][i1] * sinf(lat);
      }
    }

  } else {

  for(int i2 = 0; i2 < n2; i2++) {
      for(int i1 = 0; i1 < n1; i1++) {
        x[i2][i1] = sf_o(a2) + i2 * sf_d(a2);
        y[i2][i1] = sf_o(a1) + i1 * sf_d(a1);
        z[i2][i1] = map[i2][i1];
      }
    }

  }

  /*------------------------------------------------------------*/
  for(int i2 = 0; i2 < n2; i2++) {
    for(int i1 = 0; i1 < n1; i1++) {

      /* compute vector a */
      if(i2 < n2 - 1) {
        ax = x[i2 + 1][i1] - x[i2][i1];
        ay = y[i2 + 1][i1] - y[i2][i1];
        az = z[i2 + 1][i1] - z[i2][i1];
      } else {
        ax = x[i2  ][i1] - x[i2 - 1][i1];
        ay = y[i2  ][i1] - y[i2 - 1][i1];
        az = z[i2  ][i1] - z[i2 - 1][i1];
      }

      /* compute vector b */
      if(i1 < n1 - 1) {
        bx = x[i2][i1 + 1] - x[i2][i1];
        by = y[i2][i1 + 1] - y[i2][i1];
        bz = z[i2][i1 + 1] - z[i2][i1];
      } else {
        bx = x[i2][i1  ] - x[i2][i1 - 1];
        by = y[i2][i1  ] - y[i2][i1 - 1];
        bz = z[i2][i1  ] - z[i2][i1 - 1];
      }

      /* compute the normal */
      nx = ay*bz - by*az;
      ny = az*bx - bz*ax;
      nz = ax*by - bx*ay;
      nn = sqrtf(nx*nx + ny*ny + nz*nz);

      /* output cloud point */
      dou[0] = x[i2][i1];
      dou[1] = y[i2][i1];
      dou[2] = z[i2][i1];

      dou[3] = - nx/nn;
      dou[4] = - ny/nn;
      dou[5] = - nz/nn;

      dou[6] = 0.0;
      dou[7] = 0.0;
      dou[8] = 0.0;

      if( mask == true ) {
        if( msk[i2][i1] != 0 ) {
          sf_floatwrite(dou, NCO, Fou);
        }
      } else {
        sf_floatwrite(dou, NCO, Fou);
      }

    }
  }

  /*------------------------------------------------------------*/
  free(map);
  free(  *x); free(x);
  free(  *y); free(y);
  free(  *z); free(z);

  exit (0);
}
