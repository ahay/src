/*
in-memory transpose 12
Paul Sava
Copyright (C) 2022 Colorado School of Mines
*/
#include <rsf.h>

/*------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  bool verb, isreal;

  /* I/O files */
  sf_file Fin = NULL; /* input */
  sf_file Fou = NULL; /* output */

  sf_axis a1, a2; /* cube axes */
  size_t n1,n2,n12;
  
  float      *dR, *wR; /* data real */
  sf_complex *dC, *wC; /* data complex */

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc, argv);

  /* default behavior */
  if (!sf_getbool("verb", &verb)) verb = false; /* verbosity  */

  /* setup i/o and auxiliary file */
  Fin = sf_input("in");
  Fou = sf_output("out");

  if (SF_FLOAT == sf_gettype(Fin)) isreal = true;
  else                             isreal = false;

  /* coordinate axes */
  a1 = sf_iaxa(Fin,1); n1 = sf_n(a1);
  a2 = sf_iaxa(Fin,2); n2 = sf_n(a2);
  if (verb) {
    sf_raxa(a1);
    sf_raxa(a2);
  }
  sf_oaxa(Fou, a2, 1);
  sf_oaxa(Fou, a1, 2);

  n12 = n1 * n2;

  /*------------------------------------------------------------*/

  if(isreal) dR =   sf_floatalloc( n12 );
  else       dC = sf_complexalloc( n12 );

  if(isreal) wR =   sf_floatalloc(  n2 );
  else       wC = sf_complexalloc(  n2 );

  /*------------------------------------------------------------*/

  if(isreal) {

    sf_floatread (dR, n12, Fin);

    for(size_t i1 = 0; i1 < n1; i1++) {
      for(size_t i2 = 0; i2 < n2; i2++) {
        wR[i2] = dR[i2 * n1 + i1];
      }
      sf_floatwrite(wR, n2, Fou);
    }

  } else {

    sf_complexread (dC, n12, Fin);

    for(size_t i1 = 0; i1 < n1; i1++) {
      for(size_t i2 = 0; i2 < n2; i2++) {
        wC[i2] = dC[i2 * n1 + i1];
      }
      sf_complexwrite(wC, n2, Fou);
    }

  }

  /*------------------------------------------------------------*/
  if(isreal) { free(dR); free(wR); }
  else       { free(dC); free(wC); }

  exit(0);
}
