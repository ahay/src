/*
point cloud spray
Paul Sava
Copyright (C) 2022 Colorado School of Mines
*/
#include <rsf.h>

/*------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  bool verb;
  bool isreal;
  int axis;

  /* I/O files */
  sf_file Fin = NULL; /* file to spray */
  sf_file Fou = NULL; /* file to spray into */
  sf_file Fcc = NULL; /* cloud   */

  sf_axis aa, ac; /* cube axes */

  float      *dR; /* data real */
  sf_complex *dC; /* data complex */

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc, argv);

  /* default behavior */
  if (!sf_getbool("verb", &verb)) verb = false; /* verbosity  */
  if (!sf_getint ("axis", &axis)) axis = 2;     /* spray axis */

  /* setup i/o and auxiliary file */
  Fin = sf_input("in");
  Fou = sf_output("out");
  Fcc = sf_input("cc");

  if (SF_FLOAT == sf_gettype(Fin)) isreal = true;
  else                             isreal = false;

  /* coordinate axes */
  aa = sf_iaxa(Fin,1);
  ac = sf_iaxa(Fcc,2);
  if(verb) {
    sf_raxa(aa);
    sf_raxa(ac);
  }

  if(axis == 1) {
    sf_oaxa(Fou, ac, 1);
    sf_oaxa(Fou, aa, 2);
  } else {
    sf_oaxa(Fou, aa, 1);
    sf_oaxa(Fou, ac, 2);
  }

  /*------------------------------------------------------------*/
  /* stop if the cloud is empty */
  if(sf_n(ac) == 0) {

    if(isreal) {

      if(axis == 1) {
        dR =   sf_floatalloc( 1 );
        for(int ic = 0; ic < 1; ic++) dR[ ic ] = 0.0;
        for(int ia = 0; ia < sf_n(aa); ia++)
          sf_floatwrite(dR, 1, Fou);
      } else {
        dR =   sf_floatalloc( sf_n(aa) );
        for(int ia = 0; ia < sf_n(aa); ia++) dR[ ia ] = 0.0;
        for(int ic = 0; ic < 1; ic++)
          sf_floatwrite(dR, sf_n(aa), Fou);
      }
      free(dR);

    } else {

      if(axis == 1) {
        dC = sf_complexalloc( 1 );
        for(int ic = 0; ic < 1; ic++) dC[ ic ] = 0.0;
        for(int ia = 0; ia < sf_n(aa); ia++) {
          sf_complexwrite(dC, 1, Fou);
        }
      } else {
        dC = sf_complexalloc( sf_n(aa) );
        for(int ia = 0; ia < sf_n(aa); ia++) dC[ ia ] = 0.0;
        for(int ic = 0; ic < 1; ic++)
          sf_complexwrite(dC, sf_n(aa), Fou);
      }
      free(dC);

    }
    exit(0);
  }

  /*------------------------------------------------------------*/

  if(axis == 1) {
    if(isreal) dR =   sf_floatalloc( sf_n(ac) );
    else       dC = sf_complexalloc( sf_n(ac) );
  } else {
    if(isreal) dR =   sf_floatalloc( sf_n(aa) );
    else       dC = sf_complexalloc( sf_n(aa) );
  }

  /*------------------------------------------------------------*/

  if(isreal) {

    if(axis == 1) {

      for(int ia = 0; ia < sf_n(aa); ia++) {
        sf_floatread (dR, 1, Fin);
        for(int ic = 1; ic < sf_n(ac); ic++)
          dR[ic] = dR[0];
        sf_floatwrite(dR, sf_n(ac), Fou);
      }

    } else {

      sf_floatread (dR, sf_n(aa), Fin);
      for(int ic = 0; ic < sf_n(ac); ic++)
        sf_floatwrite(dR, sf_n(aa), Fou);

    }

  } else {

    if(axis == 1) {

      for(int ia = 0; ia < sf_n(aa); ia++) {
        sf_complexread (dC, 1, Fin);
        for(int ic = 1; ic < sf_n(ac); ic++)
          dC[ic] = dC[0];
        sf_complexwrite(dC, sf_n(ac), Fou);
      }

    } else {

      sf_complexread (dC, sf_n(aa), Fin);
      for(int ic = 0; ic < sf_n(ac); ic++)
        sf_complexwrite(dC, sf_n(aa), Fou);

    }

  }

  /*------------------------------------------------------------*/
  if(isreal) free(dR);
  else       free(dC);

  exit(0);
}
