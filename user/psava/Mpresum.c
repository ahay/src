/* presum traces */
#include <rsf.h>

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  bool verb, isreal;
  int  nsum, left;

  /* I/O files */
  sf_file Fin = NULL;      /* input   cloud */
  sf_file Fou = NULL;      /* output  cloud */

  sf_axis a1,a2,b2;        /* cube axes */

  float      ** dinR = NULL;
  float       * douR = NULL;
  sf_complex ** dinC = NULL;
  sf_complex  * douC = NULL;

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  /* default behavior */
  if (!sf_getbool("verb",&verb)) verb=false; /* verbosity  */
  if (!sf_getint ("nsum",&nsum)) nsum=1;   /* number of traces to sum */

  /* setup i/o and auxiliary file */
  Fin = sf_input ( "in");
  Fou = sf_output("out");

  if (SF_FLOAT == sf_gettype(Fin)) isreal = true;
  else                             isreal = false;

  /* coordinate axes */
  a1 = sf_iaxa(Fin,1);
  a2 = sf_iaxa(Fin,2);
  if(verb) sf_raxa(a1);
  if(verb) sf_raxa(a2);

  b2 = sf_maxa( sf_n(a2)/nsum , sf_o(a2) + nsum*sf_d(a2)/2, sf_d(a2)*nsum );
  sf_oaxa(Fou, a1, 1);
  sf_oaxa(Fou, b2, 2);

  /*------------------------------------------------------------*/
  /* allocate arrays */
  if(isreal) dinR = sf_floatalloc2  ( sf_n(a1),nsum );
  else       dinC = sf_complexalloc2( sf_n(a1),nsum );

  if(isreal) douR = sf_floatalloc   ( sf_n(a1) );
  else       douC = sf_complexalloc ( sf_n(a1) );

  if(isreal) {
    for (left = sf_leftsize(Fin,1); left > 0; left -= nsum) {
      //nsum = SF_MIN(left,nsum);
      if(left < nsum) break;

      sf_floatread(dinR[0], sf_n(a1) * nsum , Fin);

      for(int i1 = 0; i1 < sf_n(a1); i1++) {
        douR[i1] = 0.0;
        for(int isum = 0; isum < nsum; isum++) {
          douR[i1] += dinR[isum][i1];
        }
        douR[i1] /= nsum;
      }

      sf_floatwrite(douR  , sf_n(a1)        , Fou);
    }
  } else {
    for (left = sf_leftsize(Fin,1); left > 0; left -= nsum) {
      //nsum = SF_MIN(left,nsum);
      if(left < nsum) break;

      sf_complexread(dinC[0], sf_n(a1) * nsum , Fin);

      for(int i1 = 0; i1 < sf_n(a1); i1++) {
        douC[i1] = 0.0;
        for(int isum = 0; isum < nsum; isum++) {
          douC[i1] += dinC[isum][i1];
        }
        douC[i1] /= nsum;
      }

      sf_complexwrite(douC  , sf_n(a1)        , Fou);
    }
  }

  /*------------------------------------------------------------*/
  /* deallocate arrays */
  if(isreal) { free(dinR); free(douR); }
  else       { free(dinC); free(douC); }

  exit (0);
}
