/* 3D CLoud WINdowing */
#include <rsf.h>

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  bool verb;
  float apt, cosapt;       /* aperture */

  /* I/O files */
  sf_file Fin = NULL;      /* input   cloud */
  sf_file Fou = NULL;      /* output  cloud */
  sf_file Fo  = NULL;      /* orbit   */

  pt3d       * oco=NULL;   /* orbit coordinates  */
  vc3d       * ono=NULL;   /* orbit normals      */
  pt3d         gco;        /* ground coordinates */

  float      * jnk=NULL;
  int          nco=9;      /* coordinates (position, normal, velocity) */

  sf_axis ao,ag,aw;        /* cube axes */
  int     io,ig;

  float R,Rx,Ry,Rz;
  float cosobl;
  int ncount=0;

  int ipass;

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  /* default behavior */
  if (!sf_getbool ("verb",&verb)) verb=false; /* verbosity  */
  if (!sf_getfloat("apt", &apt))   apt=15.0;  /* aperture (deg) */
  cosapt = cos(apt*SF_PI/180);                /* cos(aperture) */
  if(verb) {
    sf_warning("apt=%6.2f",apt);
  }

  /* setup i/o and auxiliary file */
  Fin = sf_input ( "in");
  Fou = sf_output("out");
  Fo  = sf_input ( "oo");

  /* coordinate axes */
  ao = sf_iaxa(Fo ,2); sf_setlabel(ao,"o"); /* orbit */
  ag = sf_iaxa(Fin,2); sf_setlabel(ag,"g"); /* ground input */
  if(verb) sf_raxa(ao);
  if(verb) sf_raxa(ag);

  /*------------------------------------------------------------*/
  /* orbit coordinates */
  oco = (pt3d*) sf_alloc(sf_n(ao),sizeof(*oco));
  ono = (vc3d*) sf_alloc(sf_n(ao),sizeof(*ono));
  jnk = sf_floatalloc( nco );

  for( io = 0; io < sf_n(ao); io++) {
    sf_floatread( jnk,nco,Fo);
    oco[io].x  = jnk[0];
    oco[io].y  = jnk[1];
    oco[io].z  = jnk[2];
    ono[io].dx = jnk[3];
    ono[io].dy = jnk[4];
    ono[io].dz = jnk[5];
  }

  /*------------------------------------------------------------*/
  /* ground coordinates */
  for( ipass = 0; ipass < 2; ipass++) {
    sf_seek(Fin,0,SEEK_SET); // seek to the start of the input file

    for( ig = 0; ig < sf_n(ag); ig++) {
      sf_floatread( jnk,nco,Fin);
      gco.x  = jnk[0];
      gco.y  = jnk[1];
      gco.z  = jnk[2];

      for( io = 0; io < sf_n(ao); io++ ) {

        // compute distance
        Rx = gco.x - oco[io].x;
        Ry = gco.y - oco[io].y;
        Rz = gco.z - oco[io].z;
        R = sqrtf( Rx * Rx + Ry * Ry + Rz * Rz );
        Rx /= R;
        Ry /= R;
        Rz /= R;

        cosobl = SF_ABS( Rx*ono[io].dx + Ry*ono[io].dy + Rz*ono[io].dz );
        if( cosobl > cosapt) {
          if( ipass == 0) ncount++;                    // count points
          else            sf_floatwrite( jnk,nco,Fou); // write points
          break;
        }
      }
    }

    if(ipass == 0) { // make output axis after first pass
      aw = sf_maxa(ncount,0,1);
      sf_oaxa(Fou,aw,2);
      if(verb) sf_raxa(aw);
    }

  }
  
  /*------------------------------------------------------------*/
  /* deallocate arrays */
  free(jnk);
  free(oco);
  free(ono);

  exit (0);
}
