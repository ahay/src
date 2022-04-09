/* 3D CLoud WINdowing */
#include <rsf.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

#define NCO 9

bool isInCone(pt3d *oco, vc3d *ono, pt3d *gco, float cosapt)
{
  float Rx,Ry,Rz,R, cosobl;

  // cone axis vector
  Rx = gco->x - oco->x;
  Ry = gco->y - oco->y;
  Rz = gco->z - oco->z;
  R = sqrtf( Rx * Rx + Ry * Ry + Rz * Rz );
  Rx /= R;
  Ry /= R;
  Rz /= R;

  cosobl = SF_ABS( Rx*ono->dx + Ry*ono->dy + Rz*ono->dz );
  if( cosobl > cosapt) return true;
  else                 return false;
}


/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  bool verb, fast;
  float apt, cosapt;       /* aperture */

  /* I/O files */
  sf_file Fin = NULL;      /* input   cloud */
  sf_file Fou = NULL;      /* output  cloud */
  sf_file Fo  = NULL;      /* orbit   */

  sf_axis ao,ag,aw;        /* cube axes */
  size_t  no,ng,nw;
  size_t  io,ig,iw;
  size_t ncount = 0;
  int ico;

  pt3d       * oco=NULL;   /* orbit coordinates  */
  vc3d       * ono=NULL;   /* orbit normals      */
  pt3d         gco;        /* ground coordinate  */
  float      * jnk=NULL;
  float      * din=NULL;
  float      * dou=NULL;

  bool * fin;

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  /* OMP init */
  #ifdef _OPENMP
    omp_init();
  #endif

  /* default behavior */
  if (!sf_getbool ("verb",&verb)) verb=false; /* verbosity  */
  if (!sf_getbool ("fast",&fast)) fast=true;  /* in-core windowing  */
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

  no = sf_n(ao);
  ng = sf_n(ag);

  /*------------------------------------------------------------*/
  /* orbit coordinates */
  oco = (pt3d*) sf_alloc(no,sizeof(*oco));
  ono = (vc3d*) sf_alloc(no,sizeof(*ono));
  jnk = sf_floatalloc( NCO );

  for( io = 0; io < no; io++) {
    sf_floatread( jnk,NCO,Fo);
    oco[io].x  = jnk[0];
    oco[io].y  = jnk[1];
    oco[io].z  = jnk[2];
    ono[io].dx = jnk[3];
    ono[io].dy = jnk[4];
    ono[io].dz = jnk[5];
  }

  /* ground coordinates */
  if(fast) {
    /*------------------------------------------------------------*/
    /* FAST */
    /*------------------------------------------------------------*/

    // read input
    if(verb) sf_warning("read points");
    din = sf_floatalloc( ng * NCO );
    sf_floatread(din, ng * NCO, Fin);

    // flag window points
    if(verb) sf_warning("flag window points");
    fin = sf_boolalloc ( ng );
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) \
  private( ig, io, gco ) \
  shared(  ng, no, din, fin, oco,ono,cosapt )
#endif
    for( ig = 0; ig < ng; ig++) {
      fin[ig] = 0;

      gco.x  = din[ig * NCO];
      gco.y  = din[ig * NCO + 1];
      gco.z  = din[ig * NCO + 2];

      for( io = 0; io < no; io++ ) {
        if( isInCone( &oco[io], &ono[io], &gco, cosapt) ) {
          fin[ig] = 1;
          break;
        }
      }
    }

    // count window points
    if(verb) sf_warning("count window points");
    nw = 0;
#ifdef _OPENMP
  #pragma omp parallel for reduction(+ : nw)
#endif
    for( ig = 0; ig < ng; ig++) nw += fin[ig];

    // avoid empty window
    nw = SF_MAX(nw,1);

    // write window header
    if(verb) sf_warning("write window header");
    aw = sf_maxa(nw,0,1);
    sf_oaxa(Fou,aw,2);
    if(verb) sf_raxa(aw);

    // index window points
    if(verb) sf_warning("index window points");
    dou = sf_floatalloc( nw * NCO );
    iw = 0;
    for( ig = 0; ig < ng; ig++) {
      if( fin[ig] ) {
        dou[ iw * NCO ] = ig;
        iw++;
      }
    }

    // move window points
    if(verb) sf_warning("move window points");
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) \
  private( iw, ig, ico) \
  shared(  nw, dou )
#endif
    for( iw = 0; iw < nw; iw++ ) {
      ig = dou[iw * NCO];
      for(ico = 0; ico < NCO; ico++) {
        dou[iw * NCO + ico] = din[ig * NCO + ico];
      }
    }

    // write window points
    if(verb) sf_warning("write window points");
    sf_floatwrite(dou, nw * NCO, Fou);

    // deallocate arrays
    free(din);
    free(dou);

  } else {
    /*------------------------------------------------------------*/
    /* SLOW */
    /*------------------------------------------------------------*/

    for( int ipass = 0; ipass < 2; ipass++) {
      sf_seek(Fin,0,SEEK_SET); // seek to the start of the input file

      for( ig = 0; ig < sf_n(ag); ig++) {

        sf_floatread( jnk,NCO,Fin);
        gco.x  = jnk[0];
        gco.y  = jnk[1];
        gco.z  = jnk[2];

        for( io = 0; io < sf_n(ao); io++ ) {
          if( isInCone( &oco[io], &ono[io], &gco, cosapt) ) {
            if( ipass == 0) ncount++;         // count points
            else sf_floatwrite( jnk,NCO,Fou); // write points
            break;
          }
        }

      }

      if(ipass == 0) { // make output axis after first pass
        aw = sf_maxa( SF_MAX(ncount,1),0,1);
        sf_oaxa(Fou,aw,2);
        if(verb) sf_raxa(aw);
      }

    } // ipass

  }

  /*------------------------------------------------------------*/
  /* deallocate arrays */
  free(jnk);
  free(oco);
  free(ono);

  exit (0);
}
