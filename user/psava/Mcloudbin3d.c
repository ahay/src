/*
point cloud binning
Paul Sava
Copyright (C) 2022 Colorado School of Mines
*/
#include <rsf.h>

/*------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  bool verb;

  const char **fname; // array of file names
  int nCLOUD;         //     number of point clouds

  /* I/O files */
  sf_file Fin = NULL; /*  grid file  */
  sf_file Fou = NULL; /*   bin file  */
  sf_file Fcc = NULL; /* cloud files */

  sf_axis alon,alat, ac;  /* cube axes */
  int     ilon,ilat, ic;

  float ** obin;
  float    lat; /* latitude */
  float    lon; /* longitude */

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc, argv);

  /* default behavior */
  if (!sf_getbool("verb", &verb)) verb = false; /* verbosity  */

  // array of file names
  fname = (const char**) sf_alloc( (size_t) argc, sizeof(char*) );

  // count the clouds
  nCLOUD = 0;
  for(int i = 1; i < argc; i++) {
    if (NULL != strchr(argv[i],'=')) {
      continue; /* not a file */
    } else {
      fname[nCLOUD] = argv[i];
      nCLOUD++;
    }
  }
  if( nCLOUD == 0 ) sf_error("no input");

  // list the clouds
  if(verb)
    for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++)
      sf_warning(">>> %4d %s",jCLOUD,fname[jCLOUD] );

  /*------------------------------------------------------------*/
  /* setup bin file */
  Fin = sf_input("in");
  Fou = sf_output("out");

  /* coordinate axes */
  alat = sf_iaxa(Fin,1);
  alon = sf_iaxa(Fin,2);
  if(verb) {
    sf_raxa(alat);
    sf_raxa(alon);
  }

  /*------------------------------------------------------------*/
  // allocate bin file
  obin = sf_floatalloc2(sf_n(alat),sf_n(alon));
  for( ilon = 0; ilon < sf_n(alon); ilon++) {
    for( ilat = 0; ilat < sf_n(alat); ilat++) {
        obin[ilon][ilat] = 1.0; // avoid division by 0 in bin normalization
    }
  }

  // loop over clouds
  for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++) {
    Fcc = sf_input( fname[jCLOUD] );
    ac = sf_iaxa(Fcc,2);

    for( ic = 0; ic < sf_n(ac); ic++) {
      // read cloud coordinates
      sf_floatread ( &lon, 1, Fcc);
      sf_floatread ( &lat, 1, Fcc);

      ilat = floor( (lat - sf_o(alat)) / sf_d(alat) );
      ilon = floor( (lon - sf_o(alon)) / sf_d(alon) );
      //sf_warning("%d |  %g %d      %g %d",ic, lon,ilon, lat,ilat);
      if( ilat >= 0 && ilat < sf_n(alat) &&
          ilon >= 0 && ilon < sf_n(alon)) {
        obin[ilon][ilat] += 1;
      }
    }

    sf_fileclose(Fcc);
  }

  sf_floatwrite(obin[0],sf_n(alat)*sf_n(alon),Fou);
  free(obin);
  /*------------------------------------------------------------*/

  free(fname);

  exit(0);
}
