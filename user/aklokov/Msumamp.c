/* Stack energy between two input horizons */

/*
  Copyright (C) 2013 University of Texas at Austin
  
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

int main (int argc, char* argv[]) 
{

    int n1, n2, n3;
	int n2top, n2bot, n3top, n3bot;
    float d1, o1;
    float badSample;
    int size;
    float *panel, *trace;
    float top, bottom;
    int pind, i3, i2, ind;
    float t, stack;
    sf_file dataFile, hTopFile=NULL, hBotFile=NULL, stackFile;

// Initialize RSF 
    sf_init (argc,argv);
// Input files
    dataFile = sf_input ("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dataFile) ) sf_error ("Need float input: data volume");
    /* data volume */

    if ( NULL != sf_getstring ("top") ) {
	/* top horizon */
		hTopFile = sf_input ("top");
    } else {
		sf_error ("Need float input: top horizon");
	}

    if ( NULL != sf_getstring ("bot") ) {
	/* bottom horizon */
		hBotFile = sf_input ("bot");
    } else {
		sf_error ("Need float input: bottom horizon");
	}
// Output file
    stackFile = sf_output ("out");

    // data parameters
    if ( !sf_histint   (dataFile, "n1", &n1) )  sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile, "d1", &d1) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile, "o1", &o1) )  sf_error ("Need o1= in input");
    if ( !sf_histint   (dataFile, "n2", &n2) )  sf_error ("Need n2= in input");
    if ( !sf_histint   (dataFile, "n3", &n3) )  sf_error ("Need n3= in input");

    if ( !sf_histint   (hTopFile, "n2", &n2top) )  sf_error ("Need n2= in top horizon file");
    if ( !sf_histint   (hTopFile, "n3", &n3top) )  sf_error ("Need n3= in top horizon file");	

    if ( !sf_histint   (hTopFile, "n2", &n2bot) )  sf_error ("Need n2= in bottom horizon file");
    if ( !sf_histint   (hTopFile, "n3", &n3bot) )  sf_error ("Need n3= in bottom horizon file");	

	// files-consistency checking
	if (n2bot != n2 || n3bot != n3 || n2top != n2 || n3top != n3)
		sf_error ("Horizons are not consistent with data");

    sf_putint (stackFile, "n1", 1);

    // run parameters
    if (!sf_getfloat ("badSample", &badSample)) badSample = 99999.f;
    /* non-interpreted sample in the horizons */

    size = n2 * n3;
    panel = sf_floatalloc (size);	
    trace = sf_floatalloc (n1);	

    pind = 0;
    for (i3 = 0; i3 < n3; ++i3) {
		for (i2 = 0; i2 < n2; ++i2, ++pind) {
		    panel [pind] = 0.f;

		    sf_floatread (&top, 1, hTopFile);
		    sf_floatread (&bottom, 1, hBotFile);
		    sf_floatread (trace, n1, dataFile);

		    if (fabs (top - badSample) < 1.f) continue;    // non-interpreted sample in the horizons
		    if (fabs (bottom - badSample) < 1.f) continue;			

		    t = o1;
		    ind = 0;
		    while (t < top && ind < n1) { t += d1; ++ind; }
		
		    stack = 0.f;
		    while (t < bottom && ind < n1) {
				stack += pow (trace [ind], 2);
				++ind;
				t += d1;
	    	}

	    	panel [pind] = stack;
		}
    }

    sf_floatwrite (panel, size, stackFile);

    sf_fileclose (dataFile);
    sf_fileclose (hTopFile);
    sf_fileclose (hBotFile);
    sf_fileclose (stackFile);

    free (panel);
    free (trace);

    return 0;
}
