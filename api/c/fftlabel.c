/* Fix labels and units for FFT transforms */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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
#include <string.h>

#include "fftlabel.h"

#include "file.h"
/*^*/

#include "alloc.h"

void sf_fft_unit(int axis         /* axis number */,
		 const char* unit /* input unit */,
		 sf_file out      /* output file */)
/*< Change unit to 1/unit or vice versa >*/
{
    size_t len;
    char *unit2, varname[7];

    snprintf(varname,7,"unit%d",axis);
    if (NULL != unit) {
	if (0==strcmp(unit,"s")) {
	    sf_putstring(out,varname,"Hz");
	} else if (0==strcmp(unit,"Hz")) {
	    sf_putstring(out,varname,"s");
	} else if (0==strncmp(unit,"1/",2)) {
	    sf_putstring(out,varname,unit+2);
	} else {
	    len=strlen(unit)+3;
	    unit2 = sf_charalloc(len);
	    snprintf(unit2,len,"1/%s",unit);
	    sf_putstring(out,varname,unit2);
	}
    }
}

bool sf_fft_label(int axis          /* axis number */,
		  const char* label /* input label */, 
		  sf_file out       /* output file */) 
/*< Choose an output label appropriately. 
  Returns false if the label name is not recognized. >*/
{
    char varname[8];

    snprintf(varname,8,"label%d",axis);
    if (0==strcmp(label,"Time")) {
	sf_putstring(out,varname,"Frequency");
    } else if (0==strcmp(label,"time")) {
	sf_putstring(out,varname,"frequency");
    } else if (0==strcmp(label,"Frequency")) {
	sf_putstring(out,varname,"Time");
    } else if (0==strcmp(label,"frequency")) {
	sf_putstring(out,varname,"time");
    } else {
	return false;
    }
    return true;
}

