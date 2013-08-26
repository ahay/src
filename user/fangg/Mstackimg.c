/* Stack multi-shots images 
Takes:  file0.rsf file1.rsf file2.rsf ...
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <unistd.h>

#include <rsf.h>


int main(int argc, char* argv[]) 
{

    int ii, ix, iz, ix0, iz0, nin;
    int snx, snz, lnx, lnz;
    float sox, soz, lox, loz, sdx, sdz, ldx, ldz; 
    sf_file Fins, Fin, Fout;
    
    const char **filelist;
    float **large, **small;
    
    sf_axis lax, laz, sax, saz;
    
    sf_init(argc, argv);
    Fin = sf_input("in");
    Fout = sf_output("out");
            
    filelist = (const char**) sf_alloc ( (size_t)argc, sizeof(char*));
    
    nin=0;
    for (ii=1; ii< argc; ii++) {
	if (NULL != strchr(argv[ii],'=')) continue; 
	filelist[nin] = argv[ii];
	nin++;
    }
        
    if (0==nin) sf_error ("no input");
        
    laz = sf_iaxa(Fin, 1); lnz = sf_n(laz); ldz = sf_d(laz); loz = sf_o(laz);
    lax = sf_iaxa(Fin, 2); lnx = sf_n(lax); ldx = sf_d(lax); lox = sf_o(lax);
    
    sf_oaxa(Fout, laz, 1); 
    sf_oaxa(Fout, lax, 2);

    large = sf_floatalloc2(lnz, lnx);
    
    for (ix=0; ix<lnx; ix++)
	for (iz=0; iz<lnz; iz++)
	    large[ix][iz] = 0.0;

    for (ii=0; ii <nin; ii++) {
	Fins = sf_input(filelist[ii]);
	saz = sf_iaxa(Fins, 1); snz = sf_n(saz); sdz = sf_d(saz); soz = sf_o(saz);
	sax = sf_iaxa(Fins, 2); snx = sf_n(sax); sdx = sf_d(sax); sox = sf_o(sax);
	
	if (sox < lox || sox+(snx-1)*sdx > lox+(lnx-1)*ldx ) sf_error("ox setting error !");
	if (soz < loz || soz+(snz-1)*sdz > loz+(lnz-1)*ldz ) sf_error("oz setting error !"); 
	if (sdx != ldx || sdz != ldz ) sf_error("d1, d2 setting error !");

	small =sf_floatalloc2(snz, snx);
	sf_floatread(small[0], snx*snz, Fins);

	ix0 = (int) (sox - lox)/ldx ;
	iz0 = (int) (soz - loz)/ldz ;

	for (ix=0; ix<snx; ix++) {
	    for (iz=0; iz<snz; iz++) {
		large[ix+ix0][iz+iz0] += small[ix][iz];
	    }
	}

	free(*small);
	free(small);
	sf_fileclose(Fins);
		
    }

    for (ix=0; ix<lnx; ix++) 
    sf_floatwrite(large[ix], lnz, Fout);

    exit(0);

}
    
	
    
