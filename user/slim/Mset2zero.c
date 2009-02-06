/* replaces content of RSF file with zeros

Takes: file1.rsf [file2.rsf ...]

Used to create initial guess for SLIMpy.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <rsf.h>
#define BUFFER 4096

int main (int argc, char *argv[])
{
    int i;
    sf_file fpt;
    FILE *FPT;
    char *iname;
    sf_ulargeint ecount;
    int esize;
    sf_ulargeint eleft;
    char buffer[BUFFER];
   
    sf_init(argc,argv);

    memset(buffer,0,BUFFER);
    for (i=1; i < argc; i++) {
        fpt=sf_input(argv[i]);
        iname = sf_histstring(fpt,"in");
        if (strcmp(iname,"stdin")==0) {
            fprintf(stderr,"FATAL ERROR: cannot handle packed file %s\n",argv[i]);
            exit(1);
        }
        ecount = sf_fileulargesize(fpt);
        sf_histint(fpt,"esize",&esize);
        sf_fileclose(fpt);
        eleft=ecount*esize;
        /* printf("%s %s %lu %d %lu\n",argv[i],iname,ecount,esize,eleft); */
        FPT=fopen(iname,"w");
        rewind(FPT);
        while (eleft>0) {
            size_t bsize=eleft>BUFFER?BUFFER:eleft;
            fwrite(buffer,1,bsize,FPT);
            eleft-=bsize;
        }
        assert(eleft==0);
        fclose(FPT);
    }

    exit (0);
}
