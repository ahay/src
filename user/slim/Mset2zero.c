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
#define DBG 0

int main (int argc, char *argv[])
{
    int i;
    sf_file fpt;
    char *iname;
    off_t ecount, eleft;
    int esize;
    char buffer[BUFFER];

    sf_init(argc,argv);

    memset(buffer,0,BUFFER);
    for (i=1; i < argc; i++)
    {
        char *fname = strdup(argv[i]);
        if (NULL != strchr(fname,'=')) continue;
        fpt=sf_input(fname);
        iname = sf_histstring(fpt,"in");
        ecount = sf_filesize(fpt);
        sf_histint(fpt,"esize",&esize);
        eleft=ecount*esize;
#if defined(__cplusplus) || defined(c_plusplus)
        if (DBG) sf_warning("%s %s %lu %d %lu\n",argv[i],iname,ecount,esize,(long) eleft);
#else
	if (DBG) sf_warning("%s %s %lu %d %llu\n",argv[i],iname,ecount,esize,(long long) eleft);
#endif
        if (strcmp(iname,"stdin")==0)
        {
            sf_file FPT;
            char *tfile = strcat(strdup(fname),"~");
            FPT=sf_output(tfile);
            sf_fileflush(FPT,fpt);
            while (eleft>0)
            {
                size_t bsize=eleft>BUFFER?BUFFER:eleft;
                sf_charwrite(buffer,BUFFER,FPT);
                eleft-=bsize;
            }
            sf_fileclose(fpt);
            sf_fileclose(FPT);
            FPT=sf_input(tfile);
            fpt=sf_output(fname);
            sf_cp(FPT,fpt);
            sf_fileclose(fpt);
            sf_fileclose(FPT);
            sf_rm(tfile,false,false,false);
        }
        else
        {
            FILE *FPT;
            sf_fileclose(fpt);
            FPT=fopen(iname,"w");
            rewind(FPT);
            while (eleft>0)
            {
                size_t bsize=eleft>BUFFER?BUFFER:eleft;
                fwrite(buffer,1,bsize,FPT);
                eleft-=bsize;
            }
            fclose(FPT);
        }
        assert(eleft==0);
    }

    exit (0);
}
