/* Display the size of RSF files.

Takes: file1.rsf file2.rsf ...

Prints the element size, number of elements, and number of bytes
for a list of RSF files.  Non-RSF files are ignored.
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

/* 	$Id$	 */

#include <string.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    char            *filename, *prefix=" KMGT";
    int             i, j, esize, dim;
    long long       size, bytes, total, p1024[5];
    float           size_human;
    bool            files, human;
    sf_ulargeint    n[SF_MAX_DIM];
    sf_file         file;

    sf_init (argc,argv);

    if (!sf_getbool("files",&files)) files=true;
    /* If y, print size of each file.  If n, print only total. */

    if (!sf_getbool("human",&human)) human=false;
    /* If y, print human-readable file size.  If n, print byte count. */
    
    /* Set up powers of 1024 */
    p1024[0] = 1;
    for (i=1; i<5; i++) p1024[i] = p1024[i-1]*1024;

    total = 0;

    for (i=1; i<argc; i++)      /* Loop over argument list      */
    {
        filename = argv[i];
                                /* Skip if parameter            */
        if (NULL != strchr(filename,'=')) continue;
        
                                /* Skip if not rsf file         */
        if (strcmp(&filename[strlen(filename)-4],".rsf") != 0) continue;
        
                                /* Get file properties          */
        file = sf_input(filename);
        if (sf_histint(file,"esize",&esize));
        else esize = sf_esize(file);
        dim = sf_fileulargedims(file,n);
        sf_fileclose(file);
            
        size = 1;               /* Calculate file size          */
        for (j=0; j < dim; j++) size *= n[j];    
        bytes  = size*esize;
        total += bytes;
    
        if (files)              /* Print file size              */
        {
            if (human)          /* ... in human readable format */
            {
                for (j=0; j<5; j++) if (p1024[j] > bytes) break;
                size_human = ((float) bytes)/((float) p1024[j-1]);
                printf ("%d  %10lld  %6.1f %cB  %s\n",
                        esize,size,size_human,prefix[j-1],filename);
            }
                                /* ... or simple byte count     */
            else printf ("%d  %10lld  %10lld  %s\n",esize,size,bytes,filename);
        }
    }

                                /* Print total size             */
    if (human)                  /* ... in human readable format */
    {
        for (j=0; j<5; j++) if (p1024[j] > total) break;
        size_human = ((float) total)/((float) p1024[j-1]);
        printf ("               %6.1f %cB  total\n",size_human,prefix[j-1]);
    }
                                /* ... or simple byte count     */
    else printf ("               %10lld  total\n",total);

    exit (0);
}
