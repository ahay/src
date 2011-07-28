/* Display the size of RSF files.

Takes: file1.rsf file2.rsf ...

Prints the element size, number of elements, and number of bytes
for a list of RSF files.  Non-RSF files are ignored.
*/
/*
# Copyright (C) 2009 James W. Jennings Jr.

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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <rsf.h>

int main (int argc, char* argv[])
{
    char            *filename=NULL;
    const char      *prefix=" KMGT";
    int             i, j, esize, st;
    off_t           bytes, total, p1024[5];
    float           size_human;
    bool            files, human, su;
    struct stat     buf;
    sf_file         file=NULL;

    sf_init (argc,argv);

    if (!sf_getbool("files",&files)) files=true;
    /* If y, print size of each file.  If n, print only total. */

    if (!sf_getbool("human",&human)) human=false;
    /* If y, print human-readable file size.  If n, print byte count. */

    if (!sf_getbool("su",&su)) su=false;
    /* Same for Seismic Unix */

    /* Set up powers of 1024 */
    p1024[0] = 1;
    for (i=1; i<5; i++) p1024[i] = p1024[i-1]*1024;

    total = 0;

    sf_file_error(false);

    for (i=1; i<argc; i++)      /* Loop over argument list      */
    {
        filename = argv[i];
                                /* Skip if parameter            */
        if (NULL != strchr(filename,'=')) continue;

                                /* Skip if not rsf file         */

	if (su) {
	    if (strcmp(&filename[strlen(filename)-3],".su") != 0) continue;

	    st = stat(filename,&buf);
	    if (0 != st) sf_error ("%s: cannot find file size:",__FILE__);

	    bytes = buf.st_size;
	} else {
	    if (strcmp(&filename[strlen(filename)-4],".rsf") != 0) continue;

	    /* Get file properties          */
	    file = sf_input(filename);
	    if (NULL == file) continue;

	    bytes = sf_bytes(file);
	    if (bytes < 0) { /* reading from "stdin" */
		bytes = sf_filesize (file);
		if (!sf_histint(file,"esize",&esize) || esize <= 0)
		    sf_error("%s: need esize > 0 in input",__FILE__);
		bytes *= esize;
	    }
	    sf_fileclose(file);
	}

        total += bytes;

        if (files)              /* Print file size              */
        {
            if (human)          /* ... in human readable format */
            {
                for (j=0; j<5; j++) if (p1024[j] > bytes) break;
                size_human = ((float) bytes)/((float) p1024[j-1]);
                printf ("%s: \t %6.1f %cB\n",
                        filename,size_human,prefix[j-1]);
            }
                                /* ... or simple byte count     */
#if defined(__cplusplus) || defined(c_plusplus)
            else printf ("%s: \t %10ld\n",filename,(long) bytes);
#else
	    else printf ("%s: \t %10lld\n",filename,(long long) bytes);
#endif
        }
    }

                                /* Print total size             */
    if (human)                  /* ... in human readable format */
    {
        for (j=0; j<5; j++) if (p1024[j] > total) break;
        size_human = ((float) total)/((float) p1024[j-1]);
        printf ("size= %6.1f %cB\n",size_human,prefix[j-1]);
    }
                                /* ... or simple byte count     */
#if defined(__cplusplus) || defined(c_plusplus)
    else printf ("size= %10ld\n",(long) total);
#else
    else printf ("size= %10lld\n",(long long) total);
#endif

    exit (0);
}

/* 	$Id$	 */
