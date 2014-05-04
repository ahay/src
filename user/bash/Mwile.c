/* Process data with GIMP 2.0.
   
Input samples must be in byte (unsigned char) format. Preprocess data
with sfbyte first, if they have float or other samples. Only first 2D
section (n1 x n2) is taken from input data, the rest is ignored.
Use sfwindow and/or sftransp to generate input in adequate order.
o1, o2, d1, and d2 are ignored; data are treated by GIMP as equally
spaced.

Be advised, that a greyscale image file is created; not every
GIMP plugin is capable of working with greyscale input.

Examples:

sfwile < in.rsf command="gimp-equalize TRUE" > out.rsf
sfwile < in.rsf command="plug-in-spread 5 5" > out.rsf
sfwile < in.rsf command="plug-in-gauss-rle 4.0 TRUE TRUE" > out.rsf
sfwile < in.rsf command="plug-in-sobel TRUE TRUE TRUE" > out.rsf
sfwile < in.rsf command="plug-in-neon 5.0 0.0" > out.rsf
sfwile < in.rsf command="plug-in-emboss 315.0 45.0 7 TRUE" > out.rsf
sfwile < in.rsf command="plug-in-ripple 27 2 1 0 0 TRUE TRUE" > out.rsf
sfwile < in.rsf command="plug-in-whirl-pinch 45 0.0 1.0" > out.rsf

Documentation for basic GIMP plugins can be found here:
http://docs.gimp.org/en/filters.html
Information about additional plugins is collected here:
http://registry.gimp.org/
*/

/*
  Copyright (C) 2008 University of Texas at Austin
  
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

/*

   This utility simple creates a temporary graphics file in TGA format,
   dumps input RSF data into it, creates a Script-Fu script (Scheme
   language), installs it with gimptool-2.0, invokes GIMP with the user
   specified command, uninstalls the script, outputs contents of the
   processed tmp TGA file as RSF data in the end.

   You can read about basic for the GIMP batch mode here:
   http://www.gimp.org/tutorials/Basic_Batch/
   TGA file format:
   http://en.wikipedia.org/wiki/Truevision_TGA

*/

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>

/* BSD - MAXNAMELEN, Posix - NAME_MAX */
#ifndef NAME_MAX
#ifdef MAXNAMELEN
#define	NAME_MAX MAXNAMELEN
#else
#ifdef FILENAME_MAX
#define NAME_MAX FILENAME_MAX
#endif
#endif
#endif

#include <rsf.h>

#pragma pack(1)

typedef struct {
    unsigned char ident_size;       /* size of ID field that follows 18 byte header (0 usually) */
    unsigned char colormap_type;    /* type of color map 0 = none, 1 = has palette */
    unsigned char image_type;       /* type of image 0 = none,1 = indexed, 2 = rgb, 3 = greyscale, + 8 = rle packed */

    unsigned short colormap_start;  /* first color map entry in palette */
    unsigned short colormap_length; /* number of colors in palette */
    unsigned char colormap_bits;    /* number of bits per palette entry 8, 16, 24, 32 */

    unsigned short x_start;         /* image x origin */
    unsigned short y_start;         /* image y origin */
    unsigned short width;           /* image width in pixels */
    unsigned short height;          /* image height in pixels */
    unsigned char bits;             /* image bits per pixel 8, 16, 24, 32 */
    unsigned char descriptor;       /* image descriptor bits (vh flip bits) */

    /* colormap (if defined) and pixel data follow header */
} TGA_HEADER;

#define BUFFER_SIZE 1024

int main (int argc, char* argv[]) {
    int n1, n2;
    unsigned char **buf;
    unsigned char *row;
    sf_file in, out = NULL;
    TGA_HEADER header;
    int i, j;
    FILE *gimp_file = NULL;
    FILE *scm_file = NULL;
    char *tmp_filename = NULL;
    char *tga_filename = NULL;
    char *scm_filename = NULL;
    char *eq, *cmd = NULL, *cmd_start = NULL;
    char buffer[BUFFER_SIZE + 1];
    int cmd_len = 0, arg_num = 0;

    sf_init (argc, argv);

    if (0 != system ("gimptool-2.0 --quiet"))
        sf_error ("Can not execute gimptool-2.0: make sure it is installed and reachable through $PATH");

    if (NULL == sf_getstring ("command")) sf_error ("Need command=");
    /* Command to be executed by GIMP */

    /* Search for command= statement among the command line arguments */
    for (i = 1; i < argc; i++) {
        eq  = strchr (argv[i], '=');
        if (NULL == eq)
            continue;
        if (0 == strncmp (argv[i], "command", 7)) { /* Found command= argument */
            cmd = eq;
            cmd++;
            j = 0;
            while (*cmd && isspace ((int) *cmd)) /* Skip leading spaces if any */
                cmd++;
            if (*cmd) /* Command name begins here */
                cmd_start = cmd;
            while (*cmd) { /* Start scanning command and its arguments */
                if (isspace ((int) *cmd)) { 
                    if (0 == cmd_len) /* First space - end of the name of the command */
                        cmd_len = j;
                    while (*cmd && isspace ((int) *cmd)) /* Skip spaces till the next argument */
                        cmd++;
                    if (*cmd) /* Found another argument after spaces */
                        arg_num++;
                } else {
                    j++;
                    cmd++;
                }
            }
            if (0 == cmd_len && cmd_start)
                cmd_len = j;
        }
    }

    if (NULL == cmd)
        sf_error ("Need command=");

    if (NULL == cmd_start)
        sf_error ("Need nonempty command=");

    cmd = (char*)sf_ucharalloc (cmd_len + 1);
    strncpy (cmd, cmd_start, cmd_len);
    cmd[cmd_len] = '\0';

    in = sf_input ("in");
    out = sf_output ("out");
    sf_settype (out, SF_UCHAR);

    if (SF_UCHAR != sf_gettype (in))
        sf_error ("Input data is not in byte format");

    if (!sf_histint (in, "n1", &n1))
        sf_error ("No n1= in input");
    if (n1 < 1)
        sf_error ("n1 can not be less than 1");
    if (!sf_histint (in, "n2", &n2))
        sf_error ("No n2= in input");
    if (n2 < 1)
        sf_error ("n2 can not be less than 1");

    if (0 != system ("gimp -i -b '(gimp-quit 0)'"))
        sf_error ("Can not execute gimp: make sure it is installed and reachable through $PATH");

    tmp_filename = sf_charalloc (NAME_MAX + 1);
    tga_filename = sf_charalloc (NAME_MAX + 1);
    scm_filename = sf_charalloc (NAME_MAX + 1);
    snprintf (tmp_filename, NAME_MAX, "/tmp/%sXXXXXX", sf_getprog ());
    gimp_file = fdopen (mkstemp (tmp_filename), "w+");
    if (NULL == gimp_file)
        sf_error ("%s: cannot create %s:", __FILE__, tmp_filename);

    snprintf (tga_filename, NAME_MAX, "%s.tga", tmp_filename);
    snprintf (scm_filename, NAME_MAX, "%s.scm", tmp_filename);

    scm_file = fopen (scm_filename, "w+");
    if (NULL == scm_file)
        sf_error ("%s: cannot create %s:", __FILE__, scm_filename);

    /* Create Scheme script for executing the command */
    fprintf (scm_file, "(define (wrapper-for-");
    fprintf (scm_file, "%s filename", cmd);
    for (i = 0; i < arg_num; i++) {
        fprintf (scm_file, " arg%d", i);
    }
    fprintf (scm_file, ")\n");
    fprintf (scm_file, "   (let* ((image (car (gimp-file-load RUN-NONINTERACTIVE filename filename)))\n");
    fprintf (scm_file, "          (drawable (car (gimp-image-get-active-layer image))))\n");
    fprintf (scm_file, "     (");

    /* Plugin or not plugin command? */
    if (0 == strncmp (cmd, "plug-in", 7))
        fprintf (scm_file, "%s RUN-NONINTERACTIVE image drawable", cmd);
    else
        fprintf (scm_file, "%s drawable", cmd);

    for (i = 0; i < arg_num; i++) {
        fprintf (scm_file, " arg%d", i);
    }
    fprintf (scm_file, ")\n");
    fprintf (scm_file, "     (gimp-file-save RUN-NONINTERACTIVE image drawable filename filename)\n");
    fprintf (scm_file, "     (gimp-image-delete image)))\n");
    fflush (scm_file);
    fclose (scm_file);

    /* Install the script */
    snprintf (buffer, BUFFER_SIZE, "gimptool-2.0 --quiet --install-script %s", scm_filename);
    if (0 != system (buffer))
        sf_error ("Can not install a Script-Fu script with gimptool-2.0: check if GIMP installed correctly");

    /* Read 2D section */
    buf = sf_ucharalloc2 (n1, n2);
    row = sf_ucharalloc (n2);
    sf_ucharread (buf[0], n1 * n2, in);

    /* Fill out the header */
    header.ident_size = 0;
    header.colormap_type = 0;
    header.image_type = 3; /* Greyscale */
    header.colormap_start = 0;
    header.colormap_length = 0;
    header.colormap_bits = 8; /* 8-bit colormap samples, for compatibility */
    header.x_start = 0;
    header.y_start = 0;
    header.width = n2;
    header.height = n1;
    header.bits = 8; /* 8-bit data samples */
    header.descriptor = 0;

    /* Write the header */
    fwrite (&header, sizeof (TGA_HEADER), 1, gimp_file);

    /* Write the data line by line */
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++)
            row[j] = buf[j][n1 - i - 1];
        fwrite (row, n2 * sizeof (unsigned char), 1, gimp_file);
    }
    fflush (gimp_file);
    fclose (gimp_file);

    /* Make .tga extension for the file */
    rename (tmp_filename, tga_filename);

    /* Execute the command through the script */
    snprintf (buffer, BUFFER_SIZE, "gimp -i -b '(wrapper-for-%s \"%s\" %s)' -b '(gimp-quit 0)'",
              cmd, tga_filename, &cmd_start[cmd_len]);

    if (0 != system (buffer)) {
        snprintf (buffer, BUFFER_SIZE, "gimptool-2.0 --quiet --uninstall-script %s", scm_filename);
        if (0 != system (buffer)) sf_warning("failed gimptool:");
        unlink (tga_filename);
        unlink (scm_filename);
        sf_error ("GIMP returned an error upon execution of the command: check if the command is specified correctly");
    }

    /* Read the resulted image */
    gimp_file = fopen (tga_filename, "r+");
    if (NULL == gimp_file) {
        snprintf (buffer, BUFFER_SIZE, "gimptool-2.0 --quiet --uninstall-script %s", scm_filename);
        if (0 != system (buffer)) sf_warning("failed gimptool:");
        unlink (tga_filename);
        unlink (scm_filename);
        sf_error ("%s: cannot open %s:", __FILE__, tga_filename);
    }

    /* Read the header */
    if (1 != fread (&header, sizeof (TGA_HEADER), 1, gimp_file))
	sf_error("fread error:");

    /* Dimensions might have changed */
    n2 = header.width;
    n1 = header.height;
    free (buf);
    free (row);
    buf = sf_ucharalloc2 (n1, n2);
    row = sf_ucharalloc (n2);

    /* Read the data line by line */
    for (i = 0; i < n1; i++) {
        if (n2 != fread (row, sizeof (unsigned char), n2, gimp_file))
	    sf_error("fread error:");
        for (j = 0; j < n2; j++)
            buf[j][n1 - i - 1] = row[j];
    }
    fclose (gimp_file);

    /* Output RSF data */
    sf_putint (out, "n1", n1);
    sf_putint (out, "n2", n2);
    sf_ucharwrite (buf[0], n1 * n2, out);

    /* Uninstall the script */
    snprintf (buffer, BUFFER_SIZE, "gimptool-2.0 --quiet --uninstall-script %s", scm_filename);
    if (0 != system (buffer)) sf_warning("error running \"%s\"",buffer);
    unlink (tga_filename);
    unlink (scm_filename);



    exit (0);
}
