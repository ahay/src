/*$
 
=head1 NAME

 raspen - vplot filter for Movie files (SEPlib version is preferred!)

=head1 SYNOPSIS

 raspen  -wstype=tcpr - vplot filter for Tek color plotter
 raspen  -wstype=ppm - vplot filter for ppm format output
 raspen [options] [inputfiles]

=head1 DESCRIPTION

This program creates Movie files for such things as plate tectonic movies,
or any movies based on any Vplot graphics.
This program supports EVERYTHING in Vplot to the FULL extent!

=head1 OPTIONS

B<n1=1024 n2=768 aspect=1. ppi=100. esize=1 or=n wstype=tcpr white=n>

n1 and n2 control the size of the created output raster file, with one
page for each frame in the input file.
colormult scales vplot color numbers to raster byte values. The default is
2, but for plate tectonic movies 1 might be more appropriate.
colfile=colfile gives name of color table file for color=T option of movie.
If or=y then colors are OR'd as they are drawn. IE, color 3 and color 6 will
create color 7 where they overwrite. colormult should be a power of 2
in order for this option to work. (note 1 is a power of 2)
If esize=3 then rgb triples with no color table will be written.

If wstype=ppm, then the output file is in ppm format.
The default resolution in this mode is 100dpi. If you choose
dpi=300 your output file will be big (25MB)

If wstype=tcpr, then plotting is done on the Tek color plotter.
If white=n, then the background is black instead of white.
For the Tek, the default resolution is ppi=150, half the maximum possible
resolution. The printer always stretches to fill the page, even if you ask
it not to. So just vary ppi and leave n1 and n2 alone.
Note the printer currently barfs at full resolution!

=head1 universal VPLOT filter options



in=[stdin](file)   | endpause=[y,n](d)| interact=[null](output_file)
pause=[0](secs)    | shade=[y,n]      | style=[standard,rotated]
fat=[0]            | wantras=[y,n](d) | erase=[y,n,once,literal]
fatmult=[1.0]      | overlay=[n,y]    | size=[rel,abs](d)
scale=[1.0]        | frame=[n,y]      | rotate=0 (degrees)
xscale=[1.0]       | mono=[n,y](d)    | dither=1(0=no 1=rn 2=od 3=fs 4=ht)(d)
yscale=[1.0]       | window=[y,n]     | pixc=[1.0](d)  greyc=[1.0](d)
hshift=[0.0](inch) | echo=[n,y] (d)   | xcenter,ycenter= (inch)
vshift=[0.0](inch) | break=[b,e,i]    | xwmin,ywmin,xwmax,ywmax= (inch)
txscale=[1.0]      | invras=[y,n]     | dashscale=[1.0]
txfont=[0](0-16)(d)| txsquare=[y,n]   | mkscale=[1.0]
cachepipe=[n,y]    |                  |

=head1 SEE ALSO

L<vppen>

L<pspen>

L<xtpen>

=head1 CATEGORY

B<graphics/vplot/filters>

=cut

*/

char *documentation[] = {
" ",
"NAME",
#ifdef SEP
"    Raspen - SEPlib vplot filter for Movie files",
"    Raspen wstype=tcpr - SEPlib vplot filter for Tek color plotter",
"    Raspen wstype=ppm - SEPlib vplot filter for ppm format output",
#else
"    raspen - vplot filter for Movie files (SEPlib version is preferred!)",
"    raspen wstype=tcpr - vplot filter for Tek color plotter",
"    raspen wstype=ppm - vplot filter for ppm format output",
#endif
" ",
"SYNOPSIS",
#ifdef SEP
"    Raspen [options] in=vplot-inputfile OR headerfile on standard in",
#else 
"    raspen [options] [inputfiles]",
#endif
" ",
"This program creates Movie files for such things as plate tectonic movies,",
"or any movies based on any Vplot graphics.",
"This program supports EVERYTHING in Vplot to the FULL extent!",
"",
"OPTIONS",
"n1=1024 n2=768 aspect=1. ppi=100. esize=1 or=n wstype=tcpr white=n",
"n1 and n2 control the size of the created output raster file, with one",
"page for each frame in the input file.",
"colormult scales vplot color numbers to raster byte values. The default is",
"2, but for plate tectonic movies 1 might be more appropriate.",
"colfile=colfile gives name of color table file for color=T option of movie.",
"If or=y then colors are OR'd as they are drawn. IE, color 3 and color 6 will",
"create color 7 where they overwrite. colormult should be a power of 2",
"in order for this option to work. (note 1 is a power of 2)",
"If esize=3 then rgb triples with no color table will be written.",
"",
"If wstype=ppm, then the output file is in ppm format.",
"The default resolution in this mode is 100dpi. If you choose",
"dpi=300 your output file will be big (25MB)",
"",
"If wstype=tcpr, then plotting is done on the Tek color plotter.",
"If white=n, then the background is black instead of white.",
"For the Tek, the default resolution is ppi=150, half the maximum possible",
"resolution. The printer always stretches to fill the page, even if you ask",
"it not to. So just vary ppi and leave n1 and n2 alone.",
"Note the printer currently barfs at full resolution!",
"",
#include "../include/gendoc.h"
" ",
"SEE ALSO",
"    man pen"
};
int	doclength = { sizeof documentation/sizeof documentation[0] };
