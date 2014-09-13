#!/usr/bin/env python
##   Copyright (C) 2014 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# modified from vp_annotate by Martin Karrenbach
import sys

xtpen_message='''
Move cursor to the place, where the balloon arrow 
should point to, then click left mouse button.
Fill out the label and eventually change the defaults,
then click on CONFIRM.  Repeat for more annotations.
To create the annotated file, QUIT out of xtpen.
'''

xtpen_result='''
This is the annotated vplot figure.
You might play with vpstyle=y, if you only want to 
see the original portion.
'''

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    argc = len(sys.argv)
    prog = sys.argv.pop(0)
    
    if argc < 2:
        print '''
Usage:
%s [batch=0] [vpstyle=n] file.vpl annotated.vpl

Annotates a Vplot file using sfbox
        ''' % prog
        sys.exit(2)

    interactive = 1
    vpstyle="n"
    textfile = "text_file" 

    args = []
    files = []
    for arg in sys.argv[1:]:
        if '=' in arg:
            if arg[:5] == 'batch':
                if arg[5]=='y' or arg[5]=='1':
                    interactive = 0
                else:
                    interactive = 1
            elif arg[:4] == 'text':
                textfile = arg[4:]
            elif arg[:7] == 'vpstyle':
                vpstyle = arg[7:]
            else:
                args.append(arg)
        else:
            files.append(arg)
    args = ' '.join(args)
        
def annotate(files,interactive,vpstyle,textfile)

# copy input file

$cnt=$cnt+1;
$tempfile = "temp_vplot"."$cnt";
open(INFILE,"> $tempfile.v") || die "Could not open tempfile\n";
while(<STDIN>) {print INFILE $_ ; }
close(INFILE);
	
# run xtpen in the interactive session

if ($interactive == $true) { 
  system "xtpen message=$xtpen_message $remaining interact=$textfile boxy=y <$tempfile.v";}

# digest the text file containing labels

open(TEXT,"<$textfile") || die "Could not open $textfile\n";
while(<TEXT>){
    $cnt = $cnt + 1;
    $tempfile = "temp_vplot"."$cnt" ;
    chop($_);
    system " Box $_ out=$tempfile.v head=/dev/null ";
}
close(TEXT);

$filelist="";
for ( $i=0; $i<=$cnt; $i++){ $filelist = " $filelist"."temp_vplot"."$i".".v " ;}

# run vppen

if ($interactive == $true) { 
  system "vppen $filelist  erase=once  vpstyle=$vpstyle $remaining | xtpen message=$xtpen_result"; }

open(VPPEN,"vppen $filelist  erase=once  vpstyle=$vpstyle $remaining |")||
                                                    die "Could not run vppen\n";
# write the composite vplot file to the output

while(<VPPEN>){print STDOUT $_ ;}; close(VPPEN);


system "/bin/rm $filelist" ;

exit(0);

# clean up and emergency

abort:
&abortit($cnt);
exit(-1);

sub abortit{
local($cnt);
print STDERR "Abnormal exit caught signal\n" ;
if ( $cnt >-1 ) {
   for ( $i=0; $i<=$cnt; $i++){ 
        $file = "temp_vplot"."$i".".v " ;
	if (-e $file ) { system("/bin/rm  $file") ;}
   }
}
};
