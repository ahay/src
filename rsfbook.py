##   Copyright (C) 2004 University of Texas at Austin
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
import os, string

# The following adds all SCons SConscript API to the globals of this module.
import SCons.Script.SConscript
globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

#############################################################################
# CUSTOM BUILDERS
#############################################################################

def toc(target=None,source=None,env=None):
    "Build a table of contents from a collection of papers"
    toc = open(str(target[0]),'w')
    report = string.replace(string.upper(os.path.basename(os.getcwd())),
                            'SEP','SEP--')
    toc.write('''
    \renewcommand{\REPORT}{%s} % report
    \title{\REPORT\ --- TABLE OF CONTENTS}
    \maketitle
    %% start entries

    ''')
    return 0


    
