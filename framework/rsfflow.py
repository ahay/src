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
import string, sys, os, re

import rsfdoc

top = os.environ.get('RSFROOT')
bindir = os.path.join(top,'bin')

def Flow(sources,flow,rsf=1,checkpar=False,coms=[],prefix='sf',progsuffix='',remote='',stdout=1,stdin=1,timer=''):
    'Output a command line'
    lines = string.split(str(flow),'&&')
    steps = []
    for line in lines:
        substeps = []
        sublines = string.split(line,'|')
        for subline in sublines:           
            pars = string.split(subline)
            # command is assumed to be always first in line
            command = pars.pop(0)
            # check if this command is in our list
            if rsf:
                if command[:2]==prefix:
                    # assuming prefix is two chars (sf)
                    rsfprog = command
                else:
                    rsfprog = prefix + command            
                if rsfdoc.progs.has_key(rsfprog):
                    if checkpar:
                        for par in pars:
                            if rsfdoc.progs[rsfprog].check(par):
                                sys.stderr.write('Failed on "%s"\n' % 
                                                 subline)
                                sys.exit(1)
                    command = os.path.join(bindir,rsfprog+progsuffix) 
                    sources.append(command)
                    if rsfprog not in coms:
                        coms.append(rsfprog)
            else:
                rsfprog = None
                if re.match(r'[^/]+\.exe$',command): # local program
                    command = os.path.join('.',command)
            pars.insert(0,command)
            # special rule for solvers
            if rsfprog == prefix+'conjgrad' or \
                   rsfprog == prefix+'cconjgrad':
                command = pars.pop(1)
                # check if this command is in our list
                if rsf:
                    if command[:2]==prefix:
                        # assuming prefix is two chars
                        rsfprog = command
                    else:
                        rsfprog = prefix + command            
                    if rsfdoc.progs.has_key(rsfprog):
                        command = os.path.join(bindir,rsfprog+progsuffix) 
                        sources.append(command)
                        if rsfprog not in coms:
                            coms.append(rsfprog)
                else:
                    rsfprog = None
                    if re.match(r'[^/]+\.exe$',command): # local program
                        command = os.path.join('.',command)                 
                pars.insert(1,command)
            #<- assemble the command line
            substep = remote + string.join(pars,' ')
            substeps.append(substep)
        #<-
        steps.append(string.join(substeps," | "))
    #<- assemble the pipeline
    command = string.join(steps," && ")
    if stdout==1:
        command = command + " > $TARGET"
    elif stdout==0:
        command = command + " >/dev/null"
    if stdin:
        command = "< $SOURCE " + command
    command = timer + command
        
    return command
