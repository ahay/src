# Copyright (C) 2004 University of Texas at Austin
#  
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import string, sys, os, re
import rsf.doc
import rsf.prog

def Flow(sources,flow,bindir,rsfflow=1,
         checkpar=False,coms=[],prefix='sf',progsuffix='',
         remote='', stdout=1,stdin=1,timer='',
         mpirun=None,workdir=None,batch=None,np=1,wall=''):
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
            if rsfflow:
                if command[:len(prefix)]==prefix:
                    rsfprog = command
                else:
                    rsfprog = prefix + command            
                if rsf.doc.progs.has_key(rsfprog):
                    if checkpar:
                        for par in pars:
                            if rsf.doc.progs[rsfprog].check(par):
                                sys.stderr.write('Failed on "%s"\n' % 
                                                 subline)
                                sys.exit(1)
                    command = os.path.join(bindir,rsfprog+progsuffix) 
                    sources.append(command)
                    if rsfprog not in coms:
                        coms.append(rsfprog)
                elif     rsfprog == prefix + 'mpi' or \
                         rsfprog == prefix + 'omp':
                    command = os.path.join(bindir,rsfprog+progsuffix) 
                    sources.append(command)
                elif command[:len(bindir)]==bindir:
                    sources.append(command)
            else:
                rsfprog = None
                if re.match(r'[^/]+\.exe$',command): # local program
                    command = os.path.join('.',command)
                elif command in os.listdir(bindir):
                    command = os.path.join(bindir,command) 
            pars.insert(0,command)
            # special rule for metaprograms
            if rsfprog and rsfprog[len(prefix):] in \
                    ('conjgrad','cconjgrad','mpi','omp'):
                n = 1
                command2 = pars.pop(n)
                while '=' in command2:
                    pars.insert(n,command2)
                    n += 1
                    command2 = pars.pop(n)
                # check if this command is in our list
                if rsfflow:
                    if command2[:len(prefix)]==prefix:
                        rsfprog2 = command2
                    else:
                        rsfprog2 = prefix + command2            
                    if rsf.doc.progs.has_key(rsfprog2):
                        command2 = os.path.join(bindir,rsfprog2+progsuffix) 
                        sources.append(command2)
                        if rsfprog2 not in coms:
                            coms.append(rsfprog2)
                if re.match(r'[^/]+\.exe$',command2): # local program
                    command2 = os.path.join('.',command2)                 
                pars.insert(n,command2)
            # special rule for MPI programs
            if rsfprog and rsfprog.startswith(prefix+'mpi') and mpirun:
                pars.insert(0,mpirun)
            #<- assemble the command line
            substep = remote + string.join(pars,' ')
            substeps.append(substep)
        #<-
        steps.append(string.join(substeps," | "))
    #<- assemble the pipeline
    command = string.join(steps," && ")
    if stdout==1:
        if rsfprog and rsfprog[:len(prefix)+3] == prefix+'mpi':
            command = command + " --output=$TARGET"
        else:
            command = command + " > $TARGET"
    elif stdout==0:
        command = command + " >/dev/null"
    if stdin:
        if rsfprog and rsfprog[:len(prefix)+3] == prefix+'mpi':
            command = command + " --input=$SOURCE"
        else:
            command = "< $SOURCE " + command
    if workdir:
        command = '/bin/rm -rf %s && /bin/mkdir %s && cd %s && ' % (workdir,workdir,workdir) + command
    command = timer + command

    if batch:
        command = '%s exe="%s" np=%d path="%s" ' % \
            (os.path.join(bindir,'sfbatch'),command,np,os.getcwd())
        if os.path.isfile(batch):
            command += ' batchfile="%s" ' % batch
        if wall:
            command += ' wall="%s" ' % wall
        
    return command
