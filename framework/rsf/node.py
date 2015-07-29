# Copyright (C) 2013 University of Texas at Austin
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


# next steps:
# does lock block second task?
# add change back to not running at end of task
# write code in myproj to build hosts.txt
# run test on stampede

import os, sys, random

def cpus(): 
    '''
    Returns the number of CPUs in the system
    '''
    try:
        # Tested on CentOS, Fedora and Cygwin 1.7
        import multiprocessing # module added in Python 2.6
        return multiprocessing.cpu_count()
    except: 
        # Thanks to Lawrence Oluyede on python-list 
        num = 0
        
        if sys.platform == 'win32':
            try:
                num = int(os.environ['NUMBER_OF_PROCESSORS'])
            except (ValueError, KeyError):
                pass
        elif sys.platform == 'darwin':
            try:
                num = int(os.popen('sysctl -n hw.ncpu').read())
            except ValueError:
                pass
        else:
            # A better way: parse /proc/cpuinfo for physical CPUs
            # rather than virtual CPUs from hyperthreading
            
            try:
                num = os.sysconf('SC_NPROCESSORS_ONLN')
            except (ValueError, OSError, AttributeError):
                pass
            
        if num >= 1:
            return num
        else:
            raise NotImplementedError

class Hosts(object):
    def __init__(self):
        'initialize list of nodes'
        
        # RSF_CLUSTER method
        #print "get RSF_CLUSTER in node.py"
        cluster = os.environ.get('RSF_CLUSTER')
        #print "test cluster=",cluster
        if not cluster:
            #print "get SLURM_NODELIST"
            test=os.environ.get('SLURM_NODELIST')
            #print "test=", test
            if not test:
                #print "set RSF_CLUSTER from SLURM_NODELIST"
                os.system("export RSF_CLUSTER=`slurm_nodelist2rsf_cluster `")
                cluster = os.environ.get('RSF_CLUSTER')
                #print "new cluster=",cluster
        #print "test2 cluster"
        if not cluster:        
            cluster = 'localhost %d' % cpus()
        #print "cluster=",cluster
        hosts = cluster.split()
        self.nodes = []
        self.local = True
        for i in range(1,len(hosts),2):
            nh = int(hosts[i])
            host = hosts[i-1]
            self.nodes.extend([host]*nh)
            if self.local and host != 'localhost':
                self.local = False

        # randomly shuffle the list of nodes
        random.shuffle(self.nodes)

    def local(self):
        'determine if everything is local'
        return self.local

    def action(self,target=None,source=None,env=None):
        'action for SCons'

        if not self.nodes: # empty list
            sys.stderr.write('No available nodes.\n')
            return -1
            
        host = self.nodes.pop()

        command = env.get('command')

        if host != 'localhost':
            command = 'ssh %s \"%s\" ' % (host,command)

        print command
        retcode = os.system(command)
        
        self.nodes.append(host)
        return retcode

