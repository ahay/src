"""
Backup replacement for the 'official' rsf module.
To be used if swig or numpy/numarray not present or did not work properly.
No facilities that make use of arrays are present. Only parameter reading.
Attribute noArrays allows distinguishing between the two modules.
"""
# Copyright (C) 2007 Ioan Vlad
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

import sys
from string import lower

class Par:

    def __init__(self,argv=sys.argv):
        self.noArrays = True
        self.prog = argv[0]
        self.__args = self.__argvlist2dict(argv[1:])

    def __argvlist2dict(self,argv):
        """Eliminates duplicates in argv and makes it a dictionary"""
        argv = self.__filter_equal_sign(argv)
        args = {}
        for a in argv:
            key = a.split('=')[0]
            args[key] = a.replace(key+'=','')
        return args

    def __filter_equal_sign(self,argv):
        """Eliminates "par = val", "par= val" and "par =val" mistakes."""
        argv2 = []
        # Could not use the simpler 'for elem in argv'...argv.remove because
        # lonely '=' signs are treated weirdly. Many things did not work as
        # expected -- hence long and ugly code. Test everything.
        for i in range( len(argv) ):
            if argv[i] != '=':
                if argv[i].find('=') != 0:
                    if argv[i].find('=') != -1:
                        if argv[i].find('=') != len(argv[i])-1:
                            argv2.append(argv[i])
        return argv2

    def __get(self, key, default):
        """Obtains value of argument from dictionary"""
        if key in self.__args:
            return self.__args[key]
        elif str(default):
            return default
        else:
            sys.stderr.write( '\n  Argument %s= not given to %s \n' %
                              (key, self.prog))
            sys.exit(1)

    def string(self, key, default=None):
        """Returns string argument given to program"""
        return self.__get(key, default)

    def int(self,key,default=None):
        """Returns integer argument given to program"""
        try:
            return int( self.__get(key, default) )
        except:
            sys.stderr.write( '\n  Argument %s= to %s must be integer\n' %
                              (key, self.prog))
            sys.exit(1)

    def float(self,key,default=None):
        """Returns float argument given to program"""
        try:
            return float( self.__get(key, default) )
        except:
            sys.stderr.write( '\n  Argument %s= to %s must be float\n' %
                              (key, self.prog))
            sys.exit(1)

    def bool(self,key,default=None):
        """Returns bool argument given to program"""
        val = self.__get(key, default)
        val = lower(str(val)) # No, combining with line above does not work
        if val == 'y' or val == 'true':
            return True
        elif val =='n' or val == 'false':
            return False
        else:
            msg = ('\n  Argument %s= to %s must be bool (y/n, True/False) \n' %
                   (key, self.prog))
            sys.stderr.write(msg)
            sys.exit(1)

################################################################################

class Input:

    def __create_variable_dictionary(self, header):
        'Parse RSF header into a dictionary of variables'
        self.vd={} # variable dictionary
        ilist = header.split()
        pos = 0
        squot = "'"
        dquot = '"'
        while pos < len(ilist):
            if '=' in ilist[pos]:
                tokenlist = ilist[pos].split('=')
                lhs = tokenlist[0]
                rhs = tokenlist[1]
                quotmark = None
                if rhs[0] in (squot, dquot):
                    if rhs[0] == squot:
                        quotmark = squot
                    else:
                        quotmark = dquot
                    if rhs[-1] == quotmark:
                        rhs_out = rhs.strip(quotmark)
                        pos += 1
                    else:
                        rhs_out = rhs.lstrip(quotmark)
                        while pos < len(ilist):
                            pos += 1
                            rhs_out += ' '
                            if ilist[pos][-1] == quotmark:
                                rhs_out += ilist[pos][:-1]
                                break
                            else:
                                rhs_out += ilist[pos]
                else:
                    rhs_out = rhs
                    pos += 1
                self.vd[lhs] = rhs_out
            else:
                pos += 1

    def __init__(self):
        # Temporary solution. Need to scan for \EOL\EOL\EOT, else this will
        # choke on a .HH file!
        # Also, need to add capability of reading from another file than
        # input
        self.__create_variable_dictionary(sys.stdin.read())

    def int(self, nm):
        return int(self.vd[nm])

    def float(self, nm):
        return float(self.vd[nm])
