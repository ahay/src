""" "Emergency backup" replacement for the "official" rsf module.
Installed if swig or numpy/numarray not present or did not work properly.
No facilities that make use of arrays are present. Only parameter reading.
Attribute noArrays allows distinguishing between the two modules."""

import sys
from string import lower

class Par:

    def __init__(self,argv=sys.argv):
        self.noArrays = True
        self.prog = argv[0]
        self.__args = self.__argvlist2dict(argv[1:])
        print self.__args

    def __argvlist2dict(self,argv):
        """Eliminates duplicates in argv and makes it a dictionary"""
        argv = self.__filter_equal_sign(argv)
        args = {}
        for a in argv:
            key = a.split('=')[0]
            val = a.split('=')[1]
            args[key] = val
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
        if self.__args.has_key(key):
            return self.__args[key]
        elif default:
            return default
        else:
            sys.stderr.write( '\n  Argument %s= not given to %s \n' %
                              (key, self.prog))
            sys.exit(1)

    def string(self, key):
        """Returns string argument given to program"""
        return self.__get(key, default=None)

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
        val = lower(val) # No, combining with line above does not work
        if val == 'y' or val == 'true':
            return True
        elif val =='n' or val == 'false':
            return False
        else:
            sys.stderr.write( '\n  Argument %s= to %s must be bool (y/n, True/False) \n' %
                              (key, self.prog))
            sys.exit(1)