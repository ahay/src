import os, sys

import c_vplot

class Vplot(object):
    def __init__(self,name=None):
        if name:
            # redirect stdout
            vfile = os.open(name, os.O_WRONLY | os.O_CREAT)
            os.dup2(vfile,sys.stdout.fileno())
            os.close(vfile)
        c_vplot.vp_init()
