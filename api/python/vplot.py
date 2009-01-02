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
    def uorig(self,x,y):
        c_vplot.vp_uorig(x,y)
    def orig(self,x,y):
        c_vplot.vp_orig(x,y)
    def uclip(self,xmin,ymin,xmax,ymax):
        c_vplot.vp_uclip (xmin,ymin,xmax,ymax)
    def umove(self,x,y):
        c_vplot.vp_umove(x,y)
    def udraw(self,x,y):
        c_vplot.vp_udraw(x,y)
    def move(self,x,y):
        c_vplot.vp_move(x,y)
    def draw(self,x,y):
        c_vplot.vp_draw(x,y)
    def fat(self,f):
        c_vplot.vp_fat(f)
    def color(self,col):
        c_vplot.vp_color(col)
    def penup(self):
        c_vplot.vp_penup()
    def pendn(self,x,y):
        c_vplot.vp_pendn(x,y)
    def text (self,x,y,size,orient,string):
        c_vplot.vp_text(x,y,size,orient,string)
    def utext (self,x,y,size,orient,string):
        c_vplot.vp_utext(x,y,size,orient,string)
    def scale(self,xscale,yscale):
        c_vplot.vp_scale(xscale,yscale)
    def uarrow(self,x1,y1,x,y,r):
        c_vplot.vp_uarrow(x1,y1,x,y,r)
    def tjust(self,xjust,yjust):
        c_vplot.vp_tjust(xjust,yjust)
    def clip(self,xmin,ymin,xmax,ymax):
        c_vplot.vp_clip (xmin,ymin,xmax,ymax)
