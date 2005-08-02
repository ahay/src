from rsfproj import *

def p2x2(plot,p0,p1,p2,p3,ys,xs,yc,xc):
    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2
    j3 = '_' + p3

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0 '% (ys,xs,yc)   )
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,yc,xc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=0 '% (ys,xs)      )
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=%f'% (ys,xs,   xc))

    Result(plot,[j0,j1,j2,j3],'Overlay')
