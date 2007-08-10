from rsfproj import *

def p2x2(plot,p0,p1,p2,p3,ys,xs,yc,xc):
    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2
    j3 = '_' + p3

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0 '% (ys,xs,yc   ))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,yc,xc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=0 '% (ys,xs)      )
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=%f'% (ys,xs,   xc))

    Result(plot,[j0,j1,j2,j3],'Overlay')

def p2x3(plot,p0,p1,p2,p3,p4,p5,ys,xs,yc,xc):
    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2
    j3 = '_' + p3
    j4 = '_' + p4
    j5 = '_' + p5

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0 '% (ys,xs,yc     ))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,yc,  xc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,yc,2*xc))
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=0 '% (ys,xs        ))
    Plot(j4,p4,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=%f'% (ys,xs,     xc))
    Plot(j5,p5,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=%f'% (ys,xs,   2*xc))

    Result(plot,[j0,j1,j2,j3,j4,j5],'Overlay')

def p3x2(plot,p0,p1,p2,p3,p4,p5,ys,xs,yc,xc):
    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2
    j3 = '_' + p3
    j4 = '_' + p4
    j5 = '_' + p5

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,   0-1,0-1))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,  yc-1,0-1))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,2*yc-1,0-1))

    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,   0-1,xc-1))
    Plot(j4,p4,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,  yc-1,xc-1))
    Plot(j5,p5,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,2*yc-1,xc-1))

    Result(plot,[j0,j1,j2,j3,j4,j5],'Overlay')
    
def p2x1(plot,p0,p1,ys,xs,yc):
    j0 = '_' + p0
    j1 = '_' + p1

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=0'% (ys,xs   ))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,yc))

    Plot  (plot,[j0,j1],'Overlay')
    Result(plot,[j0,j1],'Overlay')

def p3x1(plot,p0,p1,p2,ys,xs,yc):
    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=0'% (ys,xs     ))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,  yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,2*yc))

    Result(plot,[j0,j1,j2],'Overlay')

def p1x2(plot,p0,p1,ys,xs,yc):
    j0 = '_' + p0
    j1 = '_' + p1

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=0 xcenter=-1 '% (ys,xs   ))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=0 xcenter=%f'% (ys,xs,yc))

    Plot  (plot,[j0,j1],'Overlay')
    Result(plot,[j0,j1],'Overlay')
    
def p4x1(plot,p0,p1,p2,p3,ys,xs,yc):
    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2
    j3 = '_' + p3

    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=0'% (ys,xs     ))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,  yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,2*yc))
    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,3*yc))

    Result(plot,[j0,j1,j2,j3],'Overlay')
