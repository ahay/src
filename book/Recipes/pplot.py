try:    from rsf.cluster import *
except: from rsf.proj    import *
import fdmod 

# ------------------------------------------------------------
def multip(plot,allplots,ny,nx,ys,xs,yc,xc):

    alljnk=[]
    for iy in range(ny):
        for ix in range(nx):

            ii = iy*nx+ix

            plt = allplots[ii]
            jnk = '_' + plt

            Plot(jnk,
                 plt,
                 'Overlay',
                 vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f '
                 % (ys,xs,iy*yc,ix*xc))

            alljnk.append(jnk)

    Plot  (plot,alljnk,'Overlay')
    Result(plot,alljnk,'Overlay')

def animate(plot,allplots,ny,nx,ys,xs,yc,xc):
    Result(plot,allplots,'Movie')
# ------------------------------------------------------------

def p2x2(plot,p0,p1,p2,p3,ys,xs,yc,xc):
    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2
    j3 = '_' + p3

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0 '% (ys,xs,yc   ))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'% (ys,xs,yc,xc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=0 '% (ys,xs)      )
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f ycenter=0  xcenter=%f'% (ys,xs,   xc))

    Plot(plot,[j0,j1,j2,j3],'Overlay')

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

    Plot(plot,[j0,j1,j2,j3,j4,j5],'Overlay')

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

    Plot(plot,[j0,j1,j2,j3,j4,j5],'Overlay')
    
# ------------------------------------------------------------

def p2x1(plot,p0,p1,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,1*yc))

    Plot  (plot,[j0,j1],'Overlay')
#    Result(plot,[j0,j1],'Overlay')

def p3x1(plot,p0,p1,p2,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2

    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,1*yc))
    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,2*yc))

    Plot(plot,[j0,j1,j2],'Overlay')

def p4x1(plot,p0,p1,p2,p3,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2
    j3 = plot + '_' + p3

    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,0*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,1*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,2*yc))
    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,3*yc))

    Plot(plot,[j0,j1,j2,j3],'Overlay')

# ------------------------------------------------------------

def p1x2(plot,p0,p1,ys,xs,yc):
    j0 = '_' + p0
    j1 = '_' + p1

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))

    Plot  (plot,[j0,j1],'Overlay')
#    Result(plot,[j0,j1],'Overlay')
    
def p1x3(plot,p0,p1,p2,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,2*yc))

    Plot(plot,[j0,j1,j2],'Overlay')

def p1x4(plot,p0,p1,p2,p3,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2
    j3 = plot + '_' + p3

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,2*yc))
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,3*yc))

    Plot(plot,[j0,j1,j2,j3],'Overlay')

def p1x5(plot,p0,p1,p2,p3,p4,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2
    j3 = plot + '_' + p3
    j4 = plot + '_' + p4

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,2*yc))
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,3*yc))
    Plot(j4,p4,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,4*yc))

    Plot(plot,[j0,j1,j2,j3,j4],'Overlay')

def p1x6(plot,p0,p1,p2,p3,p4,p5,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2
    j3 = plot + '_' + p3
    j4 = plot + '_' + p4
    j5 = plot + '_' + p5

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,2*yc))
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,3*yc))
    Plot(j4,p4,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,4*yc))
    Plot(j5,p5,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,5*yc))

    Plot(plot,[j0,j1,j2,j3,j4,j5],'Overlay')

def p1x7(plot,p0,p1,p2,p3,p4,p5,p6,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2
    j3 = plot + '_' + p3
    j4 = plot + '_' + p4
    j5 = plot + '_' + p5
    j6 = plot + '_' + p6

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,2*yc))
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,3*yc))
    Plot(j4,p4,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,4*yc))
    Plot(j5,p5,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,5*yc))
    Plot(j6,p6,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,6*yc))

    Plot(plot,[j0,j1,j2,j3,j4,j5,j6],'Overlay')

def p1x8(plot,p0,p1,p2,p3,p4,p5,p6,p7,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2
    j3 = plot + '_' + p3
    j4 = plot + '_' + p4
    j5 = plot + '_' + p5
    j6 = plot + '_' + p6
    j7 = plot + '_' + p7

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,2*yc))
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,3*yc))
    Plot(j4,p4,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,4*yc))
    Plot(j5,p5,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,5*yc))
    Plot(j6,p6,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,6*yc))
    Plot(j7,p7,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,7*yc))

    Plot(plot,[j0,j1,j2,j3,j4,j5,j6,j7],'Overlay')

def p1x9(plot,p0,p1,p2,p3,p4,p5,p6,p7,p8,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1
    j2 = plot + '_' + p2
    j3 = plot + '_' + p3
    j4 = plot + '_' + p4
    j5 = plot + '_' + p5
    j6 = plot + '_' + p6
    j7 = plot + '_' + p7
    j8 = plot + '_' + p8

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,0*yc))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,1*yc))
    Plot(j2,p2,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,2*yc))
    Plot(j3,p3,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,3*yc))
    Plot(j4,p4,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,4*yc))
    Plot(j5,p5,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,5*yc))
    Plot(j6,p6,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,6*yc))
    Plot(j7,p7,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,7*yc))
    Plot(j8,p8,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=0'% (ys,xs,8*yc))

    Plot(plot,[j0,j1,j2,j3,j4,j5,j6,j7,j8],'Overlay')


