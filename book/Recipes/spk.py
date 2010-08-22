from rsf.proj import *

# ------------------------------------------------------------
def delt(spk,m,n,par):

    par['nzspk']=m
    par['nxspk']=n

    par['spikex']=par['nxspk']/2+1
    par['spikez']=par['nzspk']/2+1
    par['middlex']=(par['nxspk']/2)*par['dx']
    par['middlez']=(par['nzspk']/2)*par['dz']

    Flow(spk,None,
         '''
         spike nsp=1 mag=1
         n1=%(nzspk)d o1=-%(middlez)g d1=%(dz)g k1=%(spikez)g
         n2=%(nxspk)d o2=-%(middlex)g d2=%(dx)g k2=%(spikex)g |
         put
         label1=%(lz)s label2=%(lx)s
         unit1=%(uz)s unit2=%(ux)s
         ''' % par)

#    Flow(spk,None,
#         '''
#         spike nsp=1 mag=1
#         n1=%(nzspk)d o1=0 d1=1 k1=%(spikez)g
#         n2=%(nxspk)d o2=0 d2=1 k2=%(spikex)g |
#         put
#         label1=%(lz)s label2=%(lx)s
#         unit1=%(uz)s unit2=%(ux)s
#         ''' % par)




# ------------------------------------------------------------
def delt3(spk,l,m,n,par):

    par['nzspk']=l
    par['nxspk']=m
    par['nyspk']=n

    par['spikez']=par['nzspk']/2+1
    par['spikex']=par['nxspk']/2+1
    par['spikey']=par['nyspk']/2+1
#    par['middlez']=(par['nzspk']/2)*par['dz']
#    par['middlex']=(par['nxspk']/2)*par['dx']
#    par['middley']=(par['nyspk']/2)*par['dy']
    par['middlez']=0
    par['middlex']=0
    par['middley']=0
    Flow(spk,None,
         '''
         spike nsp=1 mag=1
         n1=%(nzspk)d o1=-%(middlez)g d1=%(dz)g k1=%(spikez)g
         n2=%(nxspk)d o2=-%(middlex)g d2=%(dx)g k2=%(spikex)g 
         n3=%(nyspk)d o3=-%(middley)g d3=%(dy)g k3=%(spikey)g |
         put
         label1=%(lz)s unit1=%(uz)s
         label2=%(lx)s unit2=%(ux)s
         label3=%(ly)s unit3=%(uy)s
         ''' % par)   

def deltx3(spk,l,m,n,par):

    par['nzspk']=l
    par['nxspk']=m
    par['nyspk']=n

    par['spikez']=par['nzspk']/2+1
    par['spikex']=par['nxspk']/2+1
    par['spikey']=par['nyspk']/2+1

    par['middlez']=0
    par['middlex']=0
    par['middley']=0
    Flow(spk,None,
         '''
         spike nsp=1 mag=1
         n1=%(nzspk)d k1=%(spikez)g d1=1
         n2=%(nxspk)d k2=%(spikex)g d2=1
         n3=%(nyspk)d k3=%(spikey)g d3=1|
         put
         label1=%(lz)s unit1=%(uz)s
         label2=%(lx)s unit2=%(ux)s
         label3=%(ly)s unit3=%(uy)s
         ''' % par)   



# ------------------------------------------------------------
def gaus(spk,par):

    par['spikex']=par['nx']/2
    par['spikez']=par['nz']/2
    
    par['middlex']=par['spikex']*par['dx']
    par['middlez']=par['spikez']*par['dz']
    
    Flow(spk,None,
         '''
         math output="exp(-(
         (x1-%(middlez)g)*(x1-%(middlez)g) / ( 2*%(sigz)g*%(sigz)g ) +
         (x2-%(middlex)g)*(x2-%(middlex)g) / ( 2*%(sigx)g*%(sigx)g )))"
         n1=%(nz)d o1=%(oz)g d1=%(dz)g 
         n2=%(nx)d o2=%(ox)g d2=%(dx)g |
         put label1=z unit1=km label2=x unit2=km
         ''' % par)

# ------------------------------------------------------------
def lapl(spk,par):

    par['spikex']=par['nx']/2
    par['spikez']=par['nz']/2

    Flow(spk+'-del',None,
         '''
         spike nsp=1 mag=1
         n1=%(nz)d o1=%(oz)g d1=%(dz)g k1=%(spikez)g
         n2=%(nx)d o2=%(ox)g d2=%(dx)g k2=%(spikex)g |
         put label1=%(lz)s label2=%(lx)s unit1=%(uz)s unit2=%(ux)s
         ''' % par)

    Flow(spk+'-lap',None,
         '''
         spike nsp=9 mag=0,-1,0,-1,4.000,-1,0,-1,0
         n1=3 o1=0 d1=1 k1=1,1,1,2,2,2,3,3,3
         n2=3 o2=0 d2=1 k2=1,2,3,1,2,3,1,2,3 |
         pad beg1=3 end1=3 beg2=3 end2=3
         ''')

    Flow('hlx hlx-lag',spk+'-lap',
         'afac lag=${TARGETS[1]} verb=y')

    Flow(spk,[spk+'-del','hlx'],lapldec('${SOURCES[1]}'))

# ------------------------------------------------------------
# deconvolve with laplacian
#def lapldec(flt):
#    return '''
#    helicon filt=%s div=y |
#    helicon filt=%s div=y adj=y
#    ''' %(flt,flt)

# convolve with Laplacian
#def laplcon(flt):
#    return '''
#    helicon filt=%s div=n adj=y |
#    helicon filt=%s div=n
#    ''' %(flt,flt)

# ------------------------------------------------------------
def lapldec(flt):
    return '''
    window squeeze=n
    '''
def laplcon(flt):
    return '''
    window squeeze=n
    '''
