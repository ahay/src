from rsf.proj import *
import fdmod

# ------------------------------------------------------------
# model parameters
def param():
    par = {
        'nx':3201,  'ox':10.025,'dx':0.025,  'lx':'x',    'ux':'km',
        'nz':1201,  'oz':0,     'dz':0.025,  'lz':'z',    'uz':'km',
        'ft2km':0.3048
        }

    par['ft2km']=0.3048
    
    par['ox']=par['ox']*par['ft2km']
    par['dx']=par['dx']*par['ft2km']
    par['oz']=par['oz']*par['ft2km']
    par['dz']=par['dz']*par['ft2km']

    return par

# ------------------------------------------------------------
# modeling parameters
def fdmpar(par):

    par['nt']=5000
    par['ot']=0
    par['dt']=0.001
    par['lt']='Time'
    par['ut']='s'

    par['lo']=0.25*par['ft2km']           # wavelength
    par['fo']=4.92*par['ft2km']/par['lo'] # central frequency in water
#    print  par['lo'], par['dx'], par['fo']

    par['kt']=100  # wavelet delay
    par['tpad']=100
    par['tcut']=50

    par['nb']=100
    par['jsnap']=250    

    par['vbias']=4.75*par['ft2km']
    par['wweight']=10

    par['mintplot']=1.5

    return par

# ------------------------------------------------------------
def rtmpar(par):

    par['nhz']=3
    par['nhx']=3
    par['nht']=10

    return par

# ------------------------------------------------------------
def velocity(vel,par):
    
    vstr = 'sigsbee2a_stratigraphy.sgy'
    Fetch(vstr,'sigsbee')
    
    Flow('zvstr tzvstr ./shead ./bshead',vstr,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)
    
    o1=0
    d1=0.025
    o2=10.025
    d2=0.025

    o1=o1*par['ft2km']
    d1=d1*par['ft2km']
    o2=o2*par['ft2km']
    d2=d2*par['ft2km']
    
    Flow(vel,'zvstr',
         '''
         put o1=%g d1=%g o2=%g d2=%g label1=z unit1=m label2=x unit2=m |
         scale rscale=%g
         ''' %(o1,d1,o2,d2,par['ft2km']/1000.) )

# ------------------------------------------------------------
def reflectivity(ref,par):

    vstr = 'data/sigsbee/sigsbee2a_reflection_coefficients.sgy'
#    Fetch(refl,'sigsbee')
    
    Flow('zrefl tzrefl ./rhead ./brhead',vstr,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)

    o1=0
    d1=0.025
    o2=10.025
    d2=0.025

    o1=o1*par['ft2km']
    d1=d1*par['ft2km']
    o2=o2*par['ft2km']
    d2=d2*par['ft2km']
     
    Flow(ref,'zrefl',
         '''
         put o1=%g d1=%g o2=%g d2=%g label1=z unit1=m label2=x unit2=m |
         scale rscale=%g
         ''' %(o1,d1,o2,d2,par['ft2km']/1000.) )

# ------------------------------------------------------------
def win1(vou,vin,par):

    par['nx']=1600
    par['nz']=1000
    
    par['ox']=25*par['ft2km']
    par['oz']=05*par['ft2km']

    Flow(vou,vin,'window n1=%d min1=%g n2=%d min2=%g'
         % (par['nz'],par['oz'],par['nx'],par['ox']) )

# ------------------------------------------------------------
def win2(vou,vin,par):

    par['nx']=2000
    par['nz']=1000
    
    par['ox']=11*par['ft2km']
    par['oz']=05*par['ft2km']

    Flow(vou,vin,'window n1=%d min1=%g n2=%d min2=%g'
         % (par['nz'],par['oz'],par['nx'],par['ox']) )

# ------------------------------------------------------------
def win3(vou,vin,par):

#    par['nx']=1400
#    par['nz']=1000
    par['nx']=1400
    par['nz']=1000
    
    par['ox']=30*par['ft2km']
    par['oz']=5*par['ft2km']

    Flow(vou,vin,'window n1=%d min1=%g n2=%d min2=%g'
         % (par['nz'],par['oz'],par['nx'],par['ox']) )



# ------------------------------------------------------------
def saltmask(msk,vel,par):

    Flow(msk,vel,
         '''
         mask min=%g |
         dd type=float
         ''' %(14.76*par['ft2km']) )

# ------------------------------------------------------------
def watermask(msk,vel,par):

    Flow(msk,vel,
         '''
         mask max=%g |
         dd type=float
         ''' %(4.921*par['ft2km']))

# ------------------------------------------------------------
def sedimentmask(msk,vel,par):

    Flow(msk,vel,
         '''
         mask min=%g max=%g |
         dd type=float
         ''' %(4.921*par['ft2km'],
               14.76*par['ft2km']) )

# ------------------------------------------------------------
def randomedge(noise,vel,par):

    saltmask(vel+'-salt',vel,par)
    
    Flow(vel+'-smooth',vel,
         '''
         window n2=1 |
         smooth rect1=100 |
         spray axis=2 n=%(nx)d o=%(ox)g d=%(dx)g
         ''' %par)
    
    Flow(vel+'-sediments',[vel,vel+'-salt',vel+'-smooth'],
         '''
         math m=${SOURCES[1]} s=${SOURCES[2]}
         output="input*(1-m)+s*m"
         ''')
    
    # ------------------------------------------------------------
    Flow(  vel+'-taper',vel+'-salt','smooth rect1=20 rect2=20 repeat=1')
    
    Flow(vel+'-saltbd',vel+'-taper','mask min=0.10 max=0.90 | dd type=float')
    Flow(vel+'-saltin',vel+'-taper','mask min=0.91          | dd type=float')
    
    Flow(noise,[vel,vel+'-saltbd',vel+'-saltin',vel+'-taper'],
         '''
         math output=1 |
         noise |
         smooth rect1=5 rect2=5 repeat=3 |
         mask min=1.00 |
         dd type=float |
         add mode=p ${SOURCES[1]} |
         add ${SOURCES[2]}
         ''')

# ------------------------------------------------------------

    
def xplot(x,y,custom,par):
    j1=12; j2=12
    npts=par['nz']*par['nx']/j1/j2
    Flow(x+y,[x,y],
         '''cat axis=3 ${SOURCES[1]} |
            transp plane=23 |
            transp plane=12 |
            window j2=%d j3=%d |
            put n1=2 n2=%d n3=1 '''%(j1,j2,npts))
    Plot(x+y,'''
    dd type=complex |
    graph symbol=. plotcol=5 plotfat=5 label1=%s label2=%s unit1= unit2= %s
    '''%(x,y,custom) )
