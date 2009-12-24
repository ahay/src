from rsfproj import *
import fdmod

# ------------------------------------------------------------
# model parameters
def paramwin():
    par = {
        'nx':1601, 'ox':25.000,'dx':0.025,  'lx':'x', 'ux':'km',
        'nz':601,  'oz':4.5,   'dz':0.025,  'lz':'z', 'uz':'km',
        'nt':1500, 'ot':0,     'dt':0.008,  'lt':'t', 'ut':'s'
        }
    
    par['ft2km']=0.3048
    
    par['ox']=par['ox']*par['ft2km']
    par['dx']=par['dx']*par['ft2km']
    par['oz']=par['oz']*par['ft2km']
    par['dz']=par['dz']*par['ft2km']

    par['jsnap']=300

    return par

# ------------------------------------------------------------
def getstrvelwin(velo,par):

    strvelfile = 'data/sigsbee/sigsbee2a_stratigraphy.sgy'
    strvelfile = 'sigsbee2a_stratigraphy.sgy'
    Fetch(strvelfile,'sigsbee')

    Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
         strvelfile,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)

    Flow(velo,
         velo+'-raw',
         '''
         scale rscale=0.001 |
         scale rscale=%g |
         put 
         o1=%g d1=%g label1=%s unit1=%s
         o2=%g d2=%g label2=%s unit2=%s |
         window n1=%d min1=%g n2=%d min2=%g
         ''' % (par['ft2km'],
                0.0                ,0.0250*par['ft2km'],par['lz'],par['uz'],
                10.000*par['ft2km'],0.0250*par['ft2km'],par['lx'],par['ux'],
                par['nz'],par['oz'],
                par['nx'],par['ox']
                ))
