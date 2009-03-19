## 
 # BPAIT parameters
 ##

from rsfproj import *
import fdmod

# ------------------------------------------------------------
# model parameters
def param():
    par = {
    'nt':8001,  'ot':0,      'dt':0.00050, 'lt':'t', 'ut':'s',
    'nx':5395,  'ox':2.4384, 'dx':0.01250, 'lx':'x', 'ux':'km',
    'nz':1911,  'oz':0,      'dz':0.00625, 'lz':'z', 'uz':'km'
    }

    return par

# ------------------------------------------------------------
def getmigvel(velo,par):

    migvelfile = 'data/bpait/bpaitvel.hh'
    #Fetch("bpaitvel.hh",'bpait')

    Flow(velo,migvelfile,
         '''
         dd form=native |
         put
         o1=%(oz)g d1=%(dz)g label1=%(lz)s label2=%(lx)s
         o2=%(ox)g d2=%(dx)g  unit1=%(uz)s  unit2=%(ux)s |
         scale rscale=0.001
         ''' % par)
    
